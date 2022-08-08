// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Linear
 * \brief Dumux preconditioners for iterative solvers
 */
#ifndef DUMUX_LINEAR_PRECONDITIONERS_HH
#define DUMUX_LINEAR_PRECONDITIONERS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/indices.hh>
#include <dune/common/version.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/gsetc.hh>

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

#include <dumux/assembly/coloring.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/istlsolverregistry.hh>
#include <dumux/parallel/parallel_for.hh>

#include "parallellinearoperator.hh"

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A preconditioner based on the Uzawa algorithm for saddle-point problems of the form
 * \f$
 \begin{pmatrix}
    A & B \\
    C & D
 \end{pmatrix}

 \begin{pmatrix}
    u\\
    p
 \end{pmatrix}

 =

 \begin{pmatrix}
    f\\
    g
 \end{pmatrix}
  * \f$
 *
 * This preconditioner is especially suited for solving the incompressible (Navier-)Stokes equations.
 * Here, \f$D = 0\f$ and \f$B = C^T\f$ if \f$\rho = 1\f$.
 * We do not expect good convergence if energy or mass transport is considered.
 *
 * See: Benzi, M., Golub, G. H., & Liesen, J. (2005). Numerical solution of saddle point problems. Acta numerica, 14, 1-137 \cite benzi2005 and <BR>
 *      Ho, N., Olson, S. D., & Walker, H. F. (2017). Accelerating the Uzawa algorithm. SIAM Journal on Scientific Computing, 39(5), S461-S476 \cite ho2017
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level (for compatibility reasons, unused).
 */
template<class M, class X, class Y, int l = 1>
class SeqUzawa : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "SeqUzawa expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l== 1, "SeqUzawa expects a block level of 1.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<A, U, U>;
    using Smoother = Dune::SeqSSOR<A, U, U>;
    using AMGSolverForA = Dune::Amg::AMG<LinearOperator, U, Smoother, Comm>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    /*!
     * \brief Constructor
     *
     * \param mat The matrix to operate on.
     * \param params Collection of parameters.
     */
#if DUNE_VERSION_GTE(DUNE_ISTL,2,8)
    SeqUzawa(const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op, const Dune::ParameterTree& params)
    : matrix_(op->getmat())
#else
    SeqUzawa(const M& mat, const Dune::ParameterTree& params)
    : matrix_(mat)
#endif
    , numIterations_(params.get<std::size_t>("iterations"))
    , relaxationFactor_(params.get<scalar_field_type>("relaxation"))
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    , useDirectVelocitySolverForA_(getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.DirectSolverForA", false))
    {
        const bool determineRelaxationFactor = getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.DetermineRelaxationFactor", true);

        // AMG is needed for determination of omega
        if (determineRelaxationFactor || !useDirectVelocitySolverForA_)
            initAMG_(params);

        if (useDirectVelocitySolverForA_)
            initUMFPack_();

        if (determineRelaxationFactor)
            relaxationFactor_ = estimateOmega_();
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    virtual void pre(X& x, Y& b) {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     */
    virtual void apply(X& update, const Y& currentDefect)
    {
        using namespace Dune::Indices;

        auto& A = matrix_[_0][_0];
        auto& B = matrix_[_0][_1];
        auto& C = matrix_[_1][_0];
        auto& D = matrix_[_1][_1];

        const auto& f = currentDefect[_0];
        const auto& g = currentDefect[_1];
        auto& u = update[_0];
        auto& p = update[_1];

        // incorporate Dirichlet cell values
        // TODO: pass Dirichlet constraint handler from outside
        for (std::size_t i = 0; i < D.N(); ++i)
        {
            const auto& block = D[i][i];
            for (auto rowIt = block.begin(); rowIt != block.end(); ++rowIt)
                if (Dune::FloatCmp::eq<scalar_field_type>(rowIt->one_norm(), 1.0))
                    p[i][rowIt.index()] = g[i][rowIt.index()];
        }

        // the actual Uzawa iteration
        for (std::size_t k = 0; k < numIterations_; ++k)
        {
            // u_k+1 = u_k + Q_A^−1*(f − (A*u_k + B*p_k)),
            auto uRhs = f;
            A.mmv(u, uRhs);
            B.mmv(p, uRhs);
            auto uIncrement = u;
            applySolverForA_(uIncrement, uRhs);
            u += uIncrement;

            // p_k+1 = p_k + omega*(g - C*u_k+1 - D*p_k)
            auto pIncrement = g;
            C.mmv(u, pIncrement);
            D.mmv(p, pIncrement);
            pIncrement *= relaxationFactor_;
            p += pIncrement;

            if (verbosity_ > 1)
            {
                std::cout << "Uzawa iteration " << k
                << ", residual: " << uRhs.two_norm() + pIncrement.two_norm()/relaxationFactor_ << std::endl;
            }
        }
    }

    /*!
     * \brief Clean up.
     */
    virtual void post(X& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:

    void initAMG_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;
        auto linearOperator = std::make_shared<LinearOperator>(matrix_[_0][_0]);
        amgSolverForA_ = std::make_unique<AMGSolverForA>(linearOperator, params);
    }

    void initUMFPack_()
    {
#if HAVE_UMFPACK
            using namespace Dune::Indices;
            umfPackSolverForA_ = std::make_unique<Dune::UMFPack<A>>(matrix_[_0][_0]);
#else
            DUNE_THROW(Dune::InvalidStateException, "UMFPack not available. Use LinearSolver.Preconditioner.DirectVelocitySolver = false.");
#endif
    }

    /*!
     * \brief Estimate the relaxation factor omega
     *
     * The optimal relaxation factor is omega = 2/(lambdaMin + lambdaMax), where lambdaMin and lambdaMax are the smallest and largest
     * eigenvalues of the Schur complement -C*Ainv*B (assuming D = 0).
     * lambdaMax can be easily determined using the power iteration algorithm (https://en.wikipedia.org/wiki/Power_iteration) and lambdaMin
     * could be estimated in a similar manner. We do not consider lambdaMin because for certain cases, e.g., when C contains some rows of zeroes only,
     * this estimate will fail.
     *
     * Instead we assume that lambdaMin is sufficiently close to lambdaMax such that omega = 1/lambdaMax.
     * This seems to work rather well for various applications.
     * We will underestimate omega by a factor of 2 in the worst case (i.e, lambdaMin = 0).
     *
     * When facing convergence issues, you may set LinearSolver.Preconditioner.Verbosity = 1 to see the estimate of lambdaMax.
     * In a new simulation run, you can then set LinearSolver.Preconditioner.DetermineRelaxationFactor = false and set some other value
     * for LinearSolver.Preconditioner.Relaxation based on the estimate of lambdaMax.
     *
     * See: Benzi, M., Golub, G. H., & Liesen, J. (2005). Numerical solution of saddle point problems. Acta numerica, 14, 1-137.
     */
    scalar_field_type estimateOmega_() const
    {
        using namespace Dune::Indices;
        auto& A = matrix_[_0][_0];
        auto& B = matrix_[_0][_1];
        auto& C = matrix_[_1][_0];

        U x(A.M());
        x = 1.0;

        scalar_field_type omega = 0.0;
        scalar_field_type lambdaMax = 0.0;

        static const auto iterations = Dumux::getParamFromGroup<std::size_t>(paramGroup_, "LinearSolver.Preconditioner.PowerLawIterations", 5);

        // apply power iteration x_k+1 = M*x_k/|M*x_k| for the matrix M = -C*Ainv*B
        for (std::size_t i = 0; i < iterations; ++i)
        {
            // bx = B*x
            U bx(x.size());
            B.mv(x, bx);

            // ainvbx = Ainv*(B*x)
            auto ainvbx = x;
            applySolverForA_(ainvbx, bx);

            // v = M*x = -C*(Ainv*B*x)
            U v(x.size());
            C.mv(ainvbx, v);
            v *= -1.0;

            // eigenvalue lambdaMax = xt*M*x/(xt*x) = xt*v/(xt*x);
            lambdaMax = x.dot(v)/(x.dot(x));

            // relaxation factor omega = 1/lambda;
            omega = 1.0/lambdaMax;

            // new iterate x = M*x/|M*x| = v/|v|
            x = v;
            x /= v.two_norm();
        }

        if (verbosity_ > 0)
        {
            std::cout << "\n*** Uzawa Preconditioner ***" << std::endl;
            std::cout << "Estimating relaxation factor based on Schur complement" << std::endl;
            std::cout << "Largest estimated eigenvalue lambdaMax = " << lambdaMax << std::endl;
            std::cout << "Relaxation factor omega = 1/lambdaMax = " << omega << std::endl;
        }

        return omega;
    }

    template<class Sol, class Rhs>
    void applySolverForA_(Sol& sol, Rhs& rhs) const
    {
        if (useDirectVelocitySolverForA_)
        {
#if HAVE_UMFPACK
            Dune::InverseOperatorResult res;
            umfPackSolverForA_->apply(sol, rhs, res);
#endif
        }
        else
        {
            amgSolverForA_->pre(sol, rhs);
            amgSolverForA_->apply(sol, rhs);
            amgSolverForA_->post(sol);
        }
    }

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The number of steps to do in apply
    const std::size_t numIterations_;
    //! \brief The relaxation factor to use
    scalar_field_type relaxationFactor_;
    //! \brief The verbosity level
    const int verbosity_;

    std::unique_ptr<AMGSolverForA> amgSolverForA_;
#if HAVE_UMFPACK
    std::unique_ptr<Dune::UMFPack<A>> umfPackSolverForA_;
#endif
    const std::string paramGroup_;
    const bool useDirectVelocitySolverForA_;
};

DUMUX_REGISTER_PRECONDITIONER("uzawa", Dumux::MultiTypeBlockMatrixPreconditionerTag, Dune::defaultPreconditionerBlockLevelCreator<Dumux::SeqUzawa, 1>());


/*! \brief Multihreaded Jacobi preconditioner.
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 *  \tparam l The block level to invert. Default is 1
 */
template<class M, class X, class Y, int l=1>
class ParMTJac : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;
    //! \brief real scalar type underlying the field_type
    typedef typename Dune::FieldTraits<scalar_field_type>::real_type real_field_type;

    //! \brief Constructor.
    ParMTJac (const M& A, int n, real_field_type w)
    : _A_(A), numSteps_(n), relaxationFactor_(w)
    { Dune::CheckIfDiagonalPresent<M,l>::check(_A_); }

    ParMTJac (const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& A, const Dune::ParameterTree& configuration)
    : ParMTJac(A->getmat(), configuration)
    {}

    ParMTJac (const M& A, const Dune::ParameterTree& configuration)
    : ParMTJac(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    {}

    void pre (X&, Y&) override {}

    void apply (X& update, const Y& defect) override
    {
        X xOld(update);
        const auto& A = _A_;

        for (int k=0; k<numSteps_; k++)
        {
            if (k > 0)
                xOld = update;

            Dumux::parallelFor(A.N(), [&](const std::size_t i)
            {
                auto& row = A[i];
                auto v = update[i];
                auto rhs = defect[i];
                const auto endColIt = row.end();
                auto colIt = row.begin();

                for (; colIt.index()<i; ++colIt)
                    colIt->mmv(xOld[colIt.index()], rhs);
                const auto diag = colIt;
                for (; colIt != endColIt; ++colIt)
                    colIt->mmv(xOld[colIt.index()], rhs);

                Dune::algmeta_itsteps<l-1, typename M::block_type>::dbjac(
                    *diag, v, rhs, relaxationFactor_
                );

                update[i].axpy(relaxationFactor_, v);
            });
        }
    }

    void post (X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    const M& _A_;
    int numSteps_;
    real_field_type relaxationFactor_;
};

namespace Detail {

// compute coloring for parallel sweep
template<class M>
void computeColorsForMatrixSweep_(const M& matrix, std::deque<std::vector<std::size_t>>& coloredIndices)
{
    // allocate some temporary memory
    std::vector<int> colors(matrix.N(), -1);
    std::vector<int> neighborColors; neighborColors.reserve(30);
    std::vector<bool> colorUsed; colorUsed.reserve(30);

    // a row that has my row index in their column set cannot have the same color
    // TODO this assumes a symmetric matrix pattern
    for (std::size_t i = 0; i < colors.size(); ++i)
    {
        neighborColors.clear();
        auto& row = matrix[i];
        const auto endColIt = row.end();
        for (auto colIt = row.begin(); colIt != endColIt; ++colIt)
            neighborColors.push_back(colors[colIt.index()]);

        const auto c = Detail::smallestAvailableColor(neighborColors, colorUsed);
        colors[i] = c;

        // add element to the set of elements with the same color
        if (c < coloredIndices.size())
            coloredIndices[c].push_back(i);
        else
            coloredIndices.emplace_back(1, i);
    }
}

// parallel SOR kernel (relaxed Gauss-Seidel)
template<bool forward, int l, class M, class X, class Y, class K>
void parallelBlockSOR_(const M& A, X& update, const Y& b, const K& relaxationFactor,
                       const std::deque<std::vector<std::size_t>>& coloredIndices)
{
    for (int color = 0; color < coloredIndices.size(); ++color)
    {
        const auto c = forward ? color : coloredIndices.size()-1-color;
        const auto& rowIndices = coloredIndices[c];
        Dumux::parallelFor(rowIndices.size(), [&](const std::size_t k)
        {
            const auto i = rowIndices[k];
            auto& row = A[i];
            auto v = update[i];
            auto rhs = b[i];
            const auto endColIt = row.end();
            auto colIt = row.begin();

            for (; colIt.index()<i; ++colIt)
                colIt->mmv(update[colIt.index()], rhs);
            const auto diag = colIt;
            for (; colIt != endColIt; ++colIt)
                colIt->mmv(update[colIt.index()], rhs);

            if constexpr (forward)
                Dune::algmeta_itsteps<l-1,typename M::block_type>::bsorf(
                    *diag, v, rhs, relaxationFactor
                );
            else
                Dune::algmeta_itsteps<l-1,typename M::block_type>::bsorb(
                    *diag, v, rhs, relaxationFactor
                );

            update[i].axpy(relaxationFactor, v);
        });
    }
}

} // end namespace Detail

/*! \brief Multihreaded SOR preconditioner
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 *  \tparam l The block level to invert. Default is 1
 */
template<class M, class X, class Y, int l=1>
class ParMTSOR : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;
    //! \brief real scalar type underlying the field_type
    typedef typename Dune::FieldTraits<scalar_field_type>::real_type real_field_type;

    //! \brief Constructor.
    ParMTSOR (const M& A, int n, real_field_type w)
    : _A_(A), numSteps_(n), relaxationFactor_(w)
    {
        Dune::CheckIfDiagonalPresent<M,l>::check(_A_);
        Detail::computeColorsForMatrixSweep_(_A_, colors_);
    }

    ParMTSOR (const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& A, const Dune::ParameterTree& configuration)
    : ParMTSOR(A->getmat(), configuration)
    {}

    ParMTSOR (const M& A, const Dune::ParameterTree& configuration)
    : ParMTSOR(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    {}

    void pre (X&, Y&) override {}

    void apply (X& v, const Y& d) override
    {
        this->template apply<true>(v,d);
    }

    template<bool forward>
    void apply(X& v, const Y& d)
    {
        for (int i=0; i<numSteps_; i++)
            Detail::parallelBlockSOR_<forward, l>(_A_, v, d, relaxationFactor_, colors_);
    }

    void post (X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    const M& _A_;
    int numSteps_;
    real_field_type relaxationFactor_;

    //! for each color a vector of row indices that can be dealt with in parallel
    std::deque<std::vector<std::size_t>> colors_;
};

/*! \brief Multithreaded SSOR preconditioner
 *
 *  \tparam M The matrix type to operate on
 *  \tparam X Type of the update
 *  \tparam Y Type of the defect
 *  \tparam l The block level to invert. Default is 1
 */
template<class M, class X, class Y, int l=1>
class ParMTSSOR : public Dune::Preconditioner<X,Y>
{
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;
    //! \brief real scalar type underlying the field_type
    typedef typename Dune::FieldTraits<scalar_field_type>::real_type real_field_type;

    //! \brief Constructor.
    ParMTSSOR (const M& A, int n, real_field_type w)
    : _A_(A), numSteps_(n), relaxationFactor_(w)
    {
        Dune::CheckIfDiagonalPresent<M,l>::check(_A_);
        Detail::computeColorsForMatrixSweep_(_A_, colors_);
    }

    ParMTSSOR (const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& A, const Dune::ParameterTree& configuration)
    : ParMTSSOR(A->getmat(), configuration)
    {}

    ParMTSSOR (const M& A, const Dune::ParameterTree& configuration)
    : ParMTSSOR(A, configuration.get<int>("iterations",1), configuration.get<real_field_type>("relaxation",1.0))
    {}

    void pre (X&, Y&) override {}

    void apply (X& v, const Y& d) override
    {
        for (int i=0; i<numSteps_; i++)
        {
            Detail::parallelBlockSOR_<true, l>(_A_, v, d, relaxationFactor_, colors_);
            Detail::parallelBlockSOR_<false, l>(_A_, v, d, relaxationFactor_, colors_);
        }
    }

    void post (X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    const M& _A_;
    int numSteps_;
    real_field_type relaxationFactor_;

    //! for each color a vector of row indices that can be dealt with in parallel
    std::deque<std::vector<std::size_t>> colors_;
};

} // end namespace Dumux

namespace Dune::Amg {

// make it possible to use Dumux::ParMTJac as AMG smoother
template<class M, class X, class Y, int l>
struct ConstructionTraits<Dumux::ParMTJac<M,X,Y,l> >
{
    using Arguments = DefaultConstructionArgs<SeqJac<M,X,Y,l> >;
    static inline std::shared_ptr<Dumux::ParMTJac<M,X,Y,l>> construct(Arguments& args)
    {
        return std::make_shared<Dumux::ParMTJac<M,X,Y,l>>(
            args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor
        );
    }
};

// make it possible to use Dumux::ParMTSOR as AMG smoother
template<class M, class X, class Y, int l>
struct ConstructionTraits<Dumux::ParMTSOR<M,X,Y,l> >
{
    using Arguments = DefaultConstructionArgs<SeqSOR<M,X,Y,l> >;
    static inline std::shared_ptr<Dumux::ParMTSOR<M,X,Y,l>> construct(Arguments& args)
    {
        return std::make_shared<Dumux::ParMTSOR<M,X,Y,l>>(
            args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor
        );
    }
};

template<class M, class X, class Y, int l>
struct SmootherApplier<Dumux::ParMTSOR<M,X,Y,l> >
{
    typedef Dumux::ParMTSOR<M,X,Y,l> Smoother;
    typedef typename Smoother::range_type Range;
    typedef typename Smoother::domain_type Domain;

    static void preSmooth(Smoother& smoother, Domain& v, Range& d)
    { smoother.template apply<true>(v,d); }

    static void postSmooth(Smoother& smoother, Domain& v, Range& d)
    { smoother.template apply<false>(v,d); }
};

// make it possible to use Dumux::ParMTSSOR as AMG smoother
template<class M, class X, class Y, int l>
struct ConstructionTraits<Dumux::ParMTSSOR<M,X,Y,l> >
{
    using Arguments = DefaultConstructionArgs<SeqSSOR<M,X,Y,l> >;
    static inline std::shared_ptr<Dumux::ParMTSSOR<M,X,Y,l>> construct(Arguments& args)
    {
        return std::make_shared<Dumux::ParMTSSOR<M,X,Y,l>>(
            args.getMatrix(), args.getArgs().iterations, args.getArgs().relaxationFactor
        );
    }
};

} // end namespace Dune::Amg

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A preconditioner based on the preconditioned Uzawa algorithm for saddle-point problems of the form
 * \f$
 \begin{pmatrix} A & B \\ C & 0 \end{pmatrix}
 \begin{pmatrix} u \\ p \end{pmatrix} =
 \begin{pmatrix} f \\ g \end{pmatrix}
 * \f$
 *
 * This preconditioner is especially suited for solving the incompressible (Navier-)Stokes equations.
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level (for compatibility reasons, unused).
 */
template<class M, class X, class Y, class Comm, int l = 1>
class IncompressibleStokesPreconditioner : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "SeqUzawa expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l== 1, "IncompressibleStokesPreconditioner expects a block level of 1.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using P = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_1][Dune::Indices::_1])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;
    //! \brief the type of the pressure operator
    using PressureLinearOperator = Dune::OverlappingSchwarzOperator<P,V,V,Comm>;

    /*!
     * \brief Constructor
     *
     * \param mat The matrix to operate on.
     * \param params Collection of paramters.
     */
    IncompressibleStokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op,
        const std::shared_ptr<const Dune::AssembledLinearOperator<P,V,V>>& pop,
        const std::array<std::shared_ptr<const Comm>, 2>& comms,
        const Dune::ParameterTree& params
    )
    : IncompressibleStokesPreconditioner(op, pop, pop, comms, params)
    {
        useMultiOperator_ = false;
    }

    IncompressibleStokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op,
        const std::shared_ptr<const Dune::AssembledLinearOperator<P,V,V>>& pop,
        const std::shared_ptr<const Dune::AssembledLinearOperator<P,V,V>>& pop2,
        const std::array<std::shared_ptr<const Comm>, 2>& comms,
        const Dune::ParameterTree& params
    )
    : matrix_(op->getmat())
    , pmatrix_(pop->getmat())
    , p2matrix_(pop2->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , comms_(comms)
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    {
        initPreconditioner_(params);
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    virtual void pre(X& x, Y& b) {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     *
     * The currentDefect has be be in a consistent representation,
     * Definition 2.3 Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112
     * The update is initially zero. At exit the update has to be
     * in a consistent representation. This usually requires communication.
     */
    virtual void apply(X& update, const Y& currentDefect)
    {
        using namespace Dune::Indices;

        const auto& B = matrix_[_0][_1];
        const auto& C = matrix_[_1][_0];
        const auto& f = currentDefect[_0];
        const auto& g = currentDefect[_1];
        auto& u = update[_0];
        auto& p = update[_1];

        // The preconditioned Uzawa iteration
        //
        // in comparison to the classical SeqUzawa, we precondition the
        // pressure update. Also, here we assume D = 0 and that the initial
        // update is zero.
        // Depending on the problem, the pressure preconditioner might be an
        // inverse weighted mass matrix or a (compatible) Laplacian, or
        // a weighted combination, see https://doi.org/10.1002/fld.1650080802

        // u_k+1 = u_k + Q_A^−1*f,
        auto uRhs = f;
        applyPreconditionerForA_(u, uRhs);

        // p_k+1 = p_k + Q_P^-1*(g - C*u_k+1)
        auto pRhs = g;
        C.mmv(u, pRhs);
        applyPreconditionerForP_(p, pRhs);

        // full Schur complement style update (symmetric)
        // would require another application of Q_A^−1
        // // u_k+2 = u_k+1 - Q_A^−1*(Bp_k+1)
        // uRhs = 0;
        // B.mv(p, uRhs);
        // uIncrement = 0;
        // applyPreconditionerForA_(uIncrement, uRhs);
        // u -= uIncrement; // uIncrement = -(u_k+2 - u_k+1)

        // instead we only do a projection
        // project velocity into divergence-free space
        // u_k+2 = u_k+1 + Bp_k+1
        auto uHalf = u;
        // when using this in parallel:
        // uHalf is already consistent
        // we have to make u consistent after applying B
        // i.e. sum up over the border entities
        u = 0.0;
        B.umv(p, u);
        comms_[0]->addOwnerCopyToOwnerCopy(u, u);
        u += uHalf;

        // // update pressure
        // // p_k+2 = p_k+1 + Cu_k+1 - g
        // C.umv(uHalf, p);
        // p -= g;
    }

    /*!
     * \brief Clean up.
     */
    virtual void post(X& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        // return Dune::SolverCategory::sequential;
        return Dune::SolverCategory::overlapping;
    }

private:
    void initPreconditioner_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.DirectSolverForVelocity", false))
        {
            // directSolver_ = std::make_shared<Dune::UMFPack<A>>(matrix_[_0][_0], verbosity_);
            // using Wrap = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<U, U>>;
            // preconditionerForA_ = std::make_shared<Wrap>(*directSolver_);
        }
        else
        {
            // using VelLinearOperator = Dune::MatrixAdapter<A, U, U>;
            // auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0]);
            // using AMGForA = Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSOR<A, U, U>, Comm>;
            using VelLinearOperator = Dune::NonoverlappingSchwarzOperator<A, U, U, Comm>;
            auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0], *comms_[0]);
            auto creator = Dune::AMGCreator();
            preconditionerForA_ = creator.makeAMG(lopV, "sor", params);
        }

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.AmgForPressure", false))
        {
            using PressLinearOperator = Dune::OverlappingSchwarzOperator<P, V, V, Comm>;
            auto lopP = std::make_shared<PressLinearOperator>(pmatrix_, *comms_[1]);
            auto creator = Dune::AMGCreator();
            preconditionerForP_ = creator.makeAMG(lopP, "ssor", params);

            if (useMultiOperator_)
            {
                auto lopP2 = std::make_shared<PressLinearOperator>(p2matrix_, *comms_[1]);
                using PressJacobi = Dune::SeqJac<P, V, V>;
                using ParPressJacobi = Dune::BlockPreconditioner<V, V, Comm, PressJacobi>;
                auto seqPre = std::make_shared<PressJacobi>(lopP2, params);
                preconditionerForP2_ = std::make_shared<ParPressJacobi>(seqPre, *comms_[1]);
            }
        }
        else
        {
            using PressLinearOperator = Dune::OverlappingSchwarzOperator<P, V, V, Comm>;
            auto lopP = std::make_shared<PressLinearOperator>(pmatrix_, *comms_[1]);
            using PressJacobi = Dune::SeqJac<P, V, V>;
            using ParPressJacobi = Dune::BlockPreconditioner<V, V, Comm, PressJacobi>;
            auto seqPre = std::make_shared<PressJacobi>(lopP, params);
            preconditionerForP_ = std::make_shared<ParPressJacobi>(seqPre, *comms_[1]);

            if (useMultiOperator_)
            {
                auto lopP2 = std::make_shared<PressLinearOperator>(p2matrix_, *comms_[1]);
                using PressJacobi = Dune::SeqJac<P, V, V>;
                using ParPressJacobi = Dune::BlockPreconditioner<V, V, Comm, PressJacobi>;
                auto seqPre2 = std::make_shared<PressJacobi>(lopP2, params);
                preconditionerForP2_ = std::make_shared<ParPressJacobi>(seqPre2, *comms_[1]);
            }
        }
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForA_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForA_->pre(sol, rhs);
        preconditionerForA_->apply(sol, rhs);
        preconditionerForA_->post(sol);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForP_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForP_->pre(sol, rhs);
        preconditionerForP_->apply(sol, rhs);
        preconditionerForP_->post(sol);

        if (useMultiOperator_)
        {
            auto sol2 = sol; sol2 = 0.0;
            preconditionerForP2_->pre(sol2, rhs);
            preconditionerForP2_->apply(sol2, rhs);
            preconditionerForP2_->post(sol2);
            sol += sol2;
        }
    }

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The matrix we operate on.
    const P& pmatrix_;
    //! \brief The matrix we operate on.
    const P& p2matrix_;
    //! \brief The verbosity level
    const int verbosity_;

    std::array<std::shared_ptr<const Comm>, 2> comms_;

    std::shared_ptr<Dune::Preconditioner<U, U>> preconditionerForA_;
    std::shared_ptr<Dune::Preconditioner<V, V>> preconditionerForP_;
    std::shared_ptr<Dune::Preconditioner<V, V>> preconditionerForP2_;
    std::shared_ptr<Dune::InverseOperator<U, U>> directSolver_;
    const std::string paramGroup_;
    bool useMultiOperator_ = true;
};

} // end namespace Dumux

#endif
