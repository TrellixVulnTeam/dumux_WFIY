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
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes
 *        with multiple domains
 */
#ifndef DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FV_ASSEMBLER_HH

#if HAVE_TBB
#include "tbb/parallel_for.h"
// Maybe change this to true at some point,
// together with some heuristic based on number of elements, scheme,...
// #define DUMUX_MULTITHREADING_ASSEMBLY_DEFAULT true
#endif

#ifndef DUMUX_MULTITHREADING_ASSEMBLY_DEFAULT
#define DUMUX_MULTITHREADING_ASSEMBLY_DEFAULT false
#endif

#include <type_traits>
#include <tuple>

#include <dune/common/hybridutilities.hh>
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/discretization/method.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/linear/parallelhelpers.hh>

#include "couplingjacobianpattern.hh"
#include "subdomaincclocalassembler.hh"
#include "subdomainboxlocalassembler.hh"
#include "subdomainstaggeredlocalassembler.hh"
#include "subdomainfclocalassembler.hh"

#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup Assembly
 * \brief A linear system assembler (residual and Jacobian) for finite volume schemes (box, tpfa, mpfa, ...)
 *        with multiple domains
 * \tparam MDTraits the multidimension traits
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam useImplicitAssembly if to use an implicit or explicit time discretization
 */
template<class MDTraits, class CMType, DiffMethod diffMethod, bool useImplicitAssembly = true>
class MultiDomainFVAssembler
{
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomain<id>::TypeTag;

public:
    using Traits = MDTraits;

    using Scalar = typename MDTraits::Scalar;

    //! TODO get rid of this GetPropType
    template<std::size_t id>
    using LocalResidual = GetPropType<SubDomainTypeTag<id>, Properties::LocalResidual>;

    template<std::size_t id>
    using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;

    template<std::size_t id>
    using GridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

    template<std::size_t id>
    using ElementSeed = typename GridGeometry<id>::GridView::Grid::template Codim<0>::EntitySeed;

    template<std::size_t id>
    using ElementSets = std::deque<std::vector<ElementSeed<id>>>;

    template<std::size_t id>
    using Problem = typename MDTraits::template SubDomain<id>::Problem;

    using JacobianMatrix = typename MDTraits::JacobianMatrix;
    using SolutionVector = typename MDTraits::SolutionVector;
    using ResidualType = SolutionVector;

    using CouplingManager = CMType;

    /*!
     * \brief Returns true if the assembler considers implicit assembly.
     */
    static constexpr bool isImplicit()
    { return useImplicitAssembly; }

private:

    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using GridGeometryTuple = typename MDTraits::template TupleOfSharedPtrConst<GridGeometry>;
    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;
    using ElementSetsTuple = typename MDTraits::template Tuple<ElementSets>;

    using TimeLoop = TimeLoopBase<Scalar>;
    using ThisType = MultiDomainFVAssembler<MDTraits, CouplingManager, diffMethod, isImplicit()>;

    template<DiscretizationMethod discMethod, std::size_t id>
    struct SubDomainAssemblerType;

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::cctpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::ccmpfa, id>
    {
        using type = SubDomainCCLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::box, id>
    {
        using type = SubDomainBoxLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::staggered, id>
    {
        using type = SubDomainStaggeredLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    struct SubDomainAssemblerType<DiscretizationMethod::fcstaggered, id>
    {
        using type = SubDomainFaceCenteredLocalAssembler<id, SubDomainTypeTag<id>, ThisType, diffMethod, isImplicit()>;
    };

    template<std::size_t id>
    using SubDomainAssembler = typename SubDomainAssemblerType<GridGeometry<id>::discMethod, id>::type;

public:


    /*!
     * \brief The constructor for stationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple&& problem,
                           GridGeometryTuple&& gridGeometry,
                           GridVariablesTuple&& gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager)
    : couplingManager_(couplingManager)
    , problemTuple_(problem)
    , gridGeometryTuple_(gridGeometry)
    , gridVariablesTuple_(gridVariables)
    , timeLoop_()
    , isStationaryProblem_(true)
    , warningIssued_(false)
    {
        static_assert(isImplicit(), "Explicit assembler for stationary problem doesn't make sense!");
        std::cout << "Instantiated assembler for a stationary problem." << std::endl;
        enableMultithreading_ = getParam<bool>("Assembly.Multithreading", DUMUX_MULTITHREADING_ASSEMBLY_DEFAULT);
    }

    /*!
     * \brief The constructor for instationary problems
     * \note the grid variables might be temporarily changed during assembly (if caching is enabled)
     *       it is however guaranteed that the state after assembly will be the same as before
     */
    MultiDomainFVAssembler(ProblemTuple&& problem,
                           GridGeometryTuple&& gridGeometry,
                           GridVariablesTuple&& gridVariables,
                           std::shared_ptr<CouplingManager> couplingManager,
                           std::shared_ptr<const TimeLoop> timeLoop,
                           const SolutionVector& prevSol)
    : couplingManager_(couplingManager)
    , problemTuple_(problem)
    , gridGeometryTuple_(gridGeometry)
    , gridVariablesTuple_(gridVariables)
    , timeLoop_(timeLoop)
    , prevSol_(&prevSol)
    , isStationaryProblem_(false)
    , warningIssued_(false)
    {
        std::cout << "Instantiated assembler for an instationary problem." << std::endl;
        enableMultithreading_ = getParam<bool>("Assembly.Multithreading", DUMUX_MULTITHREADING_ASSEMBLY_DEFAULT);
    }

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        checkAssemblerState_();
        resetJacobian_();
        resetResidual_();

        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainId)
        {
            auto& jacRow = (*jacobian_)[domainId];
            auto& subRes = (*residual_)[domainId];
            this->assembleJacobianAndResidual_(domainId, jacRow, subRes, curSol);

            const auto gridGeometry = std::get<domainId>(gridGeometryTuple_);
            enforcePeriodicConstraints_(domainId, jacRow, subRes, *gridGeometry, curSol[domainId]);
        });
    }

    //! compute the residuals using the internal residual
    void assembleResidual(const SolutionVector& curSol)
    {
        resetResidual_();
        assembleResidual(*residual_, curSol);
    }

    //! assemble a residual r
    void assembleResidual(ResidualType& r, const SolutionVector& curSol)
    {
        r = 0.0;

        checkAssemblerState_();

        // update the grid variables for the case of active caching
        updateGridVariables(curSol);

        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(r)), [&](const auto domainId)
        {
            auto& subRes = r[domainId];
            this->assembleResidual_(domainId, subRes, curSol);
        });
    }

    //! compute the residual and return it's vector norm
    Scalar residualNorm(const SolutionVector& curSol)
    {
        ResidualType residual;
        setResidualSize(residual);
        assembleResidual(residual, curSol);

        // calculate the squared norm of the residual
        Scalar resultSquared = 0.0;

        // for box communicate the residual with the neighboring processes
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(residual)), [&](const auto domainId)
        {
            const auto gridGeometry = std::get<domainId>(gridGeometryTuple_);
            const auto& gridView = gridGeometry->gridView();

            if (gridView.comm().size() > 1 && gridView.overlapSize(0) == 0)
            {
                if constexpr (GridGeometry<domainId>::discMethod == DiscretizationMethod::box)
                {
                    using GV = typename GridGeometry<domainId>::GridView;
                    using DM = typename GridGeometry<domainId>::VertexMapper;
                    using PVHelper = ParallelVectorHelper<GV, DM, GV::dimension>;

                    PVHelper vectorHelper(gridView, gridGeometry->vertexMapper());

                    vectorHelper.makeNonOverlappingConsistent(residual[domainId]);
                }
            }
            else if (!warningIssued_)
            {
                if (gridView.comm().size() > 1 && gridView.comm().rank() == 0)
                    std::cout << "\nWarning: norm calculation adds entries corresponding to\n"
                              << "overlapping entities multiple times. Please use the norm\n"
                              << "function provided by a linear solver instead." << std::endl;

                warningIssued_ = true;
            }

            Scalar localNormSquared = residual[domainId].two_norm2();

            if (gridView.comm().size() > 1)
            {
                localNormSquared = gridView.comm().sum(localNormSquared);
            }

            resultSquared += localNormSquared;
        });

        using std::sqrt;
        return sqrt(resultSquared);
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern(*jacobian_);
        setResidualSize(*residual_);

        if (enableMultithreading_)
            buildElementSets_();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler if you don't need them outside this class
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();
        residual_ = std::make_shared<SolutionVector>();

        setJacobianBuildMode(*jacobian_);
        setJacobianPattern(*jacobian_);
        setResidualSize(*residual_);

        if (enableMultithreading_)
            buildElementSets_();
    }

    /*!
     * \brief Sets the jacobian build mode
     */
    void setJacobianBuildMode(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto i)
        {
            forEach(jac[i], [&](auto& jacBlock)
            {
                using BlockType = std::decay_t<decltype(jacBlock)>;
                if (jacBlock.buildMode() == BlockType::BuildMode::unknown)
                    jacBlock.setBuildMode(BlockType::BuildMode::random);
                else if (jacBlock.buildMode() != BlockType::BuildMode::random)
                    DUNE_THROW(Dune::NotImplemented, "Only BCRS matrices with random build mode are supported at the moment");
            });
        });
    }

    /*!
     * \brief Sets the jacobian sparsity pattern.
     */
    void setJacobianPattern(JacobianMatrix& jac) const
    {
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto domainI)
        {
            forEach(integralRange(Dune::Hybrid::size(jac[domainI])), [&](const auto domainJ)
            {
                const auto pattern = this->getJacobianPattern_(domainI, domainJ);
                pattern.exportIdx(jac[domainI][domainJ]);
            });
        });
    }

    /*!
     * \brief Resizes the residual
     */
    void setResidualSize(SolutionVector& res) const
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(res)), [&](const auto domainId)
        { res[domainId].resize(this->numDofs(domainId)); });
    }

    /*!
     * \brief Updates the grid variables with the given solution
     */
    void updateGridVariables(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).update(curSol[domainId]); });
    }

    /*!
     * \brief Resets the grid variables to the last time step
     */
    void resetTimeStep(const SolutionVector& curSol)
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(gridVariablesTuple_)), [&](const auto domainId)
        { this->gridVariables(domainId).resetTimeStep(curSol[domainId]); });
    }

    //! the number of dof locations of domain i
    template<std::size_t i>
    std::size_t numDofs(Dune::index_constant<i> domainId) const
    { return std::get<domainId>(gridGeometryTuple_)->numDofs(); }

    //! the problem of domain i
    template<std::size_t i>
    const auto& problem(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(problemTuple_); }

    //! the finite volume grid geometry of domain i
    template<std::size_t i>
    const auto& gridGeometry(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridGeometryTuple_); }

    //! the grid view of domain i
    template<std::size_t i>
    const auto& gridView(Dune::index_constant<i> domainId) const
    { return gridGeometry(domainId).gridView(); }

    //! the grid variables of domain i
    template<std::size_t i>
    GridVariables<i>& gridVariables(Dune::index_constant<i> domainId)
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the grid variables of domain i
    template<std::size_t i>
    const GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariablesTuple_); }

    //! the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! the full Jacobian matrix
    JacobianMatrix& jacobian()
    { return *jacobian_; }

    //! the full residual vector
    SolutionVector& residual()
    { return *residual_; }

    //! the solution of the previous time step
    const SolutionVector& prevSol() const
    { return *prevSol_; }

    /*!
     * \brief Set time loop for instationary problems
     * \note calling this turns this into a stationary assembler
     */
    void setTimeManager(std::shared_ptr<const TimeLoop> timeLoop)
    { timeLoop_ = timeLoop; isStationaryProblem_ = !(static_cast<bool>(timeLoop)); }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { prevSol_ = &u; }

    /*!
     * \brief Whether we are assembling a stationary or instationary problem
     */
    bool isStationaryProblem() const
    { return isStationaryProblem_; }

    /*!
     * \brief Create a local residual object (used by the local assembler)
     */
    template<std::size_t i>
    LocalResidual<i> localResidual(Dune::index_constant<i> domainId) const
    { return LocalResidual<i>(std::get<domainId>(problemTuple_).get(), timeLoop_.get()); }

protected:
    //! the coupling manager coupling the sub domains
    std::shared_ptr<CouplingManager> couplingManager_;

private:
    // reset the residual vector to 0.0
    void resetResidual_()
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize(*residual_);
        }

        (*residual_) = 0.0;
    }

    // reset the jacobian vector to 0.0
    void resetJacobian_()
    {
        if(!jacobian_)
        {
            jacobian_ = std::make_shared<JacobianMatrix>();
            setJacobianBuildMode(*jacobian_);
            setJacobianPattern(*jacobian_);
        }

       (*jacobian_)  = 0.0;
    }

    // check if the assembler is in a correct state for assembly
    void checkAssemblerState_() const
    {
        if (!isStationaryProblem_ && !prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        if (isStationaryProblem_ && prevSol_)
            DUNE_THROW(Dune::InvalidStateException, "Assembling stationary problem but a previous solution was set."
                                                    << " Did you forget to set the timeLoop to make this problem instationary?");
    }

    template<std::size_t i, class JacRow, class SubRes>
    void assembleJacobianAndResidual_(Dune::index_constant<i> domainId, JacRow& jacRow, SubRes& subRes,
                                      const SolutionVector& curSol)
    {
        assemble_(domainId, [&](const auto& element)
        {
            SubDomainAssembler<i> subDomainAssembler(*this, element, curSol, *couplingManager_);
            subDomainAssembler.assembleJacobianAndResidual(jacRow, subRes, gridVariablesTuple_);
        });
    }

    template<std::size_t i, class SubRes>
    void assembleResidual_(Dune::index_constant<i> domainId, SubRes& subRes,
                           const SolutionVector& curSol)
    {
        assemble_(domainId, [&](const auto& element)
        {
            SubDomainAssembler<i> subDomainAssembler(*this, element, curSol, *couplingManager_);
            subDomainAssembler.assembleResidual(subRes);
        });
    }

    /*!
     * \brief A method assembling something per element
     * \note Handles exceptions for parallel runs
     * \throws NumericalProblem on all processes if something throwed during assembly
     * TODO: assemble in parallel
     */
    template<std::size_t i, class AssembleElementFunc>
    void assemble_(Dune::index_constant<i> domainId, AssembleElementFunc&& assembleElement) const
    {
        // let the local assembler add the element contributions
        if (enableMultithreading_)
        {
#if HAVE_TBB
            for (const auto& elements : std::get<i>(elementSets_))
            {
                tbb::parallel_for(tbb::blocked_range<std::size_t>(0, elements.size()),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (std::size_t i=r.begin(); i<r.end(); ++i)
                    {
                        const auto element = gridView(domainId).grid().entity(elements[i]);
                        assembleElement(element);
                    }
                });
            }
#else
            DUNE_THROW(Dune::Exception,
                "Multithread assembly has been explicitly requested but TBB has not been found!");
#endif // HAVE_TBB
        }
        else
        {
            for (const auto& element : elements(gridView(domainId)))
                assembleElement(element);
        }
    }

    // get diagonal block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i==j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        const auto& gg = gridGeometry(domainI);
        auto pattern = getJacobianPattern<isImplicit()>(gg);
        couplingManager_->extendJacobianPattern(domainI, pattern);
        return pattern;
    }

    // get coupling block pattern
    template<std::size_t i, std::size_t j, typename std::enable_if_t<(i!=j), int> = 0>
    Dune::MatrixIndexSet getJacobianPattern_(Dune::index_constant<i> domainI,
                                             Dune::index_constant<j> domainJ) const
    {
        return getCouplingJacobianPattern<isImplicit()>(*couplingManager_,
                                                        domainI, gridGeometry(domainI),
                                                        domainJ, gridGeometry(domainJ));
    }

    // build periodic contraints into the system matrix
    template<std::size_t i, class JacRow, class Sol, class GG>
    void enforcePeriodicConstraints_(Dune::index_constant<i> domainI, JacRow& jacRow, Sol& res, const GG& gridGeometry, const Sol& curSol)
    {
        if constexpr (GG::discMethod == DiscretizationMethod::box || GG::discMethod == DiscretizationMethod::fcstaggered)
        {
            for (const auto& m : gridGeometry.periodicVertexMap())
            {
                if (m.first < m.second)
                {
                    auto& jac = jacRow[domainI];

                    // add the second row to the first
                    res[m.first] += res[m.second];

                    // enforce the solution of the first periodic DOF to the second one
                    res[m.second] = curSol[m.second] - curSol[m.first];

                    const auto end = jac[m.second].end();
                    for (auto it = jac[m.second].begin(); it != end; ++it)
                        jac[m.first][it.index()] += (*it);

                    // enforce constraint in second row
                    for (auto it = jac[m.second].begin(); it != end; ++it)
                        (*it) = it.index() == m.second ? 1.0 : it.index() == m.first ? -1.0 : 0.0;

                    using namespace Dune::Hybrid;
                    forEach(makeIncompleteIntegerSequence<JacRow::size(), domainI>(), [&](const auto couplingDomainId)
                    {
                        auto& jacCoupling = jacRow[couplingDomainId];

                        for (auto it = jacCoupling[m.second].begin(); it != jacCoupling[m.second].end(); ++it)
                            jacCoupling[m.first][it.index()] += (*it);

                        for (auto it = jacCoupling[m.second].begin(); it != jacCoupling[m.second].end(); ++it)
                            (*it) = 0.0;
                    });
                }
            }
        }
    }

    //! build element sets for lock-free multithreaded assembly
    void buildElementSets_(int verbosity = 1)
    {
        // for multidomain, we need to make sure to consider the entire row
        // hence all conflicts arising from writing into the coupling blocks too
        // - we assume for now that the coupling manager creates local contexts per color so there is no conflict there
        // - when writing into the coupling block there is no conflict as we only touch rows related to our own element
        // - however, when deflecting the solution vector, there is a conflict with all other elements that couple
        //   with the same dof, so we need to consider the coupling stencil in the set
        // - in the case of caching deflecting the solution may require updating caches in the local stencil of the other domain
        // - if there is an extended own stencil we have added restrictions for the own domain too
        // - all coupling is binary, because we still assemble domain after domain sequentially
        //
        // working theory: if we construct a dofToElement map (dof J to element I) for every subdomain J, we can iterate over
        // all domains and then find conflicting elements in I by checking if they are in the dofToElement map of
        // any of the dofs that we interact with (own stencil + extended stencil + coupling stencil)
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<JacobianMatrix::N()>(), [&](const auto i)
        {
            Dune::Timer timer;

            constexpr auto domainI = Dune::index_constant<i>;
            const auto& gg = gridGeometry(domainI);
            const auto dofToElement = [&]
            {
                std::array<std::vector<std::vector<std::size_t>>, Traits::numSubDomains> dofToElement;
                forEach(std::make_index_sequence<JacobianMatrix::M()>(), [&](const auto j)
                {
                    dofToElement[j] = Detail::computeDofToElementMap(
                        domainI, domainJ, couplingManager(), gg
                    );
                });
            }();

            std::vector<int> colors(gg.gridView().size(0), -1);

            // reserve some memory for helper arrays
            std::vector<int> neighborColors; neighborColors.reserve(50);
            std::vector<bool> notAssigned; notAssigned.reserve(50);

            for (const auto& element : elements(gg.gridView()))
            {
                // compute neighbor colors based on discretization-dependent stencil
                neighborColors.clear();

                // neighbors due to conflicts in own domain (regular)
                Detail::addNeighborColors(gg, element, colors, dofToElement[i], neighborColors);
                // neighbors due to conflicts in own domain (extended stencil)
                Detail::addNeighborColors(gg, element, colors, dofToElement[i], neighborColors);
                // neighbors due to conflicts in other domains (coupling stencils)
                forEach(std::make_index_sequence<JacobianMatrix::M()>(), [&](const auto j){
                    Detail::addNeighborColors(gg, element, colors, dofToElement[i], neighborColors);
                });

                // find smallest color (positive integer) not in neighborColors
                const auto color = Detail::smallestAvailableColor(neighborColors, notAssigned);

                // assign color to element
                colors[gg.elementMapper().index(element)] = color;

                // add element to the set of elements with the same color
                auto& elementSets = std::get<i>(elementSets_);
                if (color < elementSets.size())
                    elementSets[color].push_back(element.seed());
                else
                    elementSets.push_back(std::vector<ElementSeed>{ element.seed() });
            }

            if (verbosity > 0)
                std::cout << Fmt::format("Colored {} elements of domain {} with {} colors in {} seconds.\n",
                                        gridView(domainI).size(0), i, elementSets.size(), timer.elapsed());
        });
    }

    //! pointer to the problem to be solved
    ProblemTuple problemTuple_;

    //! the finite volume geometry of the grid
    GridGeometryTuple gridGeometryTuple_;

    //! the variables container for the grid
    GridVariablesTuple gridVariablesTuple_;

    //! the time loop for instationary problem assembly
    std::shared_ptr<const TimeLoop> timeLoop_;

    //! an observing pointer to the previous solution for instationary problems
    const SolutionVector* prevSol_ = nullptr;

    //! if this assembler is assembling an instationary problem
    bool isStationaryProblem_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    //! Issue a warning if the calculation is used in parallel with overlap. This could be a static local variable if it wasn't for g++7 yielding a linker error.
    bool warningIssued_;

    //! element sets for parallel assembly
    bool enableMultithreading_ = false;
    ElementSetsTuple elementSets_;
};

} // end namespace Dumux

#endif
