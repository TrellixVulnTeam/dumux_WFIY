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
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (pq1bubble method)
 */
#ifndef DUMUX_ASSEMBLY_PQ1BUBBLE_LOCAL_ASSEMBLER_HH
#define DUMUX_ASSEMBLY_PQ1BUBBLE_LOCAL_ASSEMBLER_HH

#include <dune/grid/common/gridenums.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>
#include <dumux/assembly/entitycolor.hh>
#include <dumux/assembly/boxlocalassembler.hh>
#include <dumux/discretization/pq1bubble/elementsolution.hh>
#include <dumux/assembly/partialreassembler.hh>

#include "volvardeflectionhelper_.hh"

namespace Dumux {

namespace Detail::PQ1Bubble {

struct NoOperator
{
    template<class... Args>
    constexpr void operator()(Args&&...) const {}
};

template<class X, class Y>
using Impl = std::conditional_t<!std::is_same_v<X, void>, X, Y>;

} // end namespace Detail

/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief A base class for all local pq1bubble assemblers
 * \tparam TypeTag The TypeTag
 * \tparam Assembler The assembler type
 * \tparam Implementation The actual implementation
 * \tparam implicit Specifies whether the time discretization is implicit or not (i.e. explicit)
 */
template<class TypeTag, class Assembler, class Implementation, bool implicit>
class PQ1BubbleLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler, Implementation, implicit>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr int dim = GridView::dimension;

public:

    using ParentType::ParentType;

    void bindLocalViews()
    {
        ParentType::bindLocalViews();
        this->elemBcTypes().update(this->asImp_().problem(), this->element(), this->fvGeometry());
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template <class PartialReassembler = DefaultPartialReassembler, class CouplingFunction = Detail::PQ1Bubble::NoOperator>
    void assembleJacobianAndResidual(JacobianMatrix& jac, SolutionVector& res, GridVariables& gridVariables,
                                     const PartialReassembler* partialReassembler,
                                     const CouplingFunction& maybeAssembleCouplingBlocks = CouplingFunction())
    {
        static_assert(!std::decay_t<decltype(this->asImp_().problem())>::enableInternalDirichletConstraints(), "Internal Dirichlet constraints are currently not implemented for face-centered staggered models!");

        this->asImp_().bindLocalViews();
        const auto& gridGeometry = this->asImp_().problem().gridGeometry();
        const auto eIdxGlobal = gridGeometry.elementMapper().index(this->element());
        if (partialReassembler
            && partialReassembler->elementColor(eIdxGlobal) == EntityColor::green)
        {
            const auto residual = this->asImp_().evalLocalResidual(); // forward to the internal implementation
            for (const auto& scv : scvs(this->fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else if (!this->elementIsGhost())
        {
            const auto residual = this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables, partialReassembler); // forward to the internal implementation
            for (const auto& scv : scvs(this->fvGeometry()))
                res[scv.dofIndex()] += residual[scv.localDofIndex()];

            // assemble the coupling blocks for coupled models (does nothing if not coupled)
            maybeAssembleCouplingBlocks(residual);
        }
        else
        {
            assert(this->elementIsGhost());

            // handle vertex dofs
            const auto numVerticesLocal = this->element().subEntities(dim);
            for (int i = 0; i < numVerticesLocal; ++i)
            {
                // do not change the non-ghost vertices
                auto vertex = this->element().template subEntity<dim>(i);
                if (vertex.partitionType() == Dune::InteriorEntity || vertex.partitionType() == Dune::BorderEntity)
                    continue;

                // set main diagonal entries for the vertex
                const auto dofIndex = gridGeometry.dofMapper().index(vertex);

                // this might be a vector-valued dof
                typedef typename JacobianMatrix::block_type BlockType;
                BlockType &J = jac[dofIndex][dofIndex];
                for (int j = 0; j < BlockType::rows; ++j)
                    J[j][j] = 1.0;

                // set residual for the ghost vertex dof
                res[dofIndex] = 0;
            }

            // handle element dof
            const auto elemDofIndex = gridGeometry.dofMapper().index(this->element());

            // this might be a vector-valued dof
            typedef typename JacobianMatrix::block_type BlockType;
            BlockType &J = jac[elemDofIndex][elemDofIndex];
            for (int j = 0; j < BlockType::rows; ++j)
                J[j][j] = 1.0;

            // set residual for the ghost element dof
            res[elemDofIndex] = 0;
        }

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];

            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;


            // if a periodic dof has Dirichlet values also apply the same Dirichlet values to the other dof TODO periodic
            // if (this->assembler().gridGeometry().dofOnPeriodicBoundary(scvI.dofIndex()))
            // {
            //     const auto periodicDof = this->assembler().gridGeometry().periodicallyMappedDof(scvI.dofIndex());
            //     res[periodicDof][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
            //     const auto end = jac[periodicDof].end();
            //     for (auto it = jac[periodicDof].begin(); it != end; ++it)
            //         (*it) = periodicDof != it.index() ? 0.0 : 1.0;
            // }
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     */
    void assembleJacobian(JacobianMatrix& jac, GridVariables& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->asImp_().assembleJacobianAndResidualImpl(jac, gridVariables); // forward to the internal implementation

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            auto& row = jac[scvI.dofIndex()];
            for (auto col = row.begin(); col != row.end(); ++col)
                row[col.index()][eqIdx] = 0.0;

            jac[scvI.dofIndex()][scvI.dofIndex()][eqIdx][pvIdx] = 1.0;
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Assemble the residual only
     */
    void assembleResidual(SolutionVector& res)
    {
        this->asImp_().bindLocalViews();
        const auto residual = this->evalLocalResidual();

        for (const auto& scv : scvs(this->fvGeometry()))
            res[scv.dofIndex()] += residual[scv.localDofIndex()];

        auto applyDirichlet = [&] (const auto& scvI,
                                   const auto& dirichletValues,
                                   const auto eqIdx,
                                   const auto pvIdx)
        {
            res[scvI.dofIndex()][eqIdx] = this->curElemVolVars()[scvI].priVars()[pvIdx] - dirichletValues[pvIdx];
        };

        this->asImp_().enforceDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Enforce Dirichlet constraints
     */
    template<typename ApplyFunction>
    void enforceDirichletConstraints(const ApplyFunction& applyDirichlet)
    {
        // enforce Dirichlet boundary conditions
        this->asImp_().evalDirichletBoundaries(applyDirichlet);
        // take care of internal Dirichlet constraints (if enabled)
        this->asImp_().enforceInternalDirichletConstraints(applyDirichlet);
    }

    /*!
     * \brief Evaluates Dirichlet boundaries
     */
    template< typename ApplyDirichletFunctionType >
    void evalDirichletBoundaries(ApplyDirichletFunctionType applyDirichlet)
    {
        // enforce Dirichlet boundaries by overwriting partial derivatives with 1 or 0
        // and set the residual to (privar - dirichletvalue)
        if (this->elemBcTypes().hasDirichlet())
        {
            for (const auto& scvI : scvs(this->fvGeometry()))
            {
                const auto bcTypes = this->elemBcTypes().get(this->fvGeometry(), scvI);
                if (bcTypes.hasDirichlet())
                {
                    const auto dirichletValues = this->asImp_().problem().dirichlet(this->element(), scvI);

                    // set the Dirichlet conditions in residual and jacobian
                    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    {
                        if (bcTypes.isDirichlet(eqIdx))
                        {
                            const auto pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                            assert(0 <= pvIdx && pvIdx < numEq);
                            applyDirichlet(scvI, dirichletValues, eqIdx, pvIdx);
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Update the coupling context for coupled models.
     * \note This does nothing per default (not a coupled model).
     */
    template<class... Args>
    void maybeUpdateCouplingContext(Args&&...) {}

    /*!
     * \brief Update the additional domain derivatives for coupled models.
     * \note This does nothing per default (not a coupled model).
     */
    template<class... Args>
    void maybeEvalAdditionalDomainDerivatives(Args&&...) {}
};

/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (PQ1Bubble methods)
 * \tparam TypeTag The TypeTag
 * \tparam diffMethod The differentiation method to residual compute derivatives
 * \tparam implicit Specifies whether the time discretization is implicit or not not (i.e. explicit)
 */
template<class TypeTag, class Assembler, DiffMethod diffMethod = DiffMethod::numeric, bool implicit = true, class Implementation = void>
class PQ1BubbleLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief Face-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag, class Assembler, class Implementation>
class PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true, Implementation>
: public PQ1BubbleLocalAssemblerBase<TypeTag, Assembler,
                                        Detail::PQ1Bubble::Impl<Implementation, PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>>, true>
{
    using ThisType = PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, true, Implementation>;
    using ParentType = PQ1BubbleLocalAssemblerBase<TypeTag, Assembler, Detail::PQ1Bubble::Impl<Implementation, ThisType>, true>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr auto numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    static constexpr bool enableGridFluxVarsCache = GetPropType<TypeTag, Properties::GridVariables>::GridFluxVariablesCache::cachingEnabled;

public:

    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        // get some aliases for convenience
        const auto& problem = this->asImp_().problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();

        auto&& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to the dofs in the stencil.    //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // if all volvars in the stencil have to be updated or if it's enough to only update the
        // volVars for the scv whose associated dof has been deflected
        static const bool updateAllVolVars = getParamFromGroup<bool>(
            this->asImp_().problem().paramGroup(), "Assembly.PQ1BubbleVolVarsDependOnAllElementDofs", false
        );

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // one residual per element facet
        const auto numElementResiduals = fvGeometry.numScv();

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(numElementResiduals);

        Detail::VolVarsDeflectionHelper deflectionHelper(
            [&] (const auto& scv) -> VolumeVariables& {
                return this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            },
            fvGeometry,
            updateAllVolVars
        );

        // calculation of the derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            deflectionHelper.setCurrent(scv);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalResiduals = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    deflectionHelper.deflect(elemSol, scv, problem);
                    this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);
                    return this->evalLocalResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{problem.paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalResiduals, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (const auto& scvJ : scvs(fvGeometry))
                {
                    // TODO partial reassembly
                    for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                    {
                        // A[i][col][eqIdx][pvIdx] is the rate of change of
                        // the residual of equation 'eqIdx' at dof 'i'
                        // depending on the primary variable 'pvIdx' at dof
                        // 'col'.
                        A[scvJ.dofIndex()][dofIdx][eqIdx][pvIdx] += partialDerivs[scvJ.localDofIndex()][eqIdx];
                    }
                }

                // restore the original state of the scv's volume variables
                deflectionHelper.restore(scv);

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                this->asImp_().maybeUpdateCouplingContext(scv, elemSol, pvIdx);
            }
        }

        return origResiduals;
    }
};


/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler>
class PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/false>
: public PQ1BubbleLocalAssemblerBase<TypeTag, Assembler,
            PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>, false>
{
    using ThisType = PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::numeric, false>;
    using ParentType = PQ1BubbleLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& problem = this->asImp_().problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->asImp_().curSol();
        auto&& curElemVolVars = this->curElemVolVars();

        // get the vecor of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();
        const auto origStorageResiduals = this->evalLocalStorageResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // create the element solution
        auto elemSol = elementSolution(element, curSol, fvGeometry.gridGeometry());

        // create the vector storing the partial derivatives
        ElementResidualVector partialDerivs(fvGeometry.numScv());

        // calculation of the derivatives
        for (auto&& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scv);
            const VolumeVariables origVolVars(curVolVars);

            // calculate derivatives w.r.t to the privars at the dof at hand
            for (int pvIdx = 0; pvIdx < numEq; pvIdx++)
            {
                partialDerivs = 0.0;

                auto evalStorage = [&](Scalar priVar)
                {
                    // auto partialDerivsTmp = partialDerivs;
                    elemSol[scv.localDofIndex()][pvIdx] = priVar;
                    curVolVars.update(elemSol, problem, element, scv);
                    return this->evalLocalStorageResidual();
                };

                // derive the residuals numerically
                static const NumericEpsilon<Scalar, numEq> eps_{problem.paramGroup()};
                static const int numDiffMethod = getParamFromGroup<int>(problem.paramGroup(), "Assembly.NumericDifferenceMethod");
                NumericDifferentiation::partialDerivative(evalStorage, elemSol[scv.localDofIndex()][pvIdx], partialDerivs, origStorageResiduals,
                                                          eps_(elemSol[scv.localDofIndex()][pvIdx], pvIdx), numDiffMethod);

                // update the global stiffness matrix with the current partial derivatives
                for (int eqIdx = 0; eqIdx < numEq; eqIdx++)
                {
                    // A[i][col][eqIdx][pvIdx] is the rate of change of
                    // the residual of equation 'eqIdx' at dof 'i'
                    // depending on the primary variable 'pvIdx' at dof
                    // 'col'.
                    A[dofIdx][dofIdx][eqIdx][pvIdx] += partialDerivs[scv.localDofIndex()][eqIdx];
                }

                // restore the original state of the scv's volume variables
                curVolVars = origVolVars;

                // restore the original element solution
                elemSol[scv.localDofIndex()][pvIdx] = curSol[scv.dofIndex()][pvIdx];
                // TODO additional dof dependencies
            }
        }
        return origResiduals;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler>
class PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/true>
: public PQ1BubbleLocalAssemblerBase<TypeTag, Assembler,
            PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>, true>
{
    using ThisType = PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, true>;
    using ParentType = PQ1BubbleLocalAssemblerBase<TypeTag, Assembler, ThisType, true>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for analytic differentiation");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->asImp_().problem();
        const auto& curElemVolVars = this->curElemVolVars();
        const auto& elemFluxVarsCache = this->elemFluxVarsCache();

        // get the vecor of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            // TODO if e.g. porosity depends on all dofs in the element, we would have off-diagonal matrix entries!?
            if (!this->assembler().isStationaryProblem())
                this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                            problem,
                                                            element,
                                                            fvGeometry,
                                                            volVars,
                                                            scv);

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // add source term derivatives
            this->localResidual().addSourceDerivatives(A[dofIdx][dofIdx],
                                                       problem,
                                                       element,
                                                       fvGeometry,
                                                       volVars,
                                                       scv);
        }

        // localJacobian[scvIdx][otherScvIdx][eqIdx][priVarIdx] of the fluxes
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
            {
                // add flux term derivatives
                this->localResidual().addFluxDerivatives(A,
                                                         problem,
                                                         element,
                                                         fvGeometry,
                                                         curElemVolVars,
                                                         elemFluxVarsCache,
                                                         scvf);
            }

            // the boundary gets special treatment to simplify
            // for the user
            else
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                if (this->elemBcTypes().get(fvGeometry, insideScv).hasNeumann())
                {
                    // add flux term derivatives
                    this->localResidual().addRobinFluxDerivatives(A[insideScv.dofIndex()],
                                                                  problem,
                                                                  element,
                                                                  fvGeometry,
                                                                  curElemVolVars,
                                                                  elemFluxVarsCache,
                                                                  scvf);
                }
            }
        }

        return origResiduals;
    }
};

/*!
 * \ingroup Assembly
 * \ingroup PQ1BubbleDiscretization
 * \brief TODO docme
 */
template<class TypeTag, class Assembler>
class PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, /*implicit=*/false>
: public PQ1BubbleLocalAssemblerBase<TypeTag, Assembler,
            PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>, false>
{
    using ThisType = PQ1BubbleLocalAssembler<TypeTag, Assembler, DiffMethod::analytic, false>;
    using ParentType = PQ1BubbleLocalAssemblerBase<TypeTag, Assembler, ThisType, false>;
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();

public:
    using ParentType::ParentType;

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template <class PartialReassembler = DefaultPartialReassembler>
    ElementResidualVector assembleJacobianAndResidualImpl(JacobianMatrix& A, GridVariables& gridVariables,
                                                          const PartialReassembler* partialReassembler = nullptr)
    {
        if (partialReassembler)
            DUNE_THROW(Dune::NotImplemented, "partial reassembly for explicit time discretization");

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& problem = this->asImp_().problem();
        const auto& curElemVolVars = this->curElemVolVars();

        // get the vector of the actual element residuals
        const auto origResiduals = this->evalLocalResidual();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // calculation of the source and storage derivatives
        for (const auto& scv : scvs(fvGeometry))
        {
            // dof index and corresponding actual pri vars
            const auto dofIdx = scv.dofIndex();
            const auto& volVars = curElemVolVars[scv];

            // derivative of this scv residual w.r.t the d.o.f. of the same scv (because of mass lumping)
            // only if the problem is instationary we add derivative of storage term
            this->localResidual().addStorageDerivatives(A[dofIdx][dofIdx],
                                                        problem,
                                                        element,
                                                        fvGeometry,
                                                        volVars,
                                                        scv);
        }

        return origResiduals;
    }

}; // explicit PQ1BubbleLocalAssembler with analytic Jacobian

} // end namespace Dumux

#endif
