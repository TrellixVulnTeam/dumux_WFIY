// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup StaggeredDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (staggered method)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
#ifndef DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_STAGGERED_LOCAL_ASSEMBLER_HH

#include <dune/common/reservedvector.hh>
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dune/istl/matrixindexset.hh>

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numericdifferentiation.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvlocalassemblerbase.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A base class for all local assemblers
 * \tparam id the id of the sub domain
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation, bool isImplicit = true>
class SubDomainStaggeredLocalAssemblerBase : public FVLocalAssemblerBase<TypeTag, Assembler,Implementation, isImplicit>
{
    using ParentType = FVLocalAssemblerBase<TypeTag, Assembler,Implementation, isImplicit>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using LocalResidualValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using ElementResidualVector = typename LocalResidual::ElementResidualVector;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename Assembler::SolutionVector;
    using SubSolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    using CouplingManager = typename Assembler::CouplingManager;

    static constexpr auto numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();

public:
    static constexpr auto domainId = typename Dune::index_constant<id>();
    static constexpr auto cellCenterId = typename Dune::index_constant<0>();
    static constexpr auto faceId = typename Dune::index_constant<1>();

    static constexpr auto faceOffset = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

    using ParentType::ParentType;

    explicit SubDomainStaggeredLocalAssemblerBase(const Assembler& assembler,
                                                  const Element& element,
                                                  const SolutionVector& curSol,
                                                  CouplingManager& couplingManager)
    : ParentType(assembler,
                 element,
                 curSol,
                 localView(assembler.fvGridGeometry(domainId)),
                 localView(assembler.gridVariables(domainId).curGridVolVars()),
                 localView(assembler.gridVariables(domainId).prevGridVolVars()),
                 localView(assembler.gridVariables(domainId).gridFluxVarsCache()),
                 assembler.localResidual(domainId),
                 (element.partitionType() == Dune::GhostEntity))
    , curElemFaceVars_(localView(assembler.gridVariables(domainId).curGridFaceVars()))
    , prevElemFaceVars_(localView(assembler.gridVariables(domainId).prevGridFaceVars()))
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleJacobianAndResidual(JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        assembleJacobianAndResidualImpl_(domainId, jacRow, res, gridVariables);
    }

    /*!
     * \brief Assemble the residual only
     */
    template<class SubSol>
    void assembleResidual(SubSol& res)
    {
        this->asImp_().bindLocalViews();
        this->elemBcTypes().update(problem(), this->element(), this->fvGeometry());

        assembleResidualImpl_(domainId, res);
    }

    /*!
     * \brief Convenience function to evaluate the complete local residual for the current element. Automatically chooses the the appropriate
     *        element volume variables.
     */
    CellCenterResidualValue evalLocalResidualForCellCenter() const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            return CellCenterResidualValue(0.0);
        }

        return isImplicit ? evalLocalResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current cell center.
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    CellCenterResidualValue evalLocalResidualForCellCenter(const ElementVolumeVariables& elemVolVars,
                                                           const ElementFaceVariables& elemFaceVars) const
    {
        auto residual = evalLocalFluxAndSourceResidualForCellCenter(elemVolVars, elemFaceVars);

        if (!this->assembler().isStationaryProblem())
            residual += evalLocalStorageResidualForCellCenter();

        this->localResidual().evalBoundaryForCellCenter(residual, this->problem(), this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache());

        return residual;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    CellCenterResidualValue evalLocalFluxAndSourceResidualForCellCenter() const
    {
        return isImplicit ? evalLocalFluxAndSourceResidualForCellCenter(this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalFluxAndSourceResidualForCellCenter(this->prevElemVolVars(), this->prevElemFaceVars());
     }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element.
     *
     * \param elemVolVars The element volume variables
     * \param elemVolVars The element face variables
     */
    CellCenterResidualValue evalLocalFluxAndSourceResidualForCellCenter(const ElementVolumeVariables& elemVolVars, const ElementFaceVariables& elemFaceVars) const
    {
        return this->localResidual().evalFluxAndSourceForCellCenter(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache());
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     */
    CellCenterResidualValue evalLocalStorageResidualForCellCenter() const
    {
        return this->localResidual().evalStorageForCellCenter(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars());
    }

    /*!
     * \brief Convenience function to evaluate the  local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf) const
    {
        if (!isImplicit)
            if (this->assembler().isStationaryProblem())
                DUNE_THROW(Dune::InvalidStateException, "Using explicit jacobian assembler with stationary local residual");

        if (this->elementIsGhost())
        {
            return FaceResidualValue(0.0);
        }

        return isImplicit ? evalLocalResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the complete local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    FaceResidualValue evalLocalResidualForFace(const SubControlVolumeFace& scvf,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars) const
    {
        auto residual = evalLocalFluxAndSourceResidualForFace(scvf, elemVolVars, elemFaceVars);

        if (!this->assembler().isStationaryProblem())
            residual += evalLocalStorageResidualForFace(scvf);

        this->localResidual().evalBoundaryForFace(residual, this->problem(), this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf);

        return residual;
    }

    /*!
     * \brief Convenience function to evaluate the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current element. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf) const
    {
        return isImplicit ? evalLocalFluxAndSourceResidualForFace(scvf, this->curElemVolVars(), this->curElemFaceVars())
                          : evalLocalFluxAndSourceResidualForFace(scvf, this->prevElemVolVars(), this->prevElemFaceVars());
    }

    /*!
     * \brief Evaluates the flux and source terms (i.e, the terms without a time derivative)
     *        of the local residual for the current face.
     * \param scvf The sub control volume face
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     */
    FaceResidualValue evalLocalFluxAndSourceResidualForFace(const SubControlVolumeFace& scvf,
                                                            const ElementVolumeVariables& elemVolVars,
                                                            const ElementFaceVariables& elemFaceVars) const
    {
        return this->localResidual().evalFluxAndSourceForFace(this->element(), this->fvGeometry(), elemVolVars, elemFaceVars, this->elemBcTypes(), this->elemFluxVarsCache(), scvf);
    }

    /*!
     * \brief Convenience function to evaluate storage term (i.e, the term with a time derivative)
     *        of the local residual for the current face. Automatically chooses the the appropriate
     *        element volume and face variables.
     * \param scvf The sub control volume face
     */
    FaceResidualValue evalLocalStorageResidualForFace(const SubControlVolumeFace& scvf) const
    {
        return this->localResidual().evalStorageForFace(this->element(), this->fvGeometry(), this->prevElemVolVars(), this->curElemVolVars(), this->prevElemFaceVars(), this->curElemFaceVars(), scvf);
    }

    const Problem& problem() const
    { return this->assembler().problem(domainId); }

    //! The current element volume variables
    ElementFaceVariables& curElemFaceVars()
    { return curElemFaceVars_; }

    //! The element volume variables of the provious time step
    ElementFaceVariables& prevElemFaceVars()
    { return prevElemFaceVars_; }

    //! The current element volume variables
    const ElementFaceVariables& curElemFaceVars() const
    { return curElemFaceVars_; }

    //! The element volume variables of the provious time step
    const ElementFaceVariables& prevElemFaceVars() const
    { return prevElemFaceVars_; }

    CouplingManager& couplingManager()
    { return couplingManager_; }

private:

    //! Assembles the residuals for the cell center dofs.
    template<class SubSol>
    void assembleResidualImpl_(Dune::index_constant<0>, SubSol& res)
    {
        const auto cellCenterGlobalI = problem().fvGridGeometry().elementMapper().index(this->element());
        res[cellCenterGlobalI] = this->asImp_().assembleCellCenterResidualImpl();
    }

    //! Assembles the residuals for the face dofs.
    template<class SubSol>
    void assembleResidualImpl_(Dune::index_constant<1>, SubSol& res)
    {
        for(auto&& scvf : scvfs(this->fvGeometry()))
        {
            res[scvf.dofIndex()] +=  this->asImp_().assembleFaceResidualImpl(scvf);

        }
    }

    //! Assembles the residuals and derivatives for the cell center dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    auto assembleJacobianAndResidualImpl_(Dune::index_constant<0>, JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        const auto cellCenterGlobalI = problem().fvGridGeometry().elementMapper().index(this->element());
        const auto residual = this->asImp_().assembleCellCenterJacobianAndResidualImpl(jacRow[domainId], *std::get<domainId>(gridVariables));
        res[cellCenterGlobalI] = residual;


        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<Dune::Hybrid::size(jacRow), domainId>{};
        forEach(otherDomainIds, [&, domainId = domainId](auto&& domainJ)
        {
            this->asImp_().assembleJacobianCellCenterCoupling(domainJ, jacRow[domainJ], residual, *std::get<domainJ>(gridVariables));
        });
    }

    //! Assembles the residuals and derivatives for the face dofs.
    template<class JacobianMatrixRow, class SubSol, class GridVariablesTuple>
    void assembleJacobianAndResidualImpl_(Dune::index_constant<1>, JacobianMatrixRow& jacRow, SubSol& res, GridVariablesTuple& gridVariables)
    {
        const auto residual = this->asImp_().assembleFaceJacobianAndResidualImpl(jacRow[domainId], *std::get<domainId>(gridVariables));

        for(auto&& scvf : scvfs(this->fvGeometry()))
            res[scvf.dofIndex()] += residual[scvf.localFaceIdx()];

        // for the coupling blocks
        using namespace Dune::Hybrid;
        static constexpr auto otherDomainIds = makeIncompleteIntegerSequence<Dune::Hybrid::size(jacRow), domainId>{};
        forEach(otherDomainIds, [&, domainId = domainId](auto&& domainJ)
        {
            this->asImp_().assembleJacobianFaceCoupling(domainJ, jacRow[domainJ], residual, *std::get<domainJ>(gridVariables));
        });
    }

    ElementFaceVariables curElemFaceVars_;
    ElementFaceVariables prevElemFaceVars_;
    CouplingManager& couplingManager_; //!< the coupling manager
};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief A base class for all implicit local assemblers
 * \tparam TypeTag the TypeTag
 * \tparam Assembler the assembler type
 */
template<std::size_t id, class TypeTag, class Assembler, class Implementation>
class SubDomainStaggeredLocalAssemblerImplicitBase : public SubDomainStaggeredLocalAssemblerBase<id, TypeTag, Assembler, Implementation>
{
    using ParentType = SubDomainStaggeredLocalAssemblerBase<id, TypeTag, Assembler, Implementation>;
    static constexpr auto domainId = Dune::index_constant<id>();
public:
    using ParentType::ParentType;

    void bindLocalViews()
    {
        // get some references for convenience
        auto& couplingManager = this->couplingManager();
        const auto& element = this->element();
        const auto& curSol = this->curSol()[domainId];
        auto&& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        auto&& curElemFaceVars = this->curElemFaceVars();
        auto&& elemFluxVarsCache = this->elemFluxVarsCache();

        // bind the caches
        couplingManager.bindCouplingContext(domainId, element, this->assembler());
        fvGeometry.bind(element);
        curElemVolVars.bind(element, fvGeometry, curSol);
        curElemFaceVars.bind(element, fvGeometry, curSol);
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);
        if (!this->assembler().isStationaryProblem())
        {
            this->prevElemVolVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
            this->prevElemFaceVars().bindElement(element, fvGeometry, this->assembler().prevSol()[domainId]);
        }
    }

};

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (cell-centered methods)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler, DiffMethod DM = DiffMethod::numeric, bool implicit = true>
class SubDomainStaggeredLocalAssembler;

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Cell-centered scheme local assembler using numeric differentiation and implicit time discretization
 */
template<std::size_t id, class TypeTag, class Assembler>
class SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>
: public SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler,
            SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, true> >
{
    using ThisType = SubDomainStaggeredLocalAssembler<id, TypeTag, Assembler, DiffMethod::numeric, /*implicit=*/true>;
    using ParentType = SubDomainStaggeredLocalAssemblerImplicitBase<id, TypeTag, Assembler, ThisType>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using CellCenterResidualValue = typename LocalResidual::CellCenterResidualValue;
    using FaceResidualValue = typename LocalResidual::FaceResidualValue;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables)::LocalView;
    using FaceVariables = typename ElementFaceVariables::FaceVariables;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);

    enum { numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq() };
    enum { dim = GET_PROP_TYPE(TypeTag, GridView)::dimension };

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr int maxNeighbors = 4*(2*dim);
    static constexpr auto domainI = Dune::index_constant<id>();
    static constexpr auto cellCenterId = typename Dune::index_constant<0>();
    static constexpr auto faceId = typename Dune::index_constant<1>();

public:
    using ParentType::ParentType;

    CellCenterResidualValue assembleCellCenterResidualImpl()
    {
        return this->evalLocalResidualForCellCenter();
    }

    FaceResidualValue assembleFaceResidualImpl(const SubControlVolumeFace& scvf)
    {
        return this->evalLocalFluxAndSourceResidualForFace(scvf);
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    CellCenterResidualValue assembleCellCenterJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        assert(domainI == cellCenterId);
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        auto&& curElemVolVars = this->curElemVolVars();
        const auto& fvGridGeometry = this->problem().fvGridGeometry();
        const auto& curSol = this->curSol()[domainI];

        // get the vecor of the acutal element residuals
        // const auto origResiduals = this->evalLocalResidual();

        const auto cellCenterGlobalI = fvGridGeometry.elementMapper().index(element);
        const auto origResidual = this->evalLocalResidualForCellCenter();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

       // build derivatives with for cell center dofs w.r.t. cell center dofs

       const auto& connectivityMap = fvGridGeometry.connectivityMap();

       // create the vector storing the partial derivatives
       CellCenterResidualValue partialDeriv(0.0);

       for(const auto& globalJ : connectivityMap(cellCenterId, cellCenterId, cellCenterGlobalI))
       {
           // get the volVars of the element with respect to which we are going to build the derivative
           auto&& scvJ = fvGeometry.scv(globalJ);
           const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
           auto& curVolVars =  this->getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
           const auto origVolVars(curVolVars);

           // const auto origElemSol = elementSolution<FVElementGeometry>(curSol[globalJ]);

           for(auto pvIdx : priVarIndices_(cellCenterId))
           {
               using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
               const auto& cellCenterPriVars = curSol[globalJ];
               PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

               constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;
               partialDeriv = 0.0;

               auto evalResidual = [&](Scalar priVar)
               {
                   // update the volume variables and compute element residual
                   priVars[pvIdx + offset] = priVar;
                   auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                   curVolVars.update(elemSol, this->problem(), elementJ, scvJ);
                   return this->evalLocalResidualForCellCenter();
               };

               // derive the residuals numerically
               static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
               static const NumericEpsilon<Scalar, numEq> eps{this->problem().paramGroup()};
               NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, origResidual,
                                                         eps(priVars[pvIdx + offset], pvIdx + offset), numDiffMethod);


               // update the global jacobian matrix with the current partial derivatives
               updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

               // restore the original volVars
               curVolVars = origVolVars;
           }
       }

        return origResidual;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianMatrixDiagBlock, class GridVariables>
    auto assembleFaceJacobianAndResidualImpl(JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {
        assert(domainI == faceId);
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->problem().fvGridGeometry();
        const auto& curSol = this->curSol()[domainI];

        using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector); // TODO: use reserved vector
        FaceSolutionVector origResiduals;
        origResiduals.resize(fvGeometry.numScvf());
        origResiduals = 0.0;

        // get the vecor of the acutal element residuals
        // const auto origResiduals = this->evalLocalResidual();
        // treat the local residua of the face dofs:
        for(auto&& scvf : scvfs(fvGeometry))
            origResiduals[scvf.localFaceIdx()] =  this->evalLocalResidualForFace(scvf);

        //////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                              //
        // Calculate derivatives of all dofs in stencil with respect to the dofs in the element. In the //
        // neighboring elements we do so by computing the derivatives of the fluxes which depend on the //
        // actual element. In the actual element we evaluate the derivative of the entire residual.     //
        //                                                                                              //
        //////////////////////////////////////////////////////////////////////////////////////////////////

       // build derivatives with for cell center dofs w.r.t. cell center dofs

       const auto& connectivityMap = fvGridGeometry.connectivityMap();

       FaceResidualValue partialDeriv(0.0);

       for(auto&& scvf : scvfs(fvGeometry))
       {
           // set the actual dof index
           const auto faceGlobalI = scvf.dofIndex();

           using FaceSolution = typename GET_PROP_TYPE(TypeTag, StaggeredFaceSolution);
           const auto origFaceSolution = FaceSolution(scvf, curSol, fvGridGeometry);

           // build derivatives with for face dofs w.r.t. cell center dofs
           for(const auto& globalJ : connectivityMap(faceId, faceId, scvf.index()))
           {
               // get the faceVars of the face with respect to which we are going to build the derivative
            auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), this->curElemFaceVars(), scvf);
            const auto origFaceVars(faceVars);

               for(auto pvIdx : priVarIndices_(faceId))
               {
                   auto faceSolution = origFaceSolution;
                   partialDeriv = 0;

                   auto evalResidual = [&](Scalar priVar)
                   {
                       // update the volume variables and compute element residual
                       faceSolution[globalJ][pvIdx] = priVar;
                       faceVars.update(faceSolution, problem, element, fvGeometry, scvf);
                       return this->evalLocalResidualForFace(scvf);
                   };

                   // derive the residuals numerically
                   static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                   static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                   NumericDifferentiation::partialDerivative(evalResidual, faceSolution[globalJ][pvIdx], partialDeriv, origResiduals[scvf.localFaceIdx()], eps_(faceSolution[globalJ][pvIdx], pvIdx), numDiffMethod);

                   // update the global jacobian matrix with the current partial derivatives
                   updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx - ParentType::faceOffset, partialDeriv);

                   // restore the original faceVars
                   faceVars = origFaceVars;
               }
           }
       }

        return origResiduals;
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class GridVariables>
    void assembleJacobianCellCenterCoupling(Dune::index_constant<1> domainJ, JacobianBlock& A,
                                            const CellCenterResidualValue& origResidual, GridVariables& gridVariables)
    {
        // ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        // ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->problem().fvGridGeometry();
        // const auto& curSol = this->curSol()[domainI];
        // build derivatives with for cell center dofs w.r.t. cell center dofs
        const auto cellCenterGlobalI = fvGridGeometry.elementMapper().index(element);

        // create the vector storing the partial derivatives
        CellCenterResidualValue partialDeriv(0.0);

        for(const auto& scvfJ : scvfs(fvGeometry))
        {
             const auto globalJ = scvfJ.dofIndex();

            // get the faceVars of the face with respect to which we are going to build the derivative
            auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), this->curElemFaceVars(), scvfJ);
            const auto origFaceVars(faceVars);

            for(auto pvIdx : priVarIndices_(faceId))
            {
                auto facePriVars(this->curSol()[faceId][globalJ]);
                partialDeriv = 0.0;

                auto evalResidual = [&](Scalar priVar)
                {
                    // update the volume variables and compute element residual
                    facePriVars[pvIdx] = priVar;
                    faceVars.updateOwnFaceOnly(facePriVars);
                    return this->evalLocalResidualForCellCenter();
                };

                // derive the residuals numerically
                static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                NumericDifferentiation::partialDerivative(evalResidual, facePriVars[pvIdx], partialDeriv, origResidual, eps_(facePriVars[pvIdx], pvIdx), numDiffMethod);

                // update the global jacobian matrix with the current partial derivatives
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx - ParentType::faceOffset, partialDeriv);

                // restore the original faceVars
                faceVars = origFaceVars;
            }
        }

    }

    template<std::size_t otherId, class JacobianBlock, class GridVariables>
    void assembleJacobianCellCenterCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                            const CellCenterResidualValue& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& element = this->element();

        // get stencil informations
        const auto& stencil = this->couplingManager().couplingStencil(domainI, element, domainJ);

        if(stencil.empty())
            return;

        for (const auto globalJ : stencil)
        {
            const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
            const auto& curSol = this->curSol()[domainJ];
            const auto origPriVarsJ = curSol[globalJ];

            for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
            {
                auto evalCouplingResidual = [&](Scalar priVar)
                {
                    auto deflectedPriVarsJ = origPriVarsJ;
                    deflectedPriVarsJ[pvIdx] = priVar;
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, deflectedPriVarsJ, pvIdx);
                    return this->couplingManager().evalCouplingResidual(domainI, *this, domainJ, globalJ);
                };

                // derive the residuals numerically
                CellCenterResidualValue partialDeriv(0.0);
                NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDeriv, origResidual);

                // update the global stiffness matrix with the current partial derivatives
                const auto cellCenterGlobalI = this->problem().fvGridGeometry().elementMapper().index(element);
                updateGlobalJacobian_(A, cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

                // restore the undeflected state of the coupling context
                this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, origPriVarsJ, pvIdx);
            }
        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix.
     *
     * \return The element residual at the current solution.
     */
    template<class JacobianBlock, class ElementResidualVector, class GridVariables>
    void assembleJacobianFaceCoupling(Dune::index_constant<0> domainJ, JacobianBlock& A,
                                      const ElementResidualVector& origResiduals, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& problem = this->problem();
        const auto& fvGeometry = this->fvGeometry();
        const auto& fvGridGeometry = this->problem().fvGridGeometry();
        const auto& connectivityMap = fvGridGeometry.connectivityMap();
        FaceResidualValue partialDeriv;

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        for(auto&& scvf : scvfs(fvGeometry))
        {

            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : connectivityMap(faceId, cellCenterId, scvf.index()))
            {
                // get the volVars of the element with respect to which we are going to build the derivative
                auto&& scvJ = fvGeometry.scv(globalJ);
                const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
                auto& curVolVars = this->getVolVarAccess(gridVariables.curGridVolVars(), this->curElemVolVars(), scvJ);
                const auto origVolVars(curVolVars);
                const auto origCellCenterPriVars(this->curSol()[cellCenterId][globalJ]);

                for(auto pvIdx : priVarIndices_(cellCenterId))
                {

                    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
                    const auto& cellCenterPriVars = this->curSol()[cellCenterId][globalJ];
                    PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);

                    constexpr auto offset = PrimaryVariables::dimension - CellCenterPrimaryVariables::dimension;

                    partialDeriv = 0;

                    auto evalResidual = [&](Scalar priVar)
                    {
                        // update the volume variables and compute element residual
                        priVars[pvIdx + offset] = priVar;
                        auto elemSol = elementSolution<FVElementGeometry>(std::move(priVars));
                        curVolVars.update(elemSol, problem, elementJ, scvJ);
                        return this->evalLocalResidualForFace(scvf);
                    };

                    // derive the residuals numerically
                    static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                    static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                    NumericDifferentiation::partialDerivative(evalResidual, priVars[pvIdx + offset], partialDeriv, origResiduals[scvf.localFaceIdx()], eps_(priVars[pvIdx + offset], pvIdx + offset), numDiffMethod);

                    // update the global jacobian matrix with the current partial derivatives
                    updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the original volVars
                    curVolVars = origVolVars;
                }
            }
        }

    }

    template<std::size_t otherId, class JacobianBlock, class ElementResidualVector, class GridVariables>
    void assembleJacobianFaceCoupling(Dune::index_constant<otherId> domainJ, JacobianBlock& A,
                                      const ElementResidualVector& res, GridVariables& gridVariables)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate derivatives of all dofs in the element with respect to all dofs in the coupling stencil. //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        // get some aliases for convenience
        const auto& fvGeometry = this->fvGeometry();
        const auto& curSol = this->curSol()[domainJ];
        FaceResidualValue partialDeriv;

        // build derivatives with for cell center dofs w.r.t. cell center dofs
        for(auto&& scvf : scvfs(fvGeometry))
        {

            // set the actual dof index
            const auto faceGlobalI = scvf.dofIndex();

            // get stencil informations
            const auto& stencil = this->couplingManager().couplingStencil(domainI, scvf, domainJ);

            if(stencil.empty())
                continue;

            // build derivatives with for face dofs w.r.t. cell center dofs
            for(const auto& globalJ : stencil)
            {
                const auto origPriVarsJ = curSol[globalJ];
                const auto origResidual = this->couplingManager().evalCouplingResidual(domainI, scvf, *this, domainJ, globalJ);

                for (int pvIdx = 0; pvIdx < JacobianBlock::block_type::cols; ++pvIdx)
                {
                    auto evalCouplingResidual = [&](Scalar priVar)
                    {
                        auto deflectedPriVars = origPriVarsJ;
                        deflectedPriVars[pvIdx] = priVar;
                        this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, deflectedPriVars, pvIdx);
                        return this->couplingManager().evalCouplingResidual(domainI, scvf, *this, domainJ, globalJ);
                    };

                    // derive the residuals numerically
                    FaceResidualValue partialDeriv(0.0);
                    static const int numDiffMethod = getParamFromGroup<int>(this->problem().paramGroup(), "Assembly.NumericDifferenceMethod");
                    static const NumericEpsilon<Scalar, numEq> eps_{this->problem().paramGroup()};
                    NumericDifferentiation::partialDerivative(evalCouplingResidual, origPriVarsJ[pvIdx], partialDeriv, origResidual, eps_(origPriVarsJ[pvIdx], pvIdx), numDiffMethod);

                    // update the global stiffness matrix with the current partial derivatives
                    updateGlobalJacobian_(A, faceGlobalI, globalJ, pvIdx, partialDeriv);

                    // restore the undeflected state of the coupling context
                    this->couplingManager().updateCouplingContext(domainI, *this, domainJ, globalJ, origPriVarsJ, pvIdx);
                }
            }
        }


    }

    template<class JacobianMatrixDiagBlock, class GridVariables>
    void evalAdditionalDerivatives(const std::vector<std::size_t>& additionalDofDependencies,
                                   JacobianMatrixDiagBlock& A, GridVariables& gridVariables)
    {

    }

    /*!
     * \brief Updates the current global Jacobian matrix with the
     *        partial derivatives of all equations in regard to the
     *        primary variable 'pvIdx' at dof 'col'. Specialization for cc methods.
     */
    template<class SubMatrix, class CCOrFacePrimaryVariables>
    static void updateGlobalJacobian_(SubMatrix& matrix,
                                      const int globalI,
                                      const int globalJ,
                                      const int pvIdx,
                                      const CCOrFacePrimaryVariables &partialDeriv)
    {
        for (int eqIdx = 0; eqIdx < partialDeriv.size(); eqIdx++)
        {
            // A[i][col][eqIdx][pvIdx] is the rate of change of
            // the residual of equation 'eqIdx' at dof 'i'
            // depending on the primary variable 'pvIdx' at dof
            // 'col'.

            assert(pvIdx >= 0);
            assert(eqIdx < matrix[globalI][globalJ].size());
            assert(pvIdx < matrix[globalI][globalJ][eqIdx].size());
            matrix[globalI][globalJ][eqIdx][pvIdx] += partialDeriv[eqIdx];
        }
    }

private:
    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for cell center dofs.
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::CellCenterIdx)
    {
        constexpr auto numEqCellCenter =  GET_PROP_VALUE(TypeTag, NumEqCellCenter);

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(0, numEqCellCenter);
#else
        return IntRange(0, numEqCellCenter);
#endif
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for face dofs.
    static auto priVarIndices_(typename FVGridGeometry::DofTypeIndices::FaceIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
    #if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(numEqCellCenter, numEqCellCenter + numEqFace);
    #else
        return IntRange(numEqCellCenter, numEqCellCenter + numEqFace);
    #endif
    }

    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& gridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return elemFaceVars[scvf]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& gridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return gridFaceVariables.faceVars(scvf.index()); }
};



} // end namespace Dumux

#endif
