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
 * \ingroup StaggeredDiscretization
 * \ingroup Assembly
 * \brief An assembler for Jacobian and residual contribution per element (staggered FV method)
 */
#ifndef DUMUX_CC_LOCAL_ASSEMBLER_HH
#define DUMUX_CC_LOCAL_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/staggered/facesolution.hh>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#else
#include <dumux/common/intrange.hh>
#endif

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief An assembler for Jacobian and residual contribution per element (staggered FV method)
 * \tparam TypeTag the TypeTag
 * \tparam DM the differentiation method to residual compute derivatives
 * \tparam implicit if to use an implicit or explicit time discretization
 */
template<class TypeTag,
         DiffMethod diffMethod = DiffMethod::numeric,
         bool implicit = true>
class StaggeredLocalAssembler;

/*!
 * \ingroup Assembly
 * \brief Staggerd FV scheme local assembler using numeric differentiation and implicit time discretization
 */
template<class TypeTag>
class StaggeredLocalAssembler<TypeTag,
                       DiffMethod::numeric,
                       /*implicit=*/true>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using GridFaceVariables = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    using NumCellCenterEqVector =  typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using NumFaceEqVector = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using FaceSolutionVector = typename GET_PROP_TYPE(TypeTag, FaceSolutionVector);
    using FaceSolution = typename GET_PROP_TYPE(TypeTag, StaggeredFaceSolution);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using ElementFaceVariables = typename GET_PROP_TYPE(TypeTag, ElementFaceVariables);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr bool enableGridFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    static constexpr auto faceOffset = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

public:

    /*!
     * \brief Assemble the residual only
     */
    template<class Assembler>
    static void assembleResidual(Assembler& assembler,
                                SolutionVector& res,
                                const Element& element,
                                const SolutionVector& curSol)
    {
        using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
        typename DofTypeIndices::CellCenterIdx cellCenterIdx;
        typename DofTypeIndices::FaceIdx faceIdx;

        // get some references for convenience
        const auto& problem = assembler.problem();
        auto localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto curElemFaceVars = localView(gridVariables.curGridFaceVars());
        curElemFaceVars.bind(element, fvGeometry, curSol);

        const bool isStationary = localResidual.isStationary();
        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        auto prevElemFaceVars = localView(gridVariables.prevGridFaceVars());
        if (!isStationary)
        {
            prevElemVolVars.bindElement(element, fvGeometry, assembler.prevSol());
            prevElemFaceVars.bindElement(element, fvGeometry, assembler.prevSol());
        }

        // for compatibility with box models
        ElementBoundaryTypes elemBcTypes;

        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[cellCenterIdx][cellCenterGlobalI] = localResidual.evalCellCenter(problem,
                                                                             element,
                                                                             fvGeometry,
                                                                             prevElemVolVars,
                                                                             curElemVolVars,
                                                                             prevElemFaceVars,
                                                                             curElemFaceVars,
                                                                             elemBcTypes,
                                                                             elemFluxVarsCache);

        // treat the local residua of the face dofs:
        for(auto&& scvf : scvfs(fvGeometry))
        {
            res[faceIdx][scvf.dofIndex()] += localResidual.evalFace(problem,
                                                                    element,
                                                                    fvGeometry,
                                                                    scvf,
                                                                    prevElemVolVars,
                                                                    curElemVolVars,
                                                                    prevElemFaceVars,
                                                                    curElemFaceVars,
                                                                    elemBcTypes,
                                                                    elemFluxVarsCache);

        }
    }

    /*!
     * \brief Computes the derivatives with respect to the given element and adds them
     *        to the global matrix. The element residual is written into the right hand side.
     */
    template<class Assembler>
    static void assembleJacobianAndResidual(Assembler& assembler,
                                            JacobianMatrix& jac,
                                            SolutionVector& res,
                                            const Element& element,
                                            const SolutionVector& curSol)
    {
        using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
        typename DofTypeIndices::CellCenterIdx cellCenterIdx;
        typename DofTypeIndices::FaceIdx faceIdx;


        // get some references for convenience
        const auto& problem = assembler.problem();
        auto localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

        // prepare the local views
        auto fvGeometry = localView(assembler.fvGridGeometry());
        fvGeometry.bind(element);

        auto curElemVolVars = localView(gridVariables.curGridVolVars());
        curElemVolVars.bind(element, fvGeometry, curSol);

        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, curElemVolVars);

        auto curElemFaceVars = localView(gridVariables.curGridFaceVars());
        curElemFaceVars.bind(element, fvGeometry, curSol);

        const bool isStationary = assembler.isStationaryProblem();
        auto prevElemVolVars = localView(gridVariables.prevGridVolVars());
        auto prevElemFaceVars = localView(gridVariables.prevGridFaceVars());
        if (!isStationary)
        {
            prevElemVolVars.bindElement(element, fvGeometry, assembler.prevSol());
            prevElemFaceVars.bindElement(element, fvGeometry, assembler.prevSol());
        }

        // for compatibility with box models
        ElementBoundaryTypes elemBcTypes;

        const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);
        res[cellCenterIdx][cellCenterGlobalI] = localResidual.evalCellCenter(problem,
                                                                             element,
                                                                             fvGeometry,
                                                                             prevElemVolVars,
                                                                             curElemVolVars,
                                                                             prevElemFaceVars,
                                                                             curElemFaceVars,
                                                                             elemBcTypes,
                                                                             elemFluxVarsCache);

        // treat the local residua of the face dofs:
        // create a cache to reuse some results for the calculation of the derivatives
        FaceSolutionVector faceResidualCache;
        faceResidualCache.resize(fvGeometry.numScvf());
        faceResidualCache = 0.0;

        for(auto&& scvf : scvfs(fvGeometry))
        {
            faceResidualCache[scvf.localFaceIdx()] = localResidual.evalFace(problem,
                                                                            element,
                                                                            fvGeometry,
                                                                            scvf,
                                                                            prevElemVolVars,
                                                                            curElemVolVars,
                                                                            prevElemFaceVars,
                                                                            curElemFaceVars,
                                                                            elemBcTypes,
                                                                            elemFluxVarsCache);

            res[faceIdx][scvf.dofIndex()] += faceResidualCache[scvf.localFaceIdx()] ;
        }

        // calculate derivatives of all dofs in stencil with respect to the dofs in the element
        evalPartialDerivatives_(assembler,
                                element,
                                fvGeometry,
                                curSol,
                                prevElemVolVars,
                                curElemVolVars,
                                prevElemFaceVars,
                                curElemFaceVars,
                                elemFluxVarsCache,
                                elemBcTypes,
                                jac,
                                res[cellCenterIdx][cellCenterGlobalI],
                                faceResidualCache);

    }

protected:
    template<class Assembler>
    static void evalPartialDerivatives_(Assembler& assembler,
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SolutionVector& curSol,
                                        const ElementVolumeVariables& prevElemVolVars,
                                        ElementVolumeVariables& curElemVolVars,
                                        const ElementFaceVariables& prevElemFaceVars,
                                        ElementFaceVariables& curElemFaceVars,
                                        ElementFluxVariablesCache& elemFluxVarsCache,
                                        const ElementBoundaryTypes& elemBcTypes,
                                        JacobianMatrix& matrix,
                                        const NumCellCenterEqVector& ccResidual,
                                        const FaceSolutionVector& faceResidualCache)
    {
        // compute the derivatives of the cell center residuals with respect to cell center dofs
        dCCdCC_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the cell center residuals with respect to face dofs
        dCCdFace_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, ccResidual);

        // compute the derivatives of the face residuals with respect to cell center dofs
        dFacedCC_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);

        // compute the derivatives of the face residuals with respect to face dofs
        dFacedFace_(assembler, element, fvGeometry, curSol, prevElemVolVars, curElemVolVars, prevElemFaceVars, curElemFaceVars, elemFluxVarsCache, elemBcTypes, matrix, faceResidualCache);
    }

    /*!
     * \brief Computes the derivatives of the cell center residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dCCdCC_(Assembler& assembler,
                 const Element& element,
                 const FVElementGeometry& fvGeometry,
                 const SolutionVector& curSol,
                 const ElementVolumeVariables& prevElemVolVars,
                 ElementVolumeVariables& curElemVolVars,
                 const ElementFaceVariables& prevElemFaceVars,
                 const ElementFaceVariables& curElemFaceVars,
                 ElementFluxVariablesCache& elemFluxVarsCache,
                 const ElementBoundaryTypes& elemBcTypes,
                 JacobianMatrix& matrix,
                 const NumCellCenterEqVector& ccResidual)
    {
        using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
        typename DofTypeIndices::CellCenterIdx cellCenterIdx;

        const auto& problem = assembler.problem();
        auto localResidual = assembler.localResidual();
        auto& gridVariables = assembler.gridVariables();

       // build derivatives with for cell center dofs w.r.t. cell center dofs
       const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);

       const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();

       for(const auto& globalJ : connectivityMap(cellCenterIdx, cellCenterIdx, cellCenterGlobalI))
       {
           // get the volVars of the element with respect to which we are going to build the derivative
           auto&& scvJ = fvGeometry.scv(globalJ);
           const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
           auto& curVolVars =  getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
           VolumeVariables origVolVars(curVolVars);

           for(auto pvIdx : priVarIndices_(cellCenterIdx))
           {
               CellCenterPrimaryVariables priVars(curSol[cellCenterIdx][globalJ]);

               const Scalar eps = numericEpsilon(priVars[pvIdx], cellCenterIdx, cellCenterIdx);
               priVars[pvIdx] += eps;
               ElementSolutionVector elemSol{std::move(priVars)};
               curVolVars.update(elemSol, problem, elementJ, scvJ);

               const auto deflectedResidual = localResidual.evalCellCenter(problem, element, fvGeometry, prevElemVolVars, curElemVolVars,
                                              prevElemFaceVars, curElemFaceVars,
                                              elemBcTypes, elemFluxVarsCache);

               auto partialDeriv = (deflectedResidual - ccResidual);
               partialDeriv /= eps;

               // update the global jacobian matrix with the current partial derivatives
               updateGlobalJacobian_(matrix[cellCenterIdx][cellCenterIdx], cellCenterGlobalI, globalJ, pvIdx, partialDeriv);

               // restore the original volVars
               curVolVars = origVolVars;
           }
       }
    }

    /*!
     * \brief Computes the derivatives of the cell center residuals with respect to face dofs
     */
    template<class Assembler>
    static void dCCdFace_(Assembler& assembler,
                          const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SolutionVector& curSol,
                          const ElementVolumeVariables& prevElemVolVars,
                          const ElementVolumeVariables& curElemVolVars,
                          const ElementFaceVariables& prevElemFaceVars,
                          ElementFaceVariables& curElemFaceVars,
                          ElementFluxVariablesCache& elemFluxVarsCache,
                          const ElementBoundaryTypes& elemBcTypes,
                          JacobianMatrix& matrix,
                          const NumCellCenterEqVector& ccResidual)
    {
       // build derivatives with for cell center dofs w.r.t. face dofs
       using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
       typename DofTypeIndices::CellCenterIdx cellCenterIdx;
       typename DofTypeIndices::FaceIdx faceIdx;

       const auto& problem = assembler.problem();
       auto localResidual = assembler.localResidual();
       auto& gridVariables = assembler.gridVariables();

       // build derivatives with for cell center dofs w.r.t. cell center dofs
       const auto cellCenterGlobalI = assembler.fvGridGeometry().elementMapper().index(element);

       for(const auto& scvfJ : scvfs(fvGeometry))
       {
            const auto globalJ = scvfJ.dofIndex();

           // get the faceVars of the face with respect to which we are going to build the derivative
           auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), curElemFaceVars, scvfJ);
           const auto origFaceVars(faceVars);

           for(auto pvIdx : priVarIndices_(faceIdx))
           {
               FacePrimaryVariables facePriVars(curSol[faceIdx][globalJ]);
               const Scalar eps = numericEpsilon(facePriVars[pvIdx], cellCenterIdx, faceIdx);
               facePriVars[pvIdx] += eps;

               faceVars.updateOwnFaceOnly(facePriVars);

               const auto deflectedResidual = localResidual.evalCellCenter(problem, element, fvGeometry,
                                              prevElemVolVars, curElemVolVars,
                                              prevElemFaceVars, curElemFaceVars,
                                              elemBcTypes, elemFluxVarsCache);

               auto partialDeriv = (deflectedResidual - ccResidual);
               partialDeriv /= eps;

               // update the global jacobian matrix with the current partial derivatives
               updateGlobalJacobian_(matrix[cellCenterIdx][faceIdx], cellCenterGlobalI, globalJ, pvIdx - faceOffset, partialDeriv);

               // restore the original faceVars
               faceVars = origFaceVars;
           }
       }
    }

    /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dFacedCC_(Assembler& assembler,
                          const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const SolutionVector& curSol,
                          const ElementVolumeVariables& prevElemVolVars,
                          ElementVolumeVariables& curElemVolVars,
                          const ElementFaceVariables& prevElemFaceVars,
                          const ElementFaceVariables& curElemFaceVars,
                          ElementFluxVariablesCache& elemFluxVarsCache,
                          const ElementBoundaryTypes& elemBcTypes,
                          JacobianMatrix& matrix,
                          const FaceSolutionVector& cachedResidual)
    {
       for(auto&& scvf : scvfs(fvGeometry))
       {
           using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
           typename DofTypeIndices::CellCenterIdx cellCenterIdx;
           typename DofTypeIndices::FaceIdx faceIdx;

           const auto& problem = assembler.problem();
           auto localResidual = assembler.localResidual();
           auto& gridVariables = assembler.gridVariables();
           const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();

           // set the actual dof index
           const auto faceGlobalI = scvf.dofIndex();

           // build derivatives with for face dofs w.r.t. cell center dofs
           for(const auto& globalJ : connectivityMap(faceIdx, cellCenterIdx, scvf.index()))
           {
               // get the volVars of the element with respect to which we are going to build the derivative
               auto&& scvJ = fvGeometry.scv(globalJ);
               const auto elementJ = fvGeometry.fvGridGeometry().element(globalJ);
               auto& curVolVars = getVolVarAccess(gridVariables.curGridVolVars(), curElemVolVars, scvJ);
               VolumeVariables origVolVars(curVolVars);

               for(auto pvIdx : priVarIndices_(cellCenterIdx))
               {
                   CellCenterPrimaryVariables priVars(curSol[cellCenterIdx][globalJ]);

                   const Scalar eps = numericEpsilon(priVars[pvIdx], faceIdx, cellCenterIdx);
                   priVars[pvIdx] += eps;
                   ElementSolutionVector elemSol{std::move(priVars)};
                   curVolVars.update(elemSol, problem, elementJ, scvJ);

                   const auto deflectedResidual = localResidual.evalFace(problem, element, fvGeometry, scvf,
                                                  prevElemVolVars, curElemVolVars,
                                                  prevElemFaceVars, curElemFaceVars,
                                                  elemBcTypes, elemFluxVarsCache);

                   auto partialDeriv = (deflectedResidual - cachedResidual[scvf.localFaceIdx()]);
                   partialDeriv /= eps;

                   // update the global jacobian matrix with the current partial derivatives
                   updateGlobalJacobian_(matrix[faceIdx][cellCenterIdx], faceGlobalI, globalJ, pvIdx, partialDeriv);

                   // restore the original volVars
                   curVolVars = origVolVars;
               }
           }
       }
    }

    /*!
     * \brief Computes the derivatives of the face residuals with respect to cell center dofs
     */
    template<class Assembler>
    static void dFacedFace_(Assembler& assembler,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const SolutionVector& curSol,
                            const ElementVolumeVariables& prevElemVolVars,
                            const ElementVolumeVariables& curElemVolVars,
                            const ElementFaceVariables& prevElemFaceVars,
                            ElementFaceVariables& curElemFaceVars,
                            ElementFluxVariablesCache& elemFluxVarsCache,
                            const ElementBoundaryTypes& elemBcTypes,
                            JacobianMatrix& matrix,
                            const FaceSolutionVector& cachedResidual)
    {
        using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
        typename DofTypeIndices::FaceIdx faceIdx;

        const auto& problem = assembler.problem();
        auto localResidual = assembler.localResidual();
        const auto& connectivityMap = assembler.fvGridGeometry().connectivityMap();
        auto& gridVariables = assembler.gridVariables();

       for(auto&& scvf : scvfs(fvGeometry))
       {
           // set the actual dof index
           const auto faceGlobalI = scvf.dofIndex();

           // build derivatives with for face dofs w.r.t. cell center dofs
           for(const auto& globalJ : connectivityMap(faceIdx, faceIdx, scvf.index()))
           {
               // get the faceVars of the face with respect to which we are going to build the derivative
            auto& faceVars = getFaceVarAccess(gridVariables.curGridFaceVars(), curElemFaceVars, scvf);
            const auto origFaceVars(faceVars);

               for(auto pvIdx : priVarIndices_(faceIdx))
               {
                   auto faceSolution = FaceSolution(scvf, curSol[faceIdx], assembler.fvGridGeometry());

                   const Scalar eps = numericEpsilon(faceSolution[globalJ][pvIdx], faceIdx, faceIdx);

                   faceSolution[globalJ][pvIdx] += eps;
                   faceVars.update(faceSolution, problem, element, fvGeometry, scvf);

                   const auto deflectedResidual = localResidual.evalFace(problem, element, fvGeometry, scvf,
                                                  prevElemVolVars, curElemVolVars,
                                                  prevElemFaceVars, curElemFaceVars,
                                                  elemBcTypes, elemFluxVarsCache);

                   auto partialDeriv = (deflectedResidual - cachedResidual[scvf.localFaceIdx()]);
                   partialDeriv /= eps;

                   // update the global jacobian matrix with the current partial derivatives
                   updateGlobalJacobian_(matrix[faceIdx][faceIdx], faceGlobalI, globalJ, pvIdx - faceOffset, partialDeriv);

                   // restore the original faceVars
                   faceVars = origFaceVars;
               }
           }
       }
    }


    /*!
     * \brief Computes the epsilon used for numeric differentiation
     *        for a given value of a primary variable for
     *        an equation residual with respect to a primary variable
     *
     * \param priVar The value of the primary variable
     * \param eqIdx The equation index of the equation
     * \param pvIdx The index of the primary variable
     */
    static Scalar numericEpsilon(const Scalar priVar, const int eqIdx, const int pvIdx)
    {
        // define the base epsilon as the geometric mean of 1 and the
        // resolution of the scalar type. E.g. for standard 64 bit
        // floating point values, the resolution is about 10^-16 and
        // the base epsilon is thus approximately 10^-8.
        /*
        static const Scalar baseEps
            = Dumux::geometricMean<Scalar>(std::numeric_limits<Scalar>::epsilon(), 1.0);
        */
        using BaseEpsilon = typename GET_PROP(TypeTag, BaseEpsilon);
        const std::array<std::array<Scalar, 2>, 2> baseEps_ = BaseEpsilon::getEps();


        static const Scalar baseEps = baseEps_[eqIdx][pvIdx];
        assert(std::numeric_limits<Scalar>::epsilon()*1e4 < baseEps);
        // the epsilon value used for the numeric differentiation is
        // now scaled by the absolute value of the primary variable...
        return baseEps*(std::abs(priVar) + 1.0);
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

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for cell center dofs.
    static auto priVarIndices_(typename GET_PROP(TypeTag, DofTypeIndices)::CellCenterIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(0, numEqCellCenter);
#else
        return IntRange(0, numEqCellCenter);
#endif
    }

    //! Helper function that returns an iterable range of primary variable indices.
    //! Specialization for face dofs.
    static auto priVarIndices_(typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx)
    {
        constexpr auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        constexpr auto numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
    #if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return Dune::range(numEqCellCenter, numEqCellCenter + numEqFace);
    #else
        return IntRange(numEqCellCenter, numEqCellCenter + numEqFace);
    #endif
    }

private:
    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return elemVolVars[scv]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridVolumeVariablesCache), VolumeVariables&>::type
    getVolVarAccess(GridVolumeVariables& gridVolVars, ElementVolumeVariables& elemVolVars, const SubControlVolume& scv)
    { return gridVolVars.volVars(scv); }

    template<class T = TypeTag>
    static typename std::enable_if<!GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& GridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return elemFaceVars[scvf]; }

    template<class T = TypeTag>
    static typename std::enable_if<GET_PROP_VALUE(T, EnableGridFaceVariablesCache), FaceVariables&>::type
    getFaceVarAccess(GridFaceVariables& GridFaceVariables, ElementFaceVariables& elemFaceVars, const SubControlVolumeFace& scvf)
    { return GridFaceVariables.faceVars(scvf.index()); }
};

} // end namespace Dumux

#endif
