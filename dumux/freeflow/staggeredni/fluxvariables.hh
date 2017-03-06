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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_FREELOW_IMPLICIT_NI_FLUXVARIABLES_HH
#define DUMUX_FREELOW_IMPLICIT_NI_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include "../staggered/fluxvariables.hh"

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
}

// // forward declaration
// template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
// class FreeFlowFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
// template<class TypeTag>
// using FreeFlowFluxVariables = FreeFlowFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport),
//                                                                  GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

// specialization for immiscible, non-isothermal flow
template<class TypeTag>
class FreeFlowFluxVariablesImpl<TypeTag, false, true> : public FreeFlowFluxVariablesImpl<TypeTag, false, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables); // TODO ?
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

//    using MolecularDiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);

    static constexpr bool navierStokes = GET_PROP_VALUE(TypeTag, EnableInertiaTerms);
//    static constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

//    //! The index of the component balance equation that gets replaced with the total mass balance
//    static const int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);

    using ParentType = FreeFlowFluxVariablesImpl<TypeTag, false, false>;

    enum {
         // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        pressureIdx = Indices::pressureIdx,
        velocityIdx = Indices::velocityIdx,
        temperatureIdx = Indices::temperatureIdx, // TODO necessary?

        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        energyBalanceIdx = Indices::energyBalanceIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        phaseIdx = Indices::phaseIdx
    };

public:
    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                                        const Element &element,
                                                        const FVElementGeometry& fvGeometry,
                                                        const ElementVolumeVariables& elemVolVars,
                                                        const GlobalFaceVars& globalFaceVars,
                                                        const SubControlVolumeFace &scvf,
                                                        const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux(0.0);

        flux += advectiveFluxForCellCenter_(problem, fvGeometry, elemVolVars, globalFaceVars, scvf);
        flux += diffusiveFluxForCellCenter_(problem, fvGeometry, elemVolVars, scvf);
        return flux;
    }

private:

    CellCenterPrimaryVariables advectiveFluxForCellCenter_(const Problem& problem,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const GlobalFaceVars& globalFaceVars,
                                                          const SubControlVolumeFace &scvf)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        // TODO: catch neumann and outflow in localResidual's evalBoundary_()
        bool isOutflow = false;
        if(scvf.boundary())
        {
            const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isOutflow(momentumBalanceIdx))
                    isOutflow = true;
        }

        const auto& outsideVolVars = isOutflow ? insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        const Scalar velocity = globalFaceVars.faceVars(scvf.dofIndex()).velocity();

        const bool insideIsUpstream = sign(scvf.outerNormalScalar()) == sign(velocity);
        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;

        const Scalar upWindWeight = GET_PROP_VALUE(TypeTag, ImplicitUpwindWeight);
        const Scalar upstreamDensity = useMoles ? upstreamVolVars.molarDensity() : upstreamVolVars.density();
        const Scalar downstreamDensity = useMoles ? downstreamVolVars.molarDensity() : downstreamVolVars.density();
        const Scalar upstreamEnthalpy = upstreamVolVars.enthalpy();
        const Scalar downstreamEnthalpy = downstreamVolVars.enthalpy();

//        flux[massBalanceIdx] = TODO??
        flux[energyBalanceIdx] = (upWindWeight * upstreamDensity * upstreamEnthalpy
                                 + (1.0 - upWindWeight) * downstreamDensity * downstreamEnthalpy)
                                 * velocity;

        flux *= scvf.area() * sign(scvf.outerNormalScalar());
        return flux;
    }

    CellCenterPrimaryVariables diffusiveFluxForCellCenter_(const FluxVariables& fluxVars)
    {
        CellCenterPrimaryVariables flux(0.0);

        // compute diffusive flux --> no diffusive flux (only 1 component)

        // compute conductive flux
        computeConductiveFlux_(flux, fluxVars);

        return flux;
    }

    void computeConductiveFlux_(CellCenterPrimaryVariables& flux, FluxVariables& fluxVars)
    {
        flux[energyBalanceIdx] -= fluxVars.temperatureGrad() * fluxVars.face().normal
                   * (fluxVars.thermalConductivity() + fluxVars.thermalEddyConductivity());
    }
};

} // end namespace

#endif
