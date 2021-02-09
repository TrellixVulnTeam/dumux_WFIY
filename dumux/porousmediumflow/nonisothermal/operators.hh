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
 * \ingroup NIModel
 * \brief Sub-control entity-local evaluation of the operators
 *        involved in the system of equations non-isothermal models
 *        that assume thermal equilibrium between all phases.
 */
#ifndef DUMUX_FV_NON_ISOTHERMAL_OPERATORS_HH
#define DUMUX_FV_NON_ISOTHERMAL_OPERATORS_HH

namespace Dumux {

/*!
 * \ingroup NIModel
 * \brief Sub-control entity-local evaluation of the operators
 *        involved in the system of equations of non-isothermal models
 *        that assume thermal equilibrium between all phases.
 * \tparam FluxVariables the type that is responsible for computing the individual
 *                       flux contributions, i.e., advective, diffusive, conduvtive...
 * \tparam ElementVariables the type of element-local view on the grid variables
 */
template<class ModelTraits, class ElementVariables>
class FVNonIsothermalOperators
{
    // The variables required for the evaluation of the equation
    using GridVariables = typename ElementVariables::GridVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename ElementVariables::ElementVolumeVariables;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int energyEqIdx = ModelTraits::Indices::energyEqIdx;

public:

    /*!
     * \brief The energy storage in the fluid phase with index phaseIdx
     *
     * \param storage The mass of the component within the sub-control volume
     * \param scv The sub-control volume
     * \param volVars The volume variables
     * \param phaseIdx The phase index
     */
    template<class NumEqVector>
    static void addPhaseStorage(NumEqVector& storage,
                                const SubControlVolume& scv,
                                const VolumeVariables& volVars,
                                int phaseIdx)
    {
        storage[energyEqIdx] += volVars.porosity()
                                * volVars.density(phaseIdx)
                                * volVars.internalEnergy(phaseIdx)
                                * volVars.saturation(phaseIdx);
    }

    /*!
     * \brief The energy storage in the solid matrix
     *
     * \param storage The mass of the component within the sub-control volume
     * \param scv The sub-control volume
     * \param volVars The volume variables
     */
    template<class NumEqVector>
    static void addSolidPhaseStorage(NumEqVector& storage,
                                     const SubControlVolume& scv,
                                     const VolumeVariables& volVars)
    {
        storage[energyEqIdx] += volVars.temperature()
                                * volVars.solidHeatCapacity()
                                * volVars.solidDensity()
                                * (1.0 - volVars.porosity());
    }

    /*!
     * \brief The advective phase energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     * \param phaseIdx The phase index
     */
    template<class NumEqVector, class FluxVariables>
    static void addHeatConvectionFlux(NumEqVector& flux,
                                      FluxVariables& fluxVars,
                                      int phaseIdx)
    {
        auto upwindTerm = [phaseIdx](const auto& volVars)
        { return volVars.density(phaseIdx)*volVars.mobility(phaseIdx)*volVars.enthalpy(phaseIdx); };

        flux[energyEqIdx] += fluxVars.advectiveFlux(phaseIdx, upwindTerm);
    }

    /*!
     * \brief The diffusive energy fluxes
     *
     * \param flux The flux
     * \param fluxVars The flux variables.
     */
    template<class NumEqVector, class FluxVariables>
    static void addHeatConductionFlux(NumEqVector& flux,
                                      FluxVariables& fluxVars)
    {
        flux[energyEqIdx] += fluxVars.heatConductionFlux();
    }

    /*!
     * \brief heat transfer between the phases for nonequilibrium models
     *
     * \param source The source which ought to be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scv The sub-control volume over which we integrate the source term
     */
    template<class NumEqVector>
    static void addSourceEnergy(NumEqVector& source,
                                const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolume &scv)
    {}
};

} // end namespace Dumux

#endif
