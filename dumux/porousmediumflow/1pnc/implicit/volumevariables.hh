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
 *
 * \brief Quantities required by the single-phase, n-component box
 *        model defined on a vertex.
 */
#ifndef DUMUX_1PNC_VOLUME_VARIABLES_HH
#define DUMUX_1PNC_VOLUME_VARIABLES_HH

#include <dumux/discretization/volumevariables.hh>

#include "properties.hh"



namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(IsothermalVolumeVariables);
}

/*!
 * \ingroup OnePNCModel
 * \ingroup ImplicitVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the one-phase, n-component model.
 *
 * \note The return functions for the fluid state variables always forward to the actual
 *       fluid state using the phaseIdx from the DuMuX property system. Furthermore, the
 *       default value is not used, but is only here to enable calling these functions
 *       without handing in a phase index (as in a single-phasic context there is only one phase).
 *       This way one can use two-phase fluid systems for this one-phasic flow and transport
 *       model by specifying which phase is present through the DuMuX property system.
 */
template <class TypeTag>
class OnePNCVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    using ParentType = ImplicitVolumeVariables<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PermeabilityType = typename SpatialParams::PermeabilityType;

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        phaseIdx = Indices::phaseIdx,
        phaseCompIdx = Indices::phaseCompIdx,
        firstTransportEqIdx = Indices::firstTransportEqIdx,

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        firstMoleFracIdx = Indices::firstMoleFracIdx,

    };

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;

// new
    using DimVector = Dune::FieldVector<Scalar,dim>;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
//     typedef typename Grid::ctype CoordScalar; ?? deprecated ?

//     enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
//     enum { dofCodim = isBox ? dim : 0 };
public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \copydoc ImplicitVolumeVariables::update
     * \param priVars The primary Variables
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv
                /*bool isOldSol ??*/)
    {
        ParentType::update(elemSol, problem, element, scv);

        completeFluidState(elemSol, problem, element, scv, fluidState_/*, isOldSol*/);

        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
        // dispersivity_ = problem.spatialParams().dispersivity(element, scv, elemSol);
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);

        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState_, phaseIdx);

        int compIIdx = phaseCompIdx;

        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            diffCoeff_[compJIdx] = 0.0;
            if(compIIdx!= compJIdx)
                {
                diffCoeff_[compJIdx] = FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             compIIdx,
                                                             compJIdx);
                }
//             Valgrind::CheckDefined(diffCoeff_[compJIdx]);
        }

        // energy related quantities not contained in the fluid state
//         asImp_().callProtectedUpdateEnergy(elemSol, problem, element, fvGeometry, scvIdx, isOldSol);
    }

   /*!
    * \copydoc ImplicitModel::completeFluidState
    * \param isOldSol Specifies whether this is the previous solution or the current one
    * \param priVars The primary Variables
    */
    static void completeFluidState(const ElementSolutionVector &elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume &scv,
                                   FluidState& fluidState)

    {
//    old
//         Scalar t = IsothermalVolumeVariables::callProtectedTemperature(priVars, problem, element,
//                                                                        fvGeometry, scvIdx);
        Scalar t = ParentType::temperature(elemSol, problem, element, scv);
        fluidState.setTemperature(t);
        fluidState.setSaturation(phaseIdx, 1.);

        const auto& priVars = ParentType::extractDofPriVars(elemSol, scv); // new in next ??
        fluidState.setPressure(phaseIdx, priVars[pressureIdx]);

        // calculate the phase composition

        Dune::FieldVector<Scalar, numComponents> moleFrac;

        Scalar sumMoleFracNotWater = 0;

        for (int compIdx=firstMoleFracIdx; compIdx<numComponents; ++compIdx){
             moleFrac[compIdx] = priVars[compIdx];
             sumMoleFracNotWater +=moleFrac[compIdx];
            }
        moleFrac[0] = 1- sumMoleFracNotWater;

        // Set fluid state mole fractions
        for (int compIdx=0; compIdx<numComponents; ++compIdx)
        {
            fluidState.setMoleFraction(phaseIdx, compIdx, moleFrac[compIdx]);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
        Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);

        fluidState.setDensity(phaseIdx, rho);
        fluidState.setViscosity(phaseIdx, mu);

        // compute and set the enthalpy
        Scalar h = Implementation::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);
    }

    /*!
     * \brief Return the fluid configuration at the given primary
     *        variables
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar density(int phaseIdx) const
    { return fluidState_.density(phaseIdx); }

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ the of the fluid phase.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar molarDensity(int phaseIdx = 0) const
    { return fluidState_.molarDensity(phaseIdx); }

    /*!
     * \brief Return the saturation
     *
     * This method is here for compatibility reasons with other models. The saturation
     * is always 1.0 in a one-phasic context.
     */
    Scalar saturation(int pIdx = 0) const
    { return 1.0; }

     /*!
      * \brief Return mole fraction \f$\mathrm{[mol/mol]}\f$ of a component in the phase.
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      *
      * We always forward to the fluid state with the phaseIdx property (see class description).
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     { return fluidState_.moleFraction(phaseIdx, compIdx); }

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     { return fluidState_.massFraction(phaseIdx, compIdx); }

    /*!
     * \brief Return the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar pressure(int pIdx = 0) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Return temperature \f$\mathrm{[K]}\f$ inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(); }

    /*!
     * \brief Returns the mobility \f$\mathrm{[1/(Pa s)]}\f$.
     *
     * The term mobility is usually not employed in the one phase context.
     * The method is here for compatibility reasons with other models.
     *
     * We always forward to the fluid state with the phaseIdx property (see class description).
     */
    Scalar mobility(int pIdx = 0) const
    { return 1.0/fluidState_.viscosity(phaseIdx); }

    /*!
     * \brief Return the average porosity \f$\mathrm{[-]}\f$ within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }


//     /*!
//      * \brief Returns the binary diffusion coefficients for a phase in \f$[m^2/s]\f$.
//      */
//     Scalar diffCoeff(int compIdx) const
//     { return diffCoeff_[compIdx]; }

    /*!
     * \brief Return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ in the fluid.
     */
    Scalar diffusionCoefficient(int pIdx, int compIdx) const
    { return diffCoeff_[compIdx]; }

    /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     */
     Scalar molarity(int compIdx) const // [moles/m^3]
    { return fluidState_.molarity(phaseIdx, compIdx);}
// old    { return this->fluidState_.molarity(phaseIdx, compIdx);}
     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
//      Scalar massFraction(int compIdx) const
//      {
//         return this->fluidState_.massFraction(phaseIdx, compIdx);
//      }


   /*!
    * Circumvents the inheritance architecture of the ninisothermal model
    */
//     static Scalar callProtectedTemperature(const PrimaryVariables &priVars,
//                                            const Problem& problem,
//                                            const Element &element,
//                                            const FVElementGeometry &fvGeometry,
//                                            int scvIdx)
//     {
//          return Implementation::temperature_(priVars, problem,element, fvGeometry, scvIdx);
//     }
//
//    /*!
//     * Circumvents the inheritance architecture of the ninisothermal model
//     */
//    void callProtectedUpdateEnergy(const PrimaryVariables &priVars,
//                                    const Problem &problem,
//                                    const Element &element,
//                                    const FVElementGeometry &fvGeometry,
//                                    const int scvIdx,
//                                    bool isOldSol)
//     {
//         asImp_().updateEnergy_(priVars, problem,element, fvGeometry, scvIdx, isOldSol);
//     };

    /*!
     * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
     */
    const PermeabilityType& permeability() const
    { return permeability_; }

protected:

    static Scalar temperature_(const PrimaryVariables &priVars,
                               const Problem& problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               int scvIdx)
    {
         return problem.temperatureAtPos(fvGeometry.subContVol[scvIdx].global);
    }

//     template<class ParameterCache>
//     static Scalar enthalpy_(const FluidState& fluidState,
//                             const ParameterCache& paramCache,
//                             int phaseIdx)
//     {
//         return 0;
//     }

    /*!
        * \brief Called by update() to compute the energy related quantities
        */
//     void updateEnergy_(const PrimaryVariables &priVars,
//                         const Problem &problem,
//                         const Element &element,
//                         const FVElementGeometry &fvGeometry,
//                         const int scvIdx,
//                         bool isOldSol)
//     { };
//
    Scalar porosity_;        //!< Effective porosity within the control volume
    PermeabilityType permeability_;
    Scalar density_;
    FluidState fluidState_;
    Dune::FieldVector<Scalar, numComponents> diffCoeff_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};

} // end namespace

#endif
