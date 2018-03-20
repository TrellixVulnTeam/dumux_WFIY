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
 * \ingroup NavierStokesNCModel
 *
 * \copydoc Dumux::NavierStokesNCVolumeVariables
 */
#ifndef DUMUX_NAVIER_STOKES_NC_VOLUMEVARIABLES_HH
#define DUMUX_NAVIER_STOKES_NC_VOLUMEVARIABLES_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/volumevariables.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCModel
 * \brief Volume variables for the single-phase, multi-component Navier-Stokes model.
 */
template <class TypeTag>
class NavierStokesNCVolumeVariables : public NavierStokesVolumeVariables<TypeTag>
{
    using ParentType = NavierStokesVolumeVariables<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents(),
           numPhases = FluidSystem::numPhases,
           mainCompIdx = Indices::mainCompIdx,
           pressureIdx = Indices::pressureIdx
    };

    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr auto phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);

public:

    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    template<class ElementSolution>
    void update(const ElementSolution &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume& scv)
    {
        this->extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);

        completeFluidState(elemSol, problem, element, scv, this->fluidState_);


        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(this->fluidState_);
        int compIIdx = phaseIdx;
        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            // binary diffusion coefficents
            if(compIIdx!= compJIdx)
            {
                setDiffusionCoefficient_(phaseIdx, compJIdx,
                                         FluidSystem::binaryDiffusionCoefficient(this->fluidState_,
                                                                                 paramCache,
                                                                                 phaseIdx,
                                                                                 compIIdx,
                                                                                 compJIdx));
            }
        }
    };

    /*!
     * \brief Update the fluid state
     */
    template<class ElementSolution>
    static void completeFluidState(const ElementSolution& elemSol,
                                   const Problem& problem,
                                   const Element& element,
                                   const SubControlVolume& scv,
                                   FluidState& fluidState)
    {
        fluidState.setTemperature(ParentType::temperature(elemSol, problem, element, scv));
        fluidState.setPressure(phaseIdx, elemSol[0][Indices::pressureIdx]);

        // saturation in a single phase is always 1 and thus redundant
        // to set. But since we use the fluid state shared by the
        // immiscible multi-phase models, so we have to set it here...
        fluidState.setSaturation(phaseIdx, 1.0);

        Scalar fracMinor = 0.0;
        int transportEqIdx = 1;

        for(int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if(compIdx == mainCompIdx)
                continue;

            const Scalar moleOrMassFraction = elemSol[0][transportEqIdx++] + 1.0;
            if(useMoles)
                fluidState.setMoleFraction(phaseIdx, compIdx, moleOrMassFraction -1.0);
            else
                fluidState.setMassFraction(phaseIdx, compIdx, moleOrMassFraction -1.0);
            fracMinor += moleOrMassFraction - 1.0;
        }
        if(useMoles)
            fluidState.setMoleFraction(phaseIdx, mainCompIdx, 1.0 - fracMinor);
        else
            fluidState.setMassFraction(phaseIdx, mainCompIdx, 1.0 - fracMinor);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
        fluidState.setDensity(phaseIdx, value);

        value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        fluidState.setViscosity(phaseIdx, value);

        // compute and set the enthalpy
        const Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
        fluidState.setEnthalpy(phaseIdx, h);
    }


     /*!
      * \brief Returns the mass fraction of a component in the phase \f$\mathrm{[-]}\f$
      *
      * \param pIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int pIdx, int compIdx) const
     {
         assert(pIdx == phaseIdx);
         return this->fluidState_.massFraction(pIdx, compIdx);
     }

     /*!
      * \brief Returns the mole fraction of a component in the phase \f$\mathrm{[-]}\f$
      *
      * \param pIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar moleFraction(int pIdx, int compIdx) const
     {
         assert(pIdx == phaseIdx);
         return this->fluidState_.moleFraction(pIdx, compIdx);
     }

    /*!
     * \brief Returns the mass density of a given phase \f$\mathrm{[kg/m^3]}\f$
     *
     * \param pIdx the index of the fluid phase
     */
    Scalar molarDensity(int pIdx = phaseIdx) const
    {
        assert(pIdx == phaseIdx);
        return this->fluidState_.molarDensity(pIdx);
    }

     /*!
     * \brief Returns the diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     */
    Scalar diffusionCoefficient(int pIdx, int compIdx) const
    {
        assert(pIdx == phaseIdx);
        if (compIdx < pIdx)
            return diffCoefficient_[pIdx][compIdx];
        else if (compIdx > pIdx)
            return diffCoefficient_[pIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient called for phaseIdx = compIdx");
    }

     /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$
     */
    Scalar effectiveDiffusivity(int pIdx, int compIdx) const
    {
        return diffusionCoefficient(pIdx, compIdx);
    }

protected:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void setDiffusionCoefficient_(int pIdx, int compIdx, Scalar d)
    {
        assert(pIdx == phaseIdx);
        if (compIdx < pIdx)
            diffCoefficient_[pIdx][compIdx] = std::move(d);
        else if (compIdx > pIdx)
            diffCoefficient_[pIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coefficient for phaseIdx = compIdx doesn't exist");
    }

    std::array<std::array<Scalar, numComponents-1>, numPhases> diffCoefficient_;
};

} // end namespace Dumux

#endif
