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
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system without using
 *        any assumptions.
 */
#ifndef DUMUX_NONEQUILIBRIUM_MASS_FLUID_STATE_HH
#define DUMUX_NONEQUILIBRIUM_MASS_FLUID_STATE_HH

#include <dumux/common/valgrind.hh>
#include <dumux/material/fluidstates/nonequilibriumfluidstate.hh>
#include <dumux/common/propertysystem.hh>

#include <cmath>
#include <algorithm>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
}

/*!
 * \brief Represents all relevant thermodynamic quantities of a
 *        multi-phase, multi-component fluid system
 *
 *        This fluidstate ought to be used for the case of:
 *        - local thermal equilibrium
 *        - local chemical non-equilibrium
 */
template <class TypeTag>
class NonEquilibriumMassFluidState :
    public NonEquilibriumFluidState<typename GET_PROP_TYPE(TypeTag, Scalar),
                                    typename GET_PROP_TYPE(TypeTag, FluidSystem)>
    {
    public:
        typedef NonEquilibriumFluidState<typename GET_PROP_TYPE(TypeTag, Scalar),
                                         typename GET_PROP_TYPE(TypeTag, FluidSystem)> ParentType;

        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        enum {numEnergyEqs    = Indices::numPrimaryEnergyVars};
        enum { numPhases       = GET_PROP_VALUE(TypeTag, NumPhases)};
        enum { numComponents   = GET_PROP_VALUE(TypeTag, NumComponents)};

        NonEquilibriumMassFluidState()
        : NonEquilibriumFluidState<Scalar, FluidSystem>()
    {
    }

    /*****************************************************
     * Setter methods. Note that these are not part of the
     * generic FluidState interface but specific for each
     * implementation...
     *****************************************************/
    /*!
     * \brief Set the temperature [K] of a fluid phase
     *        Both versions of the function need to be here.
     *        Otherwise the compiler gets confused.
     *        Thus, this is just forwarding to the Parent
     *        (unclear why this is necessary).
     */
    void setTemperature(const int phaseIdx, const Scalar value)
    {
        // Unfortunately throw does not work when triggered from a constructor
        std::cout <<"file: "<< __FILE__ << ", line: " << __LINE__ <<". This is a fluidstate for *chemical* non-equilibrium, not thermal! \n" ;
        DUNE_THROW(Dune::NotImplemented, "This is a fluidstate for *chemical* non-equilibrium, not thermal!");
    }

    /*!
     * \brief Set the temperature [K] of the fluid phases.
     *        Both versions of the function need to be here.
     *        Otherwise the compiler gets confused.
     *        Thus, this is just presenting the signature to the compiler.
     */
    void setTemperature(const Scalar value)
    {
        temperature_ = value;
    }

    /*!
     * \brief Get the temperature [K] of the fluid phases.
     */
    Scalar temperature() const
    {
        return temperature_ ;
    }

    /*!
     * \brief Get the temperature [K] of the fluid phases.
     */
    Scalar temperature(const int dummy) const
    {
        return temperature_ ;
    }

    /*!
      * \brief Retrieve all parameters from an arbitrary fluid
      *        state. The assign method from the parent class cannot be used, because here, we have only one temperature.
      */
     template <class FluidState>
     void assign(const FluidState& fs)
     {
         for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
             ParentType::averageMolarMass_[phaseIdx] = 0;
             ParentType::sumMoleFractions_[phaseIdx] = 0;
             for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                 ParentType::moleFraction_[phaseIdx][compIdx] = fs.moleFraction(phaseIdx, compIdx);
                 ParentType::fugacityCoefficient_[phaseIdx][compIdx] = fs.fugacityCoefficient(phaseIdx, compIdx);
                 ParentType::averageMolarMass_[phaseIdx] += ParentType::moleFraction_[phaseIdx][compIdx]*FluidSystem::molarMass(compIdx);
                 ParentType::sumMoleFractions_[phaseIdx] += ParentType::moleFraction_[phaseIdx][compIdx];
             }

             ParentType::pressure_[phaseIdx] = fs.pressure(phaseIdx);
             ParentType::saturation_[phaseIdx] = fs.saturation(phaseIdx);
             ParentType::density_[phaseIdx] = fs.density(phaseIdx);
             ParentType::enthalpy_[phaseIdx] = fs.enthalpy(phaseIdx);
             ParentType::viscosity_[phaseIdx] = fs.viscosity(phaseIdx);
         }
         temperature_ = fs.temperature(/*phaseIdx=*/0); // in this fluidstate there is only one temperature.

     }



private:
    Scalar temperature_;
};

} // end namepace Dumux

#endif
