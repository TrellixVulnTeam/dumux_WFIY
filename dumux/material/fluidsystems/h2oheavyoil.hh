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
 * \brief A fluid system with water, gas and NAPL as phases and
 *        \f$H_2O\f$ and \f$NAPL (contaminant)\f$ as components.
 * It uses heavyoil Properties, but allows for a scaling of density
 * thus enabling an artifical DNAPL if desired
 */
#ifndef DUMUX_H2O_HEAVYOIL_FLUID_SYSTEM_HH
#define DUMUX_H2O_HEAVYOIL_FLUID_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/heavyoil.hh>

#include <dumux/material/binarycoefficients/h2o_heavyoil.hh>

#include <dumux/material/fluidsystems/base.hh>

namespace Dumux
{
namespace FluidSystems
{

/*!
 * \brief A compositional fluid with water and heavy oil
 *        components in both, the liquid and the gas phase.
 */
template <class TypeTag, class Scalar,
          class H2OType = Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> > >
class H2OHeavyOil
    : public BaseFluidSystem<Scalar, H2OHeavyOil<TypeTag, Scalar, H2OType> >
{
    typedef H2OHeavyOil<TypeTag, Scalar, H2OType> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    typedef Dumux::HeavyOil<Scalar> HeavyOil;
    typedef H2OType H2O;


    static const int numPhases = 3;
    static const int numComponents = 2;

    static const int wPhaseIdx = 0; // index of the water phase
    static const int nPhaseIdx = 1; // index of the NAPL phase
    static const int gPhaseIdx = 2; // index of the gas phase

    static const int H2OIdx = 0;
    static const int NAPLIdx = 1;

    // export component indices to indicate the main component
    // of the corresponding phase at atmospheric pressure 1 bar
    // and room temperature 20°C:
    static const int wCompIdx = H2OIdx;
    static const int nCompIdx = NAPLIdx;

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            std::cout << "Initializing tables for the H2O fluid properties ("
                      << nTemp*nPress
                      << " entries).\n";

            H2O::init(tempMin, tempMax, nTemp,
                      pressMin, pressMax, nPress);
        }
    }


    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gPhaseIdx;
    }

    static bool isIdealGas(int phaseIdx)
    { return phaseIdx == gPhaseIdx && H2O::gasIsIdeal() && HeavyOil::gasIsIdeal(); }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Raoult's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Raoult's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        if (phaseIdx == gPhaseIdx)
            return true;
        else if (phaseIdx == wPhaseIdx)
            // the water component decides for the water phase...
            return H2O::liquidIsCompressible();

        // the NAPL component decides for the napl phase...
        return HeavyOil::liquidIsCompressible();
    }

    /*!
     * \brief Return the human readable name of a phase (used in indices)
     */
    static std::string phaseName(int phaseIdx)
    {
        switch (phaseIdx) {
        case wPhaseIdx: return "w";
        case nPhaseIdx: return "n";
        case gPhaseIdx: return "g";;
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Return the human readable name of a component (used in indices)
     */
    static std::string componentName(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case NAPLIdx: return HeavyOil::name();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
        case NAPLIdx: return HeavyOil::molarMass();
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    /*!
     * \brief Given all mole fractions in a phase, return the phase
     *        density [kg/m^3].
     */
    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState, int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            // See: doctoral thesis of Steffen Ochs 2007
            // Steam injection into saturated porous media : process analysis including experimental and numerical investigations
            // http://elib.uni-stuttgart.de/bitstream/11682/271/1/Diss_Ochs_OPUS.pdf
            Scalar rholH2O = H2O::liquidDensity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
            Scalar clH2O = rholH2O/H2O::molarMass();

            // this assumes each dissolved molecule displaces exactly one
            // water molecule in the liquid
            return
                clH2O*(H2O::molarMass()*fluidState.moleFraction(wPhaseIdx, H2OIdx)
                       +
                       HeavyOil::molarMass()*fluidState.moleFraction(wPhaseIdx, NAPLIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            Scalar pressure = HeavyOil::liquidIsCompressible()?fluidState.pressure(phaseIdx):1e100;
            return HeavyOil::liquidDensity(fluidState.temperature(phaseIdx), pressure);
        }

        assert (phaseIdx == gPhaseIdx);
        Scalar pH2O =
            fluidState.moleFraction(gPhaseIdx, H2OIdx)  *
            fluidState.pressure(gPhaseIdx);
        Scalar pNAPL =
            fluidState.moleFraction(gPhaseIdx, NAPLIdx)  *
            fluidState.pressure(gPhaseIdx);
        return
            H2O::gasDensity(fluidState.temperature(phaseIdx), pH2O) +
            HeavyOil::gasDensity(fluidState.temperature(phaseIdx), pNAPL);
    }

    /*!
     * \brief Return the viscosity of a phase.
     */
    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            // assume pure water viscosity
            return H2O::liquidViscosity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            // assume pure NAPL viscosity
            return HeavyOil::liquidViscosity(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }

        assert (phaseIdx == gPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         *
         * in this case, we use a simplified version in order to avoid
         * computationally costly evaluation of sqrt and pow functions and
         * divisions
         * -- compare e.g. with Promo Class p. 32/33
         */
        const Scalar mu[numComponents] = {
            H2O::gasViscosity(fluidState.temperature(phaseIdx), H2O::vaporPressure(fluidState.temperature(phaseIdx))),
            HeavyOil::gasViscosity(fluidState.temperature(phaseIdx), HeavyOil::vaporPressure(fluidState.temperature(phaseIdx)))
        };

        return mu[H2OIdx]*fluidState.moleFraction(gPhaseIdx, H2OIdx)
                + mu[NAPLIdx]*fluidState.moleFraction(gPhaseIdx, NAPLIdx);
    }


    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
            return BinaryCoeff::H2O_HeavyOil::liquidDiffCoeff(T, p);

        // gas phase
        else if (phaseIdx == gPhaseIdx)
            return BinaryCoeff::H2O_HeavyOil::gasDiffCoeff(T, p);

        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent diffusion coefficient for phase index " << phaseIdx);
    }

     /* Henry coefficients
     */
    template <class FluidState>
    static Scalar henryCoefficient(const FluidState &fluidState,
                                   int phaseIdx,
                                   int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (compIdx == NAPLIdx && phaseIdx == wPhaseIdx)
            return Dumux::BinaryCoeff::H2O_HeavyOil::henryOilInWater(T)/p;

        else if (phaseIdx == nPhaseIdx && compIdx == H2OIdx)
            return Dumux::BinaryCoeff::H2O_HeavyOil::henryWaterInOil(T)/p;

        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent henry coefficient for phase index " << phaseIdx
                                                     << " and component index " << compIdx);
    }

     /*  partial pressures in the gas phase, taken from saturation vapor pressures
     */
    template <class FluidState>
    static Scalar partialPressureGas(const FluidState &fluidState, int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar T = fluidState.temperature(phaseIdx);
        if (compIdx == NAPLIdx)
            return HeavyOil::vaporPressure(T);
        else if (compIdx == H2OIdx)
            return H2O::vaporPressure(T);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }

     /*  inverse vapor pressures, taken from inverse saturation vapor pressures
     */
    template <class FluidState>
    static Scalar inverseVaporPressureCurve(const FluidState &fluidState,
                                            int phaseIdx,
                                            int compIdx)
    {
        assert(0 <= compIdx  && compIdx < numComponents);

        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (compIdx == NAPLIdx)
            return HeavyOil::vaporTemperature(pressure);
        else if (compIdx == H2OIdx)
            return H2O::vaporTemperature(pressure);
        else
            DUNE_THROW(Dune::InvalidStateException, "non-existent component index " << compIdx);
    }



    /*!
     * \brief Given all mole fractions in a phase, return the specific
     *        phase enthalpy [J/kg].
     */
    /*!
     *  \todo This system neglects the contribution of gas-molecules in the liquid phase.
     *        This contribution is probably not big. Somebody would have to find out the enthalpy of solution for this system. ...
     */
    using Base::enthalpy;
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == nPhaseIdx) {
            return HeavyOil::liquidEnthalpy(fluidState.temperature(phaseIdx), fluidState.pressure(phaseIdx));
        }
        else if (phaseIdx == gPhaseIdx) {  // gas phase enthalpy depends strongly on composition
            Scalar hgc = HeavyOil::gasEnthalpy(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
            Scalar hgw = H2O::gasEnthalpy(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx));

            Scalar result = 0;
            result += hgw * fluidState.massFraction(gPhaseIdx, H2OIdx);
            result += hgc * fluidState.massFraction(gPhaseIdx, NAPLIdx);

            return result;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::heatCapacity;
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "FluidSystems::H2ONAPL::heatCapacity()");
    }

    using Base::thermalConductivity;
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      int phaseIdx)
    {
        const Scalar temperature  = fluidState.temperature(phaseIdx) ;
        const Scalar pressure = fluidState.pressure(phaseIdx);
        if (phaseIdx == wPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == nPhaseIdx)
        {
            return HeavyOil::liquidThermalConductivity(temperature, pressure);
        }
        else if (phaseIdx == gPhaseIdx)
        {
            return H2O::gasThermalConductivity(temperature, pressure);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

private:

};
} // end namespace FluidSystems

#ifdef DUMUX_PROPERTIES_HH
// forward definitions of the property tags
namespace Properties {
    NEW_PROP_TAG(Scalar);
    NEW_PROP_TAG(Components);
}

/*!
 * \brief A threephase fluid system with water, heavyoil as components.
 *
 * This is an adapter to use Dumux::H2OheavyoilFluidSystem<TypeTag>, as is
 * done with most other classes in Dumux.
 *  This fluidsystem is applied by default with the tabulated version of
 *  water of the IAPWS-formulation.
 *
 *  To change the component formulation (ie to use nontabulated or
 *  incompressible water), or to switch on verbosity of tabulation,
 *  use the property system and the property "Components":
 *
 *        // Select desired version of the component
 *        SET_PROP(myApplicationProperty, Components) : public GET_PROP(TypeTag, DefaultComponents)
 *        {
 *            typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 *
 *        // Do not use the defaults !
 *        //    typedef Dumux::TabulatedComponent<Scalar, Dumux::H2O<Scalar> > H2O;
 *
 *        // Apply e.g. untabulated water:
 *        typedef Dumux::H2O<Scalar> H2O;
 *        };
 *
 *   Also remember to initialize tabulated components (FluidSystem::init()), while this
 *   is not necessary for non-tabularized ones.
 */
template<class TypeTag>
class H2OHeavyOilFluidSystem
: public FluidSystems::H2OHeavyOil<TypeTag, typename GET_PROP_TYPE(TypeTag, Scalar),
                                        typename GET_PROP(TypeTag, Components)::H2O>
{};
#endif

} // end namespace Dumux

#endif
