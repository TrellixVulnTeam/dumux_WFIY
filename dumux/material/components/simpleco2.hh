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
 * \ingroup Components
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure CO2.
 */
#ifndef DUMUX_SIMPLE_CO2_HH
#define DUMUX_SIMPLE_CO2_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/idealgas.hh>

#include <cmath>

#include <dumux/material/components/base.hh>
#include <dumux/material/components/gas.hh>

namespace Dumux::Components {

/*!
 * \ingroup Components
 * \brief A simple version of pure CO2
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class SimpleCO2
: public Components::Base<Scalar, SimpleCO2<Scalar> >
, public Components::Gas<Scalar, SimpleCO2<Scalar> >
{
    using IdealGas = Dumux::IdealGas<Scalar>;

    static const Scalar R;  // specific gas constant of CO2

    static constexpr Scalar limitPressure = 10.0e5; // [Pa]

public:
    /*!
     * \brief A human readable name for the CO2.
     */
    static std::string name()
    { return "SimpleCO2"; }

    /*!
     * \brief The mass in \f$\mathrm{[kg/mol]}\f$ of one mole of CO2.
     */
    static constexpr Scalar molarMass()
    { return 44e-3; /* [kg/mol] */ }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of CO2
     */
    static Scalar criticalTemperature()
    { return 273.15 + 30.95; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of CO2
     */
    static Scalar criticalPressure()
    { return 73.8e5; /* [Pa] */ }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at CO2's triple point.
     */
    static Scalar tripleTemperature()
    { return 273.15 - 56.35; /* [K] */ }

    /*!
     * \brief Returns the pressure \f$\mathrm{[Pa]}\f$ at CO2's triple point.
     */
    static Scalar triplePressure()
    { return 5.11e5; /* [N/m^2] */ }


    /*!
     * \brief Specific enthalpy of CO2 \f$\mathrm{[J/kg]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasEnthalpy(Scalar temperature,
                                    Scalar pressure)
    {
        if (pressure > limitPressure)
            DUNE_THROW(Dune::RangeError, "Ideal gas law not applicable at pressure " << pressure);
        static const Scalar tRef = getParam<Scalar>("SimpleCO2.ReferenceTemperature", 293.15);
        return gasHeatCapacity(temperature, pressure)*(temperature - tRef); // + vaporizationEnthalpy();
    }

    /*!
     * \brief Specific internal energy of CO2 \f$\mathrm{[J/kg]}\f$.
     *
     *        Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     *        Rearranging for internal energy yields: \f$u = h - pv\f$.
     *        Exploiting the Ideal Gas assumption (\f$pv = R_{\textnormal{specific}} T\f$)gives: \f$u = h - R / M T \f$.
     *
     *        The universal gas constant can only be used in the case of molar formulations.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        if (pressure > limitPressure)
            DUNE_THROW(Dune::RangeError, "Ideal gas law not applicable at pressure " << pressure);
        // 1/molarMass: conversion from [J/(mol K)] to [J/(kg K)]
        // R*T/molarMass: pressure *spec. volume for an ideal gas
        return gasEnthalpy(temperature, pressure)
                - 1/molarMass()*IdealGas::R*temperature;
    }

    /*!
     * \brief Returns true if the gas phase is assumed to be compressible
     */
    static constexpr bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true if the gas phase viscostiy is constant
     */
    static constexpr bool gasViscosityIsConstant()
    { return false; }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of CO2 at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        if (pressure > limitPressure)
            DUNE_THROW(Dune::RangeError, "Ideal gas law not applicable at pressure " << pressure);
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     *  \brief The molar density of CO2 in \f$\mathrm{[mol/m^3]}\f$ at a given pressure and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     */
    static Scalar gasMolarDensity(Scalar temperature, Scalar pressure)
    { return IdealGas::molarDensity(temperature, pressure); }

    /*!
     * \brief Returns true if the gas phase is assumed to be ideal
     */
    static constexpr bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure of CO2 in \f$\mathrm{[Pa]}\f$ at a given density and temperature.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        const Scalar pressure = IdealGas::pressure(temperature, density/molarMass());
        if (pressure > limitPressure)
            DUNE_THROW(Dune::RangeError, "Ideal gas law not applicable at pressure " << pressure);
        return pressure;
    }

    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of CO2.
     * Equations given in: - Vesovic et al., 1990
     *                     - Fenhour et al., 1998
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     * TODO: this does not look like a really "simple" parameterization. Can this be simplified further?
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {
        static const double a0 = 0.235156;
        static const double a1 = -0.491266;
        static const double a2 = 5.211155E-2;
        static const double a3 = 5.347906E-2;
        static const double a4 = -1.537102E-2;

        static const double d11 = 0.4071119E-2;
        static const double d21 = 0.7198037E-4;
        static const double d64 = 0.2411697E-16;
        static const double d81 = 0.2971072E-22;
        static const double d82 = -0.1627888E-22;

        static const double ESP = 251.196;

        double mu0, SigmaStar, TStar;
        double dmu, rho;
        double visco_CO2;

        if(temperature < 275.) // regularisation
        {
            temperature = 275;
            Dune::dgrave << "Temperature below 275K in viscosity function:"
                    << "Regularizing tempereature to 275K. " << std::endl;
        }


        TStar = temperature/ESP;

        /* mu0: viscosity in zero-density limit */
        using std::exp;
        using std::log;
        using std::sqrt;
        SigmaStar = exp(a0 + a1*log(TStar)
                        + a2*log(TStar)*log(TStar)
                        + a3*log(TStar)*log(TStar)*log(TStar)
                        + a4*log(TStar)*log(TStar)*log(TStar)*log(TStar) );
        mu0 = 1.00697*sqrt(temperature) / SigmaStar;

        /* dmu : excess viscosity at elevated density */
        rho = gasDensity(temperature, pressure); /* CO2 mass density [kg/m^3] */

        using Dune::power;
        dmu = d11*rho + d21*rho*rho + d64*power(rho,6)/(TStar*TStar*TStar)
            + d81*power(rho,8) + d82*power(rho,8)/TStar;

        visco_CO2 = (mu0 + dmu)/1.0E6;   /* conversion to [Pa s] */

        return visco_CO2;
    }


    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of CO2.
     *
     * Thermal conductivity of CO2 at T=20°C, see:
     * http://www.engineeringtoolbox.com/carbon-dioxide-d_1000.html
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        return 0.087;
    }

    /*!
     * \brief Specific isobaric heat capacity of CO2 \f$\mathrm{[J/(kg*K)]}\f$.
     *        source: //TODO
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasHeatCapacity(Scalar temperature, Scalar pressure)
    {
        return 850;  //TODO!! from old simpleco2.hh file.
    }

};

template <class Scalar>
const Scalar Components::SimpleCO2<Scalar>::R = Constants<Scalar>::R / Components::SimpleCO2<Scalar>::molarMass();

} // end namespace Dumux::Components

#endif
