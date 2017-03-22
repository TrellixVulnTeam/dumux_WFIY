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
 * \ingroup Components
 *
 * \brief A simple class for the air fluid properties
 */
#ifndef DUMUX_AIR_HH
#define DUMUX_AIR_HH

#include <dumux/common/exceptions.hh>
#include <dumux/material/components/component.hh>
#include <dumux/material/idealgas.hh>

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A class for the air fluid properties
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class Air : public Component<Scalar, Air<Scalar> >
{
    typedef Dumux::IdealGas<Scalar> IdealGas;

public:
    /*!
     * \brief A human readable name for Air.
     */
    static std::string name()
    { return "Air"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Air.
     *
     * Taken from constrelair.hh.
     */
    static Scalar molarMass()
    { return 0.02896; /* [kg/mol] */ }

    /*!
     * \brief Returns the critical temperature \f$\mathrm{[K]}\f$ of Air.
     */
    static Scalar criticalTemperature()
    { return 132.531 ; /* [K] */ }

    /*!
     * \brief Returns the critical pressure \f$\mathrm{[Pa]}\f$ of Air.
     */
    static Scalar criticalPressure()
    { return 37.86e5; /* [Pa] */ }

    /*!
     * \brief The density \f$\mathrm{[kg/m^3]}\f$ of Air at a given pressure and temperature.
     *
     * Ideal gas is assumed.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of phase in \f$\mathrm{[Pa]}\f$
    */
    static Scalar gasDensity(Scalar temperature, Scalar pressure)
    {
        // Assume an ideal gas
        return IdealGas::density(molarMass(), temperature, pressure);
    }

    /*!
     * \brief Returns true, the gas phase is assumed to be compressible
     */
    static bool gasIsCompressible()
    { return true; }

    /*!
     * \brief Returns true, the gas phase is assumed to be ideal
     */
    static bool gasIsIdeal()
    { return true; }

    /*!
     * \brief The pressure \f$\mathrm{[Pa]}\f$ of gaseous Air at a given density and temperature.
     *
     * Ideal gas is assumed.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param density density of component in \f$\mathrm{[kg/m^3]}\f$
     */
    static Scalar gasPressure(Scalar temperature, Scalar density)
    {
        // Assume an ideal gas
        return IdealGas::pressure(temperature, density/molarMass());
    }
    /*!
     * \brief The dynamic viscosity \f$\mathrm{[Pa*s]}\f$ of Air at a given pressure and temperature.
     *
     * Criticial specific volume calculated by \f$V_c = (R*T_c)/p_c\f$.
     *
     * Reid et al. (1987, pp 396-397, 667) \cite reid1987 <BR>
     * Poling et al. (2001, pp 9.7-9.8) \cite poling2001 <BR>
     *
     * Accentric factor taken from: <BR>
     * Adebiyi (2003) \cite adebiyi2003
     *
     * air is a non-polar substance,
     * thus dipole moment mu is zero, as well the dimensionless dipole moment mu_r
     * therefore not considered below
     * the same holds for the correction value kappa for highly polar substances
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasViscosity(Scalar temperature, Scalar pressure)
    {

        const Scalar Tc = criticalTemperature();
        const Scalar Vc = 84.525138; // critical specific volume [cm^3/mol]
        const Scalar omega = 0.078; // accentric factor
        const Scalar M = molarMass() * 1e3; // molar mas [g/mol]

        Scalar Fc = 1 - 0.2756*omega;
        Scalar Tstar = 1.2593 * temperature/Tc;
        Scalar Omega_v =
            1.16145*std::pow(Tstar, -0.14874) +
            0.52487*std::exp(- 0.77320*Tstar) +
            2.16178*std::exp(- 2.43787*Tstar);
        Scalar mu = 40.785 * Fc * std::sqrt(M * temperature)
                    / (std::cbrt(Vc * Vc) * Omega_v);

        // convertion from micro poise to Pa s
        return mu/1e6 / 10;
    }

    // simpler method, from old constrelAir.hh
    static Scalar simpleGasViscosity(Scalar temperature, Scalar pressure)
    {
        if(temperature < 273.15 || temperature > 660.)
        {
            DUNE_THROW(NumericalProblem,
                "simpleGasViscosity: Temperature out of range! (T = " << temperature << " K)");
        }
        return 1.496e-6 * std::sqrt(temperature * temperature * temperature) / (temperature + 120.);

    }

    /*!
     * \brief Specific enthalpy of Air \f$\mathrm{[J/kg]}\f$
     *        with 273.15 \f$ K \f$ as basis. <BR>
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     * Kays et al. (2005, 431ff) \cite kays2005 <BR>
     */
    static Scalar gasEnthalpy(Scalar temperature, Scalar pressure)
    {
        return 1005*(temperature-273.15);
    }

    /*!
     * \brief Specific internal energy of Air \f$\mathrm{[J/kg]}\f$.
     *
     * Definition of enthalpy: \f$h= u + pv = u + p / \rho\f$.
     * Rearranging for internal energy yields: \f$u = h - pv\f$.
     * Exploiting the Ideal Gas assumption
     * (\f$pv = R_{\textnormal{specific}} T\f$) gives: \f$u = h - R / M T \f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     */
    static const Scalar gasInternalEnergy(Scalar temperature,
                                          Scalar pressure)
    {
        return
            gasEnthalpy(temperature, pressure)
            -
            IdealGas::R * temperature // = pressure * molar volume for an ideal gas
            / molarMass(); // conversion from [J/(mol K)] to [J/(kg K)]
    }

    /*!
     * \brief Specific isobaric heat capacity \f$\mathrm{[J/(kg*K)]}\f$ of pure
     *        air.
     *
     *  This methods uses the formula for "zero-pressure" heat capacity that
     *  is only dependent on temperature, because the pressure dependence is rather small.
     *  This one should be accurate for a pressure of 1 atm.
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     * \param pressure pressure of component in \f$\mathrm{[Pa]}\f$
     *
     *  Values taken from Hollis (1996) \cite hollis1996 <BR>
     *  "Tables of Thermal Properties of Gases"
     */
    static const Scalar gasHeatCapacity(Scalar temperature,
                                        Scalar pressure)
    {
        // scale temperature with reference temp of 100K
        Scalar phi = temperature/100;

        Scalar c_p = 0.661738E+01
                -0.105885E+01 * phi
                +0.201650E+00 * std::pow(phi,2)
                -0.196930E-01 * std::pow(phi,3)
                +0.106460E-02 * std::pow(phi,4)
                -0.303284E-04 * std::pow(phi,5)
                +0.355861E-06 * std::pow(phi,6);
        c_p +=   -0.549169E+01 * std::pow(phi,-1)
                +0.585171E+01* std::pow(phi,-2)
                -0.372865E+01* std::pow(phi,-3)
                +0.133981E+01* std::pow(phi,-4)
                -0.233758E+00* std::pow(phi,-5)
                +0.125718E-01* std::pow(phi,-6);
        c_p *= IdealGas::R / (molarMass() * 1000); // in J/mol/K * mol / kg / 1000 = kJ/kg/K


        return  c_p;
    }

    /*!
     * \brief Thermal conductivity \f$\mathrm{[[W/(m*K)]}\f$ of air.
     *
     * Isobaric Properties for Nitrogen in: NIST Standard \cite NIST <BR>
     * evaluated at p=.1 MPa, T=20°C <BR>
     * Nitrogen: 0.025398 <BR>
     * Oxygen: 0.026105 <BR>
     * lambda_air is approximately 0.78*lambda_N2+0.22*lambda_O2
     *
     * \param temperature absolute temperature in \f$\mathrm{[K]}\f$
     * \param pressure of the phase in \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasThermalConductivity(Scalar temperature, Scalar pressure)
    {
        return 0.0255535;
    }
};

} // end namespace

#endif
