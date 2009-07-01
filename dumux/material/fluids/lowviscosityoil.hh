low// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Jochen Fritz, Andreas Lauser, Markus Wolff   *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_FLUID_LOW_VISCOSITY_OIL_HH
#define DUNE_FLUID_LOW_VISCOSITY_OIL_HH

#include <dumux/material/property_baseclasses.hh>

#include <dune/common/exceptions.hh>

namespace Dune
{
/*!
 * \brief LNAPL with a low viscosity.
 */
class LowViscosityOil : public Fluid
{
public:
    LowViscosityOil()
    {}

    double viscosity ( double T=283.15, double p=1e5, double X=1.0) const
    {
        return 1.1e-3;//[kg/(ms)]
    }

    double density ( double T=283.15, double p=1e5, double X=1.0) const
    {
        if (constDensity_)
            return constDensity_;
        else
            return 890.0; // [kg/m^3]
    }
    double enthalpy (double T=283.15, double p=1e5, double X = 1) const
    {
        if (constEnthalpy_)
            return constEnthalpy_;
        else {
            //            return constRelOil.enthalpy(T,p,X);
            // TODO
            DUNE_THROW(Dune::NotImplemented, "Non-constant enthalpy of oil");
        }
    }

    double intEnergy( double T=283.15, double p=1e5, double X = 1) const
    {
        double u;
        double rho_mass = density(T,p);
        double h = enthalpy(T,p);

        u = h - (p / rho_mass);
        return u;
    }

private:
    double constDensity_;
    double constViscosity_;
    double constEnthalpy_;
};
}

#endif

