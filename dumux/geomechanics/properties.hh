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
 * \ingroup Properties
 * \ingroup Geomechanics
 * \brief Defines a type tag and some properties for geomechanical DuMuX models.
 */

#ifndef DUMUX_GEOMECHANICS_PROPERTIES_HH
#define DUMUX_GEOMECHANICS_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/solidstates/inertsolidstate.hh>
#include <dumux/material/solidsystems/inertsolidphase.hh>
#include <dumux/discretization/hookeslaw.hh>

#include "stressvariablescache.hh"
#include "velocityoutput.hh"

namespace Dumux {
namespace Properties {

//! Type tag for geomechanical models
NEW_TYPE_TAG(Geomechanics, INHERITS_FROM(ModelProperties));

//! The flux variables cache class for models involving flow in porous media
SET_TYPE_PROP(Geomechanics, FluxVariablesCache, StressVariablesCache< typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                      typename GET_PROP_TYPE(TypeTag, FVGridGeometry) >);

//! The (currently empty) velocity output
SET_TYPE_PROP(Geomechanics, VelocityOutput, GeomechanicsVelocityOutput<typename GET_PROP_TYPE(TypeTag, GridVariables)>);

//! The solid state must be inert
SET_PROP(Geomechanics, SolidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SolidSystem = typename GET_PROP_TYPE(TypeTag, SolidSystem);
public:
    using type = InertSolidState<Scalar, SolidSystem>;
};

//! Per default we use one constant component in the inert solid system
SET_PROP(Geomechanics, SolidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using InertComponent = Components::Constant<1, Scalar>;
public:
    using type = SolidSystems::InertSolidPhase<Scalar, InertComponent>;
};
} // namespace Properties
} // namespace Dumux

 #endif
