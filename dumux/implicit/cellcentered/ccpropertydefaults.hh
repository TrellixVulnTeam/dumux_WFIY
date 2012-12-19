// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2010 by Bernd Flemisch                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup CCProperties
 * \ingroup CCModel
 * \file
 *
 * \brief Default properties for box models
 */
#ifndef DUMUX_CC_PROPERTY_DEFAULTS_HH
#define DUMUX_CC_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/common/implicitpropertydefaults.hh>
#include "ccassembler.hh"
#include "ccmodel.hh"
#include "ccfvelementgeometry.hh"
#include "ccelementboundarytypes.hh"
#include "cclocalresidual.hh"
#include "ccelementvolumevariables.hh"
#include "ccproperties.hh"

namespace Dumux {

// forward declaration
template<class TypeTag>
class CCModel;

namespace Properties {
//! Set the default for the FVElementGeometry
SET_TYPE_PROP(CCModel, FVElementGeometry, Dumux::CCFVElementGeometry<TypeTag>);

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCModel, ElementBoundaryTypes, Dumux::CCElementBoundaryTypes<TypeTag>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(CCModel, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCModel, BaseLocalResidual, Dumux::CCLocalResidual<TypeTag>);

//! Set the BaseModel to CCModel
SET_TYPE_PROP(CCModel, BaseModel, Dumux::CCModel<TypeTag>);

//! An array of secondary variable containers
SET_TYPE_PROP(CCModel, ElementVolumeVariables, Dumux::CCElementVolumeVariables<TypeTag>);

//! Assembler for the global jacobian matrix
SET_TYPE_PROP(CCModel, JacobianAssembler, Dumux::CCAssembler<TypeTag>);

//! indicate that this is no box discretization
SET_BOOL_PROP(CCModel, ImplicitIsBox, false);

} // namespace Properties
} // namespace Dumux

#endif
