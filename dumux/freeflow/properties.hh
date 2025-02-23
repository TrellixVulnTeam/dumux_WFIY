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
 * \ingroup FreeflowModels
 * \brief Defines a type tag and some properties for free flow models.
 */

#ifndef DUMUX_FREE_FLOW_PROPERTIES_HH
#define DUMUX_FREE_FLOW_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/flux/fourierslaw.hh>

#include "turbulencemodel.hh"
#include "spatialparams.hh"

namespace Dumux {
namespace Properties {

//! Type tag for free-flow models
// Create new type tags
namespace TTag {
struct FreeFlow { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace TTag

//! Use Fourier's Law as default heat conduction type
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::FreeFlow> { using type = FouriersLaw<TypeTag>; };

// Set the default spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::FreeFlow>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FreeFlowDefaultSpatialParams<GridGeometry, Scalar>;
};

} // namespace Properties
} // namespace Dumux

 #endif
