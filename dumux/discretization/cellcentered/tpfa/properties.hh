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
 * \ingroup CCTpfaDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with TPFA
 * \note Inherit from these properties to use a cell-centered finite volume scheme with TPFA
 */

#ifndef DUMUX_CC_TPFA_PROPERTIES_HH
#define DUMUX_CC_TPFA_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/subcontrolvolume.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/tpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/tpfa/fluxvariablescachefiller.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cell-centered tpfa scheme.
NEW_TYPE_TAG(CCTpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the default for the global finite volume geometry
SET_PROP(CCTpfaModel, FVGridGeometry)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = CCTpfaFVGridGeometry<GridView, enableCache>;
};

//! The grid volume variables vector class
SET_PROP(CCTpfaModel, GridVolumeVariables)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    using type = CCTpfaGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};

//! The grid flux variables cache vector class
SET_PROP(CCTpfaModel, GridFluxVariablesCache)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FluxVariablesCacheFiller = CCTpfaFluxVariablesCacheFiller<TypeTag>;
public:
    using type = CCTpfaGridFluxVariablesCache<Problem, FluxVariablesCache, FluxVariablesCacheFiller, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCTpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCTpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);
} // namespace Properties
} // namespace Dumux

#endif
