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
 * \ingroup CCMpfaDiscretization
 * \brief Properties for all models using cell-centered finite volume scheme with mpfa
 * \note Inherit from these properties to use a cell-centered finite volume scheme with mpfa
 */
#ifndef DUMUX_CC_MPFA_PROPERTIES_HH
#define DUMUX_CC_MPFA_PROPERTIES_HH

#include <dune/common/reservedvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/assembly/cclocalresidual.hh>

#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/discretization/cellcentered/elementboundarytypes.hh>

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometrytraits.hh>
#include <dumux/discretization/cellcentered/mpfa/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/gridvolumevariables.hh>
#include <dumux/discretization/cellcentered/mpfa/gridfluxvariablescache.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumedatahandle.hh>
#include <dumux/discretization/cellcentered/mpfa/fluxvariablescachefiller.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cell-centered mpfa scheme.
NEW_TYPE_TAG(CCMpfaModel, INHERITS_FROM(FiniteVolumeModel));

//! Set the index set type used on the dual grid nodes
SET_PROP(CCMpfaModel, DualGridNodalIndexSet)
{
private:
    using GV = typename GET_PROP_TYPE(TypeTag, GridView);
    using Traits = NodalIndexSetDefaultTraits< GV >;

public:
    using type = CCMpfaDualGridNodalIndexSet< Traits >;
};

//! Per default, we use the dynamic mpfa-o interaction volume
SET_PROP(CCMpfaModel, PrimaryInteractionVolume)
{
public:
    using type = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
};

//! Per default, we use the dynamic mpfa-o interaction volume on boundaries
SET_PROP(CCMpfaModel, SecondaryInteractionVolume)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);

    // use the default traits
    using Traits = CCMpfaODefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOInteractionVolume< Traits >;
};

//! Set the default for the global finite volume geometry
SET_PROP(CCMpfaModel, FVGridGeometry)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryIV = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using SecondaryIV = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);
    using Traits = CCMpfaFVGridGeometryTraits<GridView, NodalIndexSet, PrimaryIV, SecondaryIV>;
public:
    using type = CCMpfaFVGridGeometry<GridView, Traits, GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache)>;
};

//! The grid volume variables vector class
SET_PROP(CCMpfaModel, GridVolumeVariables)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableGridVolumeVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
public:
    using type = CCMpfaGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};

//! The grid volume variables vector class
SET_PROP(CCMpfaModel, GridFluxVariablesCache)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableGridFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;

    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);

    using PhysicsTraits = IvDataHandlePhysicsTraits<typename GET_PROP_TYPE(TypeTag, ModelTraits)>;
    using PrimaryMatVecTraits = typename PrimaryInteractionVolume::Traits::MatVecTraits;
    using SecondaryMatVecTraits = typename SecondaryInteractionVolume::Traits::MatVecTraits;

    using PrimaryIvDataHandle = InteractionVolumeDataHandle<PrimaryMatVecTraits, PhysicsTraits>;
    using SecondaryIvDataHandle = InteractionVolumeDataHandle<SecondaryMatVecTraits, PhysicsTraits>;

    using Traits = CCMpfaDefaultGridFluxVariablesCacheTraits<Problem,
                                                             FluxVariablesCache, FluxVariablesCacheFiller,
                                                             PrimaryInteractionVolume, SecondaryInteractionVolume,
                                                             PrimaryIvDataHandle, SecondaryIvDataHandle>;
public:
    using type = CCMpfaGridFluxVariablesCache<Traits, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
SET_TYPE_PROP(CCMpfaModel, ElementBoundaryTypes, CCElementBoundaryTypes);

//! Set the BaseLocalResidual to CCLocalResidual
SET_TYPE_PROP(CCMpfaModel, BaseLocalResidual, CCLocalResidual<TypeTag>);
} // namespace Properties
} // namespace Dumux

#endif
