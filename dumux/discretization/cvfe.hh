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
 * \ingroup Discretization
 * \brief Defines a type tag and some properties for models using the cvfe scheme.
 */

#ifndef DUMUX_DISCRETIZTAION_CVFE_HH
#define DUMUX_DISCRETIZTAION_CVFE_HH

#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/assembly/cvfelocalresidual.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/fvproperties.hh>

#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/discretization/cvfe/elementboundarytypes.hh>
#include <dumux/discretization/cvfe/gridfluxvariablescache.hh>
#include <dumux/discretization/cvfe/gridvolumevariables.hh>
#include <dumux/discretization/cvfe/fvgridgeometry.hh>
#include <dumux/discretization/cvfe/fluxvariablescache.hh>

#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {
namespace Properties {

//! Type tag for the cvfe scheme.
// Create new type tags
namespace TTag {
struct CvfeModel { using InheritsFrom = std::tuple<FiniteVolumeModel>; };
} // end namespace TTag

//! Set the default for the grid geometry
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::CvfeModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CvfeFVGridGeometry<Scalar, GridView, enableCache>;
};

//! The grid volume variables vector class
template<class TypeTag>
struct GridVolumeVariables<TypeTag, TTag::CvfeModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
public:
    using type = CvfeGridVolumeVariables<Problem, VolumeVariables, enableCache>;
};

//! The flux variables cache class
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::CvfeModel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = CvfeFluxVariablesCache<Scalar, GridGeometry>;
};

//! The grid flux variables cache vector class
template<class TypeTag>
struct GridFluxVariablesCache<TypeTag, TTag::CvfeModel>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridFluxVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluxVariablesCache = GetPropTypeOr<TypeTag,
        Properties::FluxVariablesCache, FluxVariablesCaching::EmptyCache<Scalar>
    >;
public:
    using type = CvfeGridFluxVariablesCache<Problem, FluxVariablesCache, enableCache>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::CvfeModel>
{
private:
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;
public:
    using type = CvfeElementBoundaryTypes<BoundaryTypes>;
};

//! Set the BaseLocalResidual to CvfeLocalResidual
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::CvfeModel> { using type = CvfeLocalResidual<TypeTag>; };

} // namespace Properties

namespace Detail {

template<class Problem>
struct ProblemTraits<Problem, DiscretizationMethods::Cvfe>
{
private:
    using GG = std::decay_t<decltype(std::declval<Problem>().gridGeometry())>;
    using Element = typename GG::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GG::SubControlVolume;
public:
    using GridGeometry = GG;
    // BoundaryTypes is whatever the problem returns from boundaryTypes(element, scv)
    using BoundaryTypes = std::decay_t<decltype(std::declval<Problem>().boundaryTypes(std::declval<Element>(), std::declval<SubControlVolume>()))>;
};

} // end namespace Detail

} // namespace Dumux

#endif
