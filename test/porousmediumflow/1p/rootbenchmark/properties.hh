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
 * \ingroup OnePTests
 * \brief Root benchmark
 */

#ifndef DUMUX_ONEP_ROOT_BENCHMARK_PROPERTIES_HH
#define DUMUX_ONEP_ROOT_BENCHMARK_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/parallelgrid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct RootBenchmark { using InheritsFrom = std::tuple<OneP>; };
struct RootBenchmarkCCTpfa { using InheritsFrom = std::tuple<RootBenchmark, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RootBenchmark>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    //using type = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<Scalar, 1>>;
    using type = Dune::ParallelGrid<Dune::FoamGrid<1, 1>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RootBenchmark>
{ using type = RootBenchmarkProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RootBenchmark>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = RootBenchmarkSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RootBenchmark>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Dumux::Properties

#endif
