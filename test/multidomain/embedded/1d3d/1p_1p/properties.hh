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
 * \ingroup EmbeddedTests
 * \brief The properties for the 1p2c tissue problem.
 */
#ifndef DUMUX_TISSUE_PROPERTIES_HH
#define DUMUX_TISSUE_PROPERTIES_HH

// default to tpfa for both domains
#ifndef BULKTYPETAG
#define BULKTYPETAG TissueCC
#endif
#ifndef LOWDIMTYPETAG
#define LOWDIMTYPETAG BloodFlowCC
#endif
#ifndef COUPLINGMODE
#define COUPLINGMODE Embedded1d3dCouplingMode::Average
#endif

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/multidomain/embedded/couplingmanager1d3d.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_tissue.hh"
#include "spatialparams_tissue.hh"
#include "spatialparams_bloodflow.hh"
#include "problem_bloodflow.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct Tissue { using InheritsFrom = std::tuple<OneP>; };
struct TissueCC { using InheritsFrom = std::tuple<Tissue, CCTpfaModel>; };
struct TissueBox { using InheritsFrom = std::tuple<Tissue, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Tissue> { using type = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 3> >; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Tissue> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::Tissue> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::Tissue> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::Tissue> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Tissue> { using type = TissueProblem<TypeTag>; };

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::Tissue> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Tissue>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Tissue>
{
    using type = TissueSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::Scalar>>;
};


// Create new type tags
namespace TTag {
struct BloodFlow { using InheritsFrom = std::tuple<OneP>; };
struct BloodFlowCC { using InheritsFrom = std::tuple<BloodFlow, CCTpfaModel>; };
struct BloodFlowBox { using InheritsFrom = std::tuple<BloodFlow, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BloodFlow> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::BloodFlow> { using type = BloodFlowProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BloodFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the problem property
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BloodFlow> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BloodFlow>
{
    using type = BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                        GetPropType<TypeTag, Properties::Scalar>>;
};


template<class Traits>
using TheCouplingManager = Embedded1d3dCouplingManager<Traits, COUPLINGMODE>;

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BULKTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<TypeTag, Properties::TTag::LOWDIMTYPETAG>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::BULKTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<0>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::BULKTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<0>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::LOWDIMTYPETAG> { using type = TheCouplingManager<MultiDomainTraits<Properties::TTag::BULKTYPETAG, TypeTag>>; };
template<class TypeTag>
struct PointSource<TypeTag, TTag::LOWDIMTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSource<1>; };
template<class TypeTag>
struct PointSourceHelper<TypeTag, TTag::LOWDIMTYPETAG> { using type = typename GetPropType<TypeTag, Properties::CouplingManager>::PointSourceTraits::template PointSourceHelper<1>; };

} // end namespace Dumux::Properties

#endif
