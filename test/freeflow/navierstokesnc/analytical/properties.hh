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
 * \ingroup NavierStokesNCTests
 * \brief The properties of the test for the staggered grid compositional Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_ANALYTICAL_NC_TEST_PROPERTIES_HH
#define DUMUX_ANALYTICAL_NC_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif


#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct AnalyticalNCTest {};
struct AnalyticalNCTestMomentum { using InheritsFrom = std::tuple<AnalyticalNCTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct AnalyticalNCTestMass {using InheritsFrom = std::tuple<AnalyticalNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::AnalyticalNCTest>
{
    using Traits = MultiDomainTraits<TTag::AnalyticalNCTestMomentum, TTag::AnalyticalNCTestMass>;
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

//! A simple fluid system with one tracer component
template<class TypeTag>
class ConstantFluidSystem : public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, ConstantFluidSystem<TypeTag>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Base = Dumux::FluidSystems::Base<Scalar, ConstantFluidSystem<TypeTag>>;

public:
    //! The number of components
    static constexpr int numComponents = 2;
    static constexpr int numPhases = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "component_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "base"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 1.0; }

    using Base::molarDensity;
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    { return 1.0; }

    using Base::binaryDiffusionCoefficient;
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState, int phaseIdx,
                                             int compIIdx, int compJIdx)
    {
        static const Scalar diffCoeff = getParam<Scalar>("Component.DiffusionCoefficient", 1e-6);
        return diffCoeff;
    }

    using Base::density;
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx = 0)
    {
        static const Scalar density = getParam<Scalar>("Component.LiquidDensity", 1.0);
        return density;
    }

    using Base::viscosity;
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const int phaseIdx)
    {
        static const Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar rho = density(fluidState,phaseIdx);
        return rho * nu;
    }

    static constexpr bool isCompressible(int phaseIdx)
    { return true; }

    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return true; }

    static constexpr bool isMiscible()
    { return false; }

    static void init()
    {}
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::AnalyticalNCTest> { static constexpr int value = 0; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::AnalyticalNCTest> { using type = ConstantFluidSystem<TypeTag>; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::AnalyticalNCTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::AnalyticalNCTest> { using type = AnalyticalNCTestProblem<TypeTag> ; };

// Use mole fraction formulation
template<class TypeTag>
struct UseMoles<TypeTag, TTag::AnalyticalNCTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::AnalyticalNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::AnalyticalNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::AnalyticalNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFaceVariablesCache<TypeTag, TTag::AnalyticalNCTest> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
