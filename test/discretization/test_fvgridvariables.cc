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
 * \brief Test for the new experimental finite volume grid variables.
 */
#include <config.h>
#include <iostream>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

// we use the 1p type tag here in order not to be obliged
// to define grid flux vars cache & vol vars cache...
#include <dumux/common/properties.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/fvgridvariables.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/spatialparams/fv1pconstant.hh>

namespace Dumux::Properties {

// new type tags
namespace TTag {
struct GridVariablesTest { using InheritsFrom = std::tuple<OneP>; };
struct GridVariablesTestBox { using InheritsFrom = std::tuple<GridVariablesTest, BoxModel>; };
} // end namespace TTag

// property definitions
template<class TypeTag>
struct Grid<TypeTag, TTag::GridVariablesTest>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::GridVariablesTest>
{ using type = FVProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::GridVariablesTest>
{
private:
    using Scalar = GetPropType<TypeTag, Scalar>;
    using GG = GetPropType<TypeTag, GridGeometry>;
public:
    using type = FVSpatialParamsOnePConstant<GG, Scalar>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::GridVariablesTest>
{
private:
    using Scalar = GetPropType<TypeTag, Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar>>;
};

// use the new experimental grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::GridVariablesTest>
{
private:
    using GVV = GetPropType<TypeTag, GridVolumeVariables>;
    using GFC = GetPropType<TypeTag, GridFluxVariablesCache>;
    using X = GetPropType<TypeTag, SolutionVector>;

public:
    using type = Dumux::Experimental::FVGridVariables<GVV, GFC, X>;
};

} // end namespace Dumux::Properties

int main (int argc, char *argv[])
{
    Dune::MPIHelper::instance(argc, argv);

    using namespace Dumux;
    Dumux::Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::GridVariablesTestBox;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridFactory = Dune::StructuredGridFactory<Grid>;
    auto grid = GridFactory::createCubeGrid({0.0, 0.0}, {1.0, 1.0}, {2, 2});

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid->leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // This constructor should fail as the problem does not implement initial()
    bool caught = false;
    try { auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry); }
    catch (...) { caught = true; }
    if (!caught)
        DUNE_THROW(Dune::Exception, "Expected construction to fail");

    // Construction with a solution
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x; x.resize(gridGeometry->numDofs()); x = 0.0;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry, x);

    // Construction from a moved solution (TODO: how to check if move succeeded?)
    gridVariables = std::make_shared<GridVariables>(problem, gridGeometry, std::move(x));

    // Construction from initializer lambda
    const auto init = [gridGeometry] (auto& x) { x.resize(gridGeometry->numDofs()); x = 2.25; };
    gridVariables = std::make_shared<GridVariables>(problem, gridGeometry, init);

    const auto& dofs = gridVariables->dofs();
    if (std::any_of(dofs.begin(), dofs.end(),
                    [] (auto d) { return Dune::FloatCmp::ne(2.25, d[0]); }))
        DUNE_THROW(Dune::Exception, "Unexpected dof value");

    std::cout << "\nAll tests passed" << std::endl;
    return 0;
}
