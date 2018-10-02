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
 *
 * \brief Test for the equal dimension boundary coupling model
 */
#include <config.h>

#include <ctime>
#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/subgridgridcreator.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/boundary/darcydarcy/couplingmanager.hh>

#include "problem.hh"

#ifndef DOMAINSPLIT
#define DOMAINSPLIT 0
#endif

namespace Dumux {
namespace Properties {

NEW_TYPE_TAG(OnePSubTypeTag, INHERITS_FROM(CCTpfaModel, OneP));
// differentiate between the two subproblems
NEW_TYPE_TAG(OnePSubTypeTag0, INHERITS_FROM(OnePSubTypeTag));
NEW_TYPE_TAG(OnePSubTypeTag1, INHERITS_FROM(OnePSubTypeTag));

// the coupling manager
SET_TYPE_PROP(OnePSubTypeTag, CouplingManager,
              DarcyDarcyBoundaryCouplingManager<MultiDomainTraits<TTAG(OnePSubTypeTag0), TTAG(OnePSubTypeTag1)>>);

// Set the grid type
#if DOMAINSPLIT==1
SET_PROP(OnePSubTypeTag, Grid)
{
    using FullDomainGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
    using type = Dune::SubGrid<FullDomainGrid::dimension, FullDomainGrid>;
};
#elif DOMAINSPLIT==0
SET_PROP(OnePSubTypeTag, Grid)
{
    using FullDomainGrid = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>;
    using type = FullDomainGrid;
};
#endif

// set the spatial params
SET_TYPE_PROP(OnePSubTypeTag, SpatialParams, OnePTestSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                                   typename GET_PROP_TYPE(TypeTag, Scalar)>);

// the fluid system
SET_PROP(OnePSubTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// differentiate between the two subproblems
SET_TYPE_PROP(OnePSubTypeTag0, Problem, OnePTestProblem<TypeTag, 0>);
SET_TYPE_PROP(OnePSubTypeTag1, Problem, OnePTestProblem<TypeTag, 1>);

} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using TypeTag = TTAG(OnePSubTypeTag);
    using SubTypeTag0 = TTAG(OnePSubTypeTag0);
    using SubTypeTag1 = TTAG(OnePSubTypeTag1);

    // create the full grid that we are gonna split for the output
    using FullDomainGrid = typename GET_PROP(TypeTag, Grid)::FullDomainGrid;
    GridManager<FullDomainGrid> gridManager;
    gridManager.init();

    // get the lens
    using GlobalPosition = typename FullDomainGrid::template Codim<0>::Geometry::GlobalCoordinate;
    const auto lensLowerLeft = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
    const auto lensUpperRight = getParam<GlobalPosition>("SpatialParams.LensUpperRight");

///////////////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for lens and the rest (using sub-grid)
///////////////////////////////////////////////////////////////////////////////////////////
#if DOMAINSPLIT==1

    // create subgrids coinciding with the lens
    auto elementSelector0 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };
    auto elementSelector1 = [&lensLowerLeft, &lensUpperRight](const auto& element)
    { return !LensSpatialParams::pointInLens(element.geometry().center(), lensLowerLeft, lensUpperRight); };

    auto subGrid0 = SubgridGridCreator<FullDomainGrid>::makeGrid(gridManager.grid(), elementSelector0);
    auto subGrid1 = SubgridGridCreator<FullDomainGrid>::makeGrid(gridManager.grid(), elementSelector1);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& gridView0 = subGrid0->leafGridView();
    const auto& gridView1 = subGrid1->leafGridView();

//////////////////////////////////////////////////////////////////////////////////
// split the domains by creating two separate grids for the lower and upper half
//////////////////////////////////////////////////////////////////////////////////
#elif DOMAINSPLIT==0

    // create an upper half and lower half grid
    GridManager<typename GET_PROP_TYPE(TypeTag, Grid)> gridManager0, gridManager1;
    gridManager0.init("1");
    gridManager1.init("2");

    // we compute on the leaf grid views
    const auto& gridView0 = gridManager0.grid().leafGridView();
    const auto& gridView1 = gridManager1.grid().leafGridView();

#endif

    ////////////////////////////////////////////////
    // run the multidomain simulation on two grids
    ////////////////////////////////////////////////

    // create the finite volume grid geometries
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry0 = std::make_shared<FVGridGeometry>(gridView0);
    auto fvGridGeometry1 = std::make_shared<FVGridGeometry>(gridView1);
    fvGridGeometry0->update();
    fvGridGeometry1->update();

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<SubTypeTag0, SubTypeTag1>;
    constexpr auto domain0Idx = Traits::template DomainIdx<0>();
    constexpr auto domain1Idx = Traits::template DomainIdx<1>();

    // the coupling manager
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (initial and boundary conditions)
    using Problem0 = typename GET_PROP_TYPE(SubTypeTag0, Problem);
    auto problem0 = std::make_shared<Problem0>(fvGridGeometry0, couplingManager, "1");
    problem0->spatialParams().setLens(lensLowerLeft, lensUpperRight);
    using Problem1 = typename GET_PROP_TYPE(SubTypeTag1, Problem);
    auto problem1 = std::make_shared<Problem1>(fvGridGeometry1, couplingManager, "2");
    problem1->spatialParams().setLens(lensLowerLeft, lensUpperRight);

    // the solution vector
    Traits::SolutionVector sol;
    sol[domain0Idx].resize(fvGridGeometry0->numDofs());
    sol[domain1Idx].resize(fvGridGeometry1->numDofs());
    problem0->applyInitialSolution(sol[domain0Idx]);
    problem1->applyInitialSolution(sol[domain1Idx]);
    auto oldSol = sol;

    // compute the coupling map
    couplingManager->init(problem0, problem1, sol);

    // the grid variables
    using GridVariables0 = typename GET_PROP_TYPE(SubTypeTag0, GridVariables);
    auto gridVariables0 = std::make_shared<GridVariables0>(problem0, fvGridGeometry0);
    using GridVariables1 = typename GET_PROP_TYPE(SubTypeTag1, GridVariables);
    auto gridVariables1 = std::make_shared<GridVariables1>(problem1, fvGridGeometry1);
    gridVariables0->init(sol[domain0Idx], oldSol[domain0Idx]);
    gridVariables1->init(sol[domain1Idx], oldSol[domain1Idx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    using SolutionVector0 = std::decay_t<decltype(sol[domain0Idx])>;
    VtkOutputModule<GridVariables0, SolutionVector0> vtkWriter0(*gridVariables0, sol[domain0Idx], problem0->name());
    GET_PROP_TYPE(SubTypeTag0, VtkOutputFields)::init(vtkWriter0);
    vtkWriter0.write(0.0);

    using SolutionVector1 = std::decay_t<decltype(sol[domain1Idx])>;
    VtkOutputModule<GridVariables1, SolutionVector1> vtkWriter1(*gridVariables1, sol[domain1Idx], problem1->name());
    GET_PROP_TYPE(SubTypeTag1, VtkOutputFields)::init(vtkWriter1);
    vtkWriter1.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(problem0, problem1),
                                                 std::make_tuple(fvGridGeometry0, fvGridGeometry1),
                                                 std::make_tuple(gridVariables0, gridVariables1),
                                                 couplingManager, timeLoop);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(oldSol);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables0->advanceTimeStep();
        gridVariables1->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter0.write(timeLoop->time());
        vtkWriter1.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(mpiHelper.getCollectiveCommunication());

    //////////////////////////////////////////////////////////////////////////
    // write out a combined vtu for vtu comparison in the testing framework
    //////////////////////////////////////////////////////////////////////////

    const auto& gridView = gridManager.grid().leafGridView();
    CCTpfaFVGridGeometry<typename FullDomainGrid::LeafGridView> fvGridGeometry(gridView);
    const auto& bBoxTree = fvGridGeometry.boundingBoxTree();
    // copy data from the subdomains to full domain data vectors
    std::vector<int> processRank(gridView.size(0), 0); // sequential simulation
    std::vector<int> pressure(gridView.size(0));

    for (const auto& element : elements(gridView0))
    {
        const auto eIdx = fvGridGeometry0->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain0Idx][eIdx][0];
    }

    for (const auto& element : elements(gridView1))
    {
        const auto eIdx = fvGridGeometry1->elementMapper().index(element);
        const auto eIdxHost = intersectingEntities(element.geometry().center(), bBoxTree)[0];
        pressure[eIdxHost] = sol[domain1Idx][eIdx][0];
    }

    Dune::VTKWriter<typename FullDomainGrid::LeafGridView> vtkWriter(gridView);
    vtkWriter.addCellData(processRank, "process rank");
    vtkWriter.addCellData(pressure, "pressure");
    const auto filename = getParam<std::string>("Problem.Name") + "_combined";
    vtkWriter.write(filename);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}

// error handler
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
