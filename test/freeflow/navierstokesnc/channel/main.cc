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
 * \brief Test for the staggered grid multi-component (Navier-)Stokes model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
//#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::ChannelNCTestMomentum;
    using MassTypeTag = Properties::TTag::ChannelNCTestMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MassTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    //using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    //auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    //using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    //using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;
    //auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    //using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    //auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    //constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    //constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    //using SolutionVector = typename Traits::SolutionVector;
    using SolutionVector = GetPropType<MassTypeTag, Properties::SolutionVector>;
    SolutionVector x;
    //x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x.resize(massGridGeometry->numDofs());
    //momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x);
    auto xOld = x;

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    massProblem->setTimeLoop(timeLoop);
    //momentumProblem->setTimeLoop(timeLoop);

    // the grid variables
    //using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    //auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // read mass data to init dirichlet for pressure
    std::ifstream fin;
    int timestep = 0;
    auto oldStaggeredVelocity = Dune::FieldVector<Scalar, 5000>();
    auto oldStaggeredMass = Dune::FieldVector<Dune::FieldVector<Scalar,2>,1250>();
    Scalar value = 0.0;
    auto scvIdx = 0;
    fin.open("old_mass_" + std::to_string(timestep) + ".txt");
    while (fin >> value)
    {
        oldStaggeredMass[scvIdx][0] = value;
        fin >> value;
        oldStaggeredMass[scvIdx][1] = value;
        ++scvIdx;
    }
    massProblem->setOldStaggered(oldStaggeredVelocity, oldStaggeredMass);
    fin.close();

    //couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    //momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x, massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    //vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.addVolumeVariable([](const auto& volVars){ return volVars.pressure() - 1.1e5; }, "deltaP");
    vtkWriter.write(0);

    //using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    using Assembler = FVAssembler<MassTypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(massProblem, massGridGeometry, massGridVariables,
            timeLoop, xOld);
    //auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
    //                                             std::make_tuple(momentumGridGeometry, massGridGeometry),
    //                                             std::make_tuple(momentumGridVariables, massGridVariables),
    //                                             couplingManager,
    //                                             timeLoop, xOld);
    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // read velocity data from old-staggered simulation
        int scvfIdx = 0;
        fin.open("old_velocity_" + std::to_string(timestep) + ".txt");
        while (fin >> value)
        {
            oldStaggeredVelocity[scvfIdx] = value;
            ++scvfIdx;
        }
        fin.close();
        scvIdx = 0;
        fin.open("old_mass_" + std::to_string(timestep) + ".txt");
        while (fin >> value)
        {
            oldStaggeredMass[scvIdx][0] = value;
            fin >> value;
            oldStaggeredMass[scvIdx][1] = value;
            ++scvIdx;
        }
        massProblem->setOldStaggered(oldStaggeredVelocity, oldStaggeredMass);
        fin.close();
        ++timestep;
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        //momentumGridVariables->advanceTimeStep();
        massGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

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
} // end main
