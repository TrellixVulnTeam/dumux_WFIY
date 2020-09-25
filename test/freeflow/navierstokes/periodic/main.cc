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
 * \ingroup NavierStokesTests
 * \brief Stationary test for the staggered grid Navier-Stokes model (Kovasznay 1948, \cite Kovasznay1948).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtkfunction.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/facecentered/staggered/fvgridgeometry.hh>

#include <dumux/assembly/staggeredfvassembler.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <dumux/io/vtk/intersectionwriter.hh>

#include "../l2error.hh"
#include "../analyticalsolution.hh"
#include "problem.hh"

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PeriodicTest>
{
private:
    using Traits = MultiDomainTraits<TTag::PeriodicTestMomentum, TTag::PeriodicTestMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::PeriodicTestMomentum;
    using MassTypeTag = Properties::TTag::PeriodicTestMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    momentumGridGeometry->update();

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    massGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = Dune::index_constant<0>();
    constexpr auto massIdx = Dune::index_constant<1>();
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());


    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);


    // intialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    // const auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
    // const auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
    // vtkWriter.addField(exactPressure, "pressureExact");
    // vtkWriter.addField(exactVelocity, "velocityExact");

    ConformingIntersectionWriter faceVtk(momentumGridGeometry->gridView());

    std::vector<std::size_t> dofIdx(x[momentumIdx].size());
        for (const auto& facet : facets(momentumGridGeometry->gridView()))
    {
        const auto idx = momentumGridGeometry->gridView().indexSet().index(facet);
        dofIdx[idx] = idx;
    }
    faceVtk.addField(dofIdx, "dofIdx");



    // faceVtk.write("facedata_" + std::to_string(gridGeometry->gridView().comm().rank()), Dune::VTK::ascii);

    // the linear solver
    // using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>>
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    std::vector<Dune::FieldVector<double,2>> faceVelocityVector(x[momentumIdx].size());
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        auto fvGeometry = localView(*momentumGridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(momentumGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, x);

        for (const auto& scv : scvs(fvGeometry))
            faceVelocityVector[scv.dofIndex()][scv.directionIndex()] = elemVolVars[scv].velocity();
    }

    faceVtk.addField(faceVelocityVector, "velocityVector");

    faceVtk.write("facedata", Dune::VTK::ascii);

    // write vtk output
    vtkWriter.write(1.0);
    timer.stop();

    // if (getParam<bool>("Problem.PrintL2Error"))
    // {
    //     auto pressureL2error = calculateL2Error(*massProblem, x[massIdx]);
    //     auto velocityL2error = calculateL2Error(*momentumProblem, x[momentumIdx]);

    //     std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
    //                     << std::setw(6) << massGridGeometry->numDofs() << " cc dofs and " << momentumGridGeometry->numDofs()
    //                     << " face dofs (total: " << massGridGeometry->numDofs() + momentumGridGeometry->numDofs() << "): "
    //                     << std::scientific
    //                     << "L2(p) = " << pressureL2error.absolute[0] << " / " << pressureL2error.relative[0]
    //                     << " , L2(vx) = " << velocityL2error.absolute[0] << " / " << velocityL2error.relative[0]
    //                     << " , L2(vy) = " << velocityL2error.absolute[1] << " / " << velocityL2error.relative[1]
    //                     << std::endl;
    // }

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

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
