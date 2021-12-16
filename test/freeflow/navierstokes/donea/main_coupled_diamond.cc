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
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>

// #include "../l2error.hh"
// #include "../analyticalsolution.hh"
#include "problem_diamond.hh"

namespace Dumux::Properties {

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DoneaTestNew>
{
private:
    using Traits = MultiDomainTraits<TTag::DoneaTestNewMomentum, TTag::DoneaTestNewMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::DoneaTestNewMomentum;
    using MassTypeTag = Properties::TTag::DoneaTestNewMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using GridManager = Dumux::GridManager<GetPropType<MomentumTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    Dune::Timer timer;
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = GetPropType<MassTypeTag, Properties::CouplingManager>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
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

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // intialize the vtk output module
    //using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    //VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    //IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    //vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    // const auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
    // const auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
    // vtkWriter.addField(exactPressure, "pressureExact");
    // vtkWriter.addField(exactVelocity, "velocityExact");

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    nonLinearSolver.solve(x);

    std::vector<Dune::FieldVector<double,2>> velocity(momentumGridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double,2>> faceVelocityVector(x[momentumIdx].size());
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        auto fvGeometry = localView(*momentumGridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(momentumGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, x[momentumIdx]);
        const auto eIdx = momentumGridGeometry->elementMapper().index(element);

        const auto geometry = element.geometry();

        const auto& localBasis = fvGeometry.feLocalBasis();
        using ShapeValue = typename Dune::FieldVector<double, 1>;
        std::vector<ShapeValue> shapeValues;

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(geometry.center());
        localBasis.evaluateFunction(ipLocal, shapeValues);

        for (const auto& scv : scvs(fvGeometry))
        {
            velocity[eIdx] += shapeValues[scv.indexInElement()][0]*elemVolVars[scv].velocity();
            faceVelocityVector[scv.dofIndex()] = elemVolVars[scv].velocity();
        }
    }

    Dune::VTKWriter<std::decay_t<decltype(momentumGridGeometry->gridView())>> writer(momentumGridGeometry->gridView());
    using Field = Vtk::template Field<std::decay_t<decltype(momentumGridGeometry->gridView())>>;
    writer.addCellData(Field(momentumGridGeometry->gridView(), momentumGridGeometry->elementMapper(), velocity,
                             "velocity",
                            /*numComp*/2, /*codim*/0).get());

    writer.addVertexData(x[massIdx],  "pressure");

    // std::vector<int> rank;
    // rank.resize(momentumGridGeometry->gridView().size(0));

    // for (const auto& element : elements(momentumGridGeometry->gridView(), Dune::Partitions::interior))
    // {
    //     const auto eIdxGlobal = momentumGridGeometry->elementMapper().index(element);
    //     rank[eIdxGlobal] = momentumGridGeometry->gridView().comm().rank();
    // }

    // writer.addCellData(rank, "rank");

    writer.write("donea_new_momentum_coupled_diamond");

    //vtkWriter.write(1.0);

    // using GridView = std::decay_t<decltype(momentumGridGeometry->gridView())>;
    // using SolVec = std::decay_t<decltype(x[momentumIdx])>;

    // class VelocityFunction
    // {
    //     public:
    //     typedef Dune::VTK::SkeletonFunctionTraits<std::decay_t<decltype(momentumGridGeometry->gridView())>, typename std::decay_t<decltype(momentumGridGeometry->gridView())>::ctype>
    //     Traits;

    //     VelocityFunction(const GridView& gv, const SolVec& s) : gv_(gv), s_(s) {}

    //     //! return number of components
    //     unsigned dimRange() const { return 1; }

    //     void evaluate(const typename Traits::Cell& c,
    //                   const typename Traits::Domain& xl,
    //                   typename Traits::Range& result) const
    //     {
    //         assert(c.conforming());


    //         const auto globalIdx = gv_.indexSet().subIndex(c.inside(), c.indexInInside(), 1);
    //         result.resize(2, s_[globalIdx]);
    //     }

    //     private:
    //     const GridView gv_;
    //     const SolVec& s_;
    // };

    // auto faceData = std::make_shared<VelocityFunction>(momentumGridGeometry->gridView(), x[momentumIdx]);
    ConformingIntersectionWriter vtk(momentumGridGeometry->gridView());
    // vtk.addCellData(faceData, "velcocityScalar");
    vtk.write("facedata", Dune::VTK::ascii);

    // if (getParam<bool>("Problem.PrintL2Error", true))
    // {
    //     const auto pressureL2error = calculateL2Error(*massProblem, x[massIdx]);
    //     const auto velocityL2error = momentumProblem->calculateL2Error(x[momentumIdx]);

    //     std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
    //                     << std::setw(6) << massGridGeometry->numDofs() << " cc dofs and " << momentumGridGeometry->numDofs()
    //                     << " face dofs (total: " << massGridGeometry->numDofs() + momentumGridGeometry->numDofs() << "): "
    //                     << std::scientific
    //                     << "L2(p) = " << pressureL2error.absolute[0] << " / " << pressureL2error.relative[0]
    //                     << " , L2(v) = " << velocityL2error
    //                     << std::endl;
    // }

    timer.stop();

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
        // Parameters::print();
        // DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
