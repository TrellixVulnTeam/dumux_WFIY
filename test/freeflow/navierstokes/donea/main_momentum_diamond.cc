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

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtk/function.hh>
#include <dumux/io/vtk/intersectionwriter.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem_diamond.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::DoneaTestNewMomentum;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using GridManager = Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x.resize(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // use the staggered FV assembler
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // linearize & solve
    nonLinearSolver.solve(x);

    std::vector<Dune::FieldVector<double,2>> velocity(gridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double,2>> faceVelocityVector(x.size());
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(gridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, x);
        const auto eIdx = gridGeometry->elementMapper().index(element);

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

    Dune::VTKWriter<std::decay_t<decltype(gridGeometry->gridView())>> writer(gridGeometry->gridView());
    using Field = Vtk::template Field<std::decay_t<decltype(gridGeometry->gridView())>>;
    writer.addCellData(Field(gridGeometry->gridView(), gridGeometry->elementMapper(), velocity,
                        "velocity",
                        /*numComp*/2, /*codim*/0).get());

    std::vector<int> rank;
    rank.resize(gridGeometry->gridView().size(0));

    for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
    {
        const auto eIdxGlobal = gridGeometry->elementMapper().index(element);
        rank[eIdxGlobal] = gridGeometry->gridView().comm().rank();
    }

    writer.addCellData(rank, "rank");
    writer.write("donea_momentum_diamond");

    ConformingIntersectionWriter faceVtk(gridGeometry->gridView());

    std::vector<std::size_t> dofIdx(x.size());
    for (const auto& facet : facets(gridGeometry->gridView()))
    {
        const auto idx = gridGeometry->gridView().indexSet().index(facet);
        dofIdx[idx] = idx;
    }
    faceVtk.addField(dofIdx, "dofIdx");
    faceVtk.addField(faceVelocityVector, "velocityVector");

    auto partionType = [&](const auto& is, const auto idx)
    {
        const auto& facet = is.inside().template subEntity <1> (is.indexInInside());
        return facet.partitionType();
    };

    faceVtk.addField(partionType, "partitionType");
    faceVtk.write("facedata_" + std::to_string(gridGeometry->gridView().comm().rank()), Dune::VTK::ascii);

    if (getParam<bool>("Problem.PrintL2Error", true))
    {
        const auto velocityL2error = problem->calculateL2Error(x);

        std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                        << std::setw(6) << gridGeometry->numDofs()  << " face dofs " << std::scientific
                        << " , L2 = " << velocityL2error
                        << std::endl;
    }


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
