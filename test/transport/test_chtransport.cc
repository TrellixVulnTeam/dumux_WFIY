// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Annika Fuchs                                      *
 *                                                                           *
 *   Copyright (C) 2009 by Yufei Cao                                         *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>

#include "dumux/material/fluids/uniform.hh"
#include "dumux/transport/ch/chtransport.hh"
#include "dumux/transport/ch/fractionalw.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

#include "simplenonlinearproblem.hh"

int main(int argc, char** argv)
{
	try{
		// define the problem dimensions
		const int dim=2;

		// time loop parameters
		const double tStart = 0;
		const double tEnd = 1.5e9;
		double dt = tEnd;
		double firstDt = dt;
		double maxDt = dt;
		int modulo = 1;
		double slLengthFactor = 10;

		typedef double Scalar;
		typedef Dune::SGrid<dim,dim> Grid;
		typedef Grid::LeafGridView GridView;
		Dune::FieldVector<int,dim> N(1); N[0] = 64;
		Dune::FieldVector<Scalar,dim> L(0);
		Dune::FieldVector<Scalar,dim> H(300); H[0] = 600;

		Grid grid(N,L,H);
		GridView gridView(grid.leafView());

		Dune::Uniform fluid;
		Dune::HomogeneousNonlinearSoil<Grid, Scalar> soil;
		Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, fluid, fluid);

		double initsat=0;
		Dune::FieldVector<double,dim> velocity(0);
		velocity[0] = 1.0/6.0*1e-6;
		typedef Dune::VariableClass<GridView, Scalar> VariableClass;
		VariableClass variables(gridView, initsat, velocity);

		Dune::SimpleNonlinearProblem<GridView, Scalar, VariableClass> problem(variables, materialLaw, L, H);
		Dune::FractionalW<GridView, Scalar, VariableClass> fractionalW(problem);

		typedef Dune::ChTransport<GridView, Scalar, VariableClass> Transport;
		Transport transport(gridView, problem, fractionalW, slLengthFactor);

		Dune::RungeKuttaStep<Grid, Transport> timeStep(1);
		Dune::TimeLoop<Grid, Transport > timeloop(tStart, tEnd, dt, "chtransport", modulo, maxDt, firstDt, timeStep);

		timeloop.execute(transport);

		printvector(std::cout, variables.saturation(), "saturation", "row", 200, 1);

		return 0;
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
