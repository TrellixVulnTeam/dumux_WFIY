#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> /*@\label{tutorial-decoupled:include-begin}@*/
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/fractionalflow/variableclass.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "tutorial_soilproperties_decoupled.hh"
#include "dumux/material/twophaserelations.hh"
#include "tutorialproblem_decoupled.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include "dumux/timedisc/timeloop.hh" /*@\label{tutorial-decoupled:include-end}@*/


int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2; /*@\label{tutorial-decoupled:dim}@*/

    // create a grid object
    typedef double NumberType; /*@\label{tutorial-decoupled:grid-begin}@*/
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(10); N[0] = 30;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H); /*@\label{tutorial-decoupled:grid-end}@*/

    // define fluid and solid properties and constitutive relationships
    Dune::Water wettingfluid; /*@\label{tutorial-decoupled:water}@*/
    Dune::Oil nonwettingfluid; /*@\label{tutorial-decoupled:oil}@*/
    Dune::TutorialSoil<GridType, NumberType> soil; /*@\label{tutorial-decoupled:soil}@*/
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wettingfluid, nonwettingfluid);/*@\label{tutorial-decoupled:twophaserelations}@*/

    // create object containing the variables
    typedef Dune::VariableClass<GridType, NumberType> VariableType;
    VariableType variables(grid);

    // create object including the problem definition
    typedef Dune::TutorialProblemDecoupled<GridType, NumberType, VariableType> Problem;
    Problem problem(variables, wettingfluid, nonwettingfluid, soil, materialLaw,L, H); /*@\label{tutorial-decoupled:problem}@*/

    // create object including the discretisation of the pressure equation
    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VariableType> DiffusionType;
    DiffusionType diffusion(grid, problem); /*@\label{tutorial-decoupled:diffusion}@*/

    // create object including the space discretisation of the saturation equation
    typedef Dune::FVTransport<GridType, NumberType, VariableType> TransportType;
    TransportType transport(grid, problem); /*@\label{tutorial-decoupled:transport}@*/

    // some parameters used in the IMPES-object
    int iterFlag = 2;
    int nIter = 30;
    double maxDefect = 1e-5;

    // create object including the IMPES (IMplicit Pressure Explicit Saturation) algorithm
    typedef Dune::IMPES<GridType, DiffusionType, TransportType, VariableType> IMPESType;
    IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect); /*@\label{tutorial-decoupled:impes}@*/

    // some parameters needed for the TimeLoop-object
    double tStart = 0; // start simulation at t = tStart
    double tEnd = 1e8; // stop simulation at t = tEnd
    const char* fileName = "tutorial_decoupled"; // name of the output files
    int modulo = 1; // define time step interval in which output files are generated
    double cFLFactor = 1; // security factor for the Courant-Friedrichs-Lewy-Criterion

    // create TimeLoop-object
    Dune::TimeLoop<GridType, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor); /*@\label{tutorial-decoupled:timeloop}@*/

    Dune::Timer timer;
    timer.reset();

    // start simulation
    timeloop.execute(impes); /*@\label{tutorial-decoupled:execute}@*/

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
