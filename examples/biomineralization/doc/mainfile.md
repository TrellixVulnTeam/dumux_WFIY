<!-- Important: This file has been automatically generated by generate_example_docs.py. Do not edit this file directly! -->


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_left: Go back to part 1](modelconcept.md) | [:arrow_right: Continue with part 3](setup.md) |
|---|---|---:|

# Main file

The main file features a time loo with check points to assure the correct timing of the injections matching an experimental setup (see [@Hommel2016], Section 4.2 under the name _column experiment D1_)
The timing of the injections is stored in an additional input file `injections_checkpoints.dat`,
which is read by `main.cc` to set the check points in the time loop.
The number of check points reached is counted and passed on to `problem.hh`, so that `problem.hh` can determine the appropriate injection type.

Additionally, there is an output of the porosity-dependent permeability.

The subsequent file documentation is structured as follows:

[[_TOC_]]

[@Hommel2016]: https://elib.uni-stuttgart.de/handle/11682/8787 "Modelling biogeochemical and mass transport processes in the subsurface: investigation of microbially induced calcite precipitation"


## The main file
This files contains the main program flow for the biomineralization example. Here we can see the programme sequence and how the system is solved using newton's method.


<details open>
<summary><b>Click to hide/show the file documentation</b> (or inspect the [source code](../main.cc))</summary>


### Included header files
<details>
These are DUNE helper classes related to parallel computations

```cpp
#include <dune/common/parallel/mpihelper.hh>
```

The following headers include functionality related to property definition or retrieval, as well as
the retrieval of input parameters specified in the input file or via the command line.

```cpp
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>
```

The following files contain the nonlinear Newtown method, the linear solver and the assembler

```cpp
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
```

The following class provides a convenient way of writing of dumux simulation results to VTK format.

```cpp
#include <dumux/io/vtkoutputmodule.hh>
```

The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
Many different Dune grid implementations are supported, of which a list can be found
in `gridmanager.hh`.

```cpp
#include <dumux/io/grid/gridmanager_yasp.hh>
```

We include the problem file which defines initial and boundary conditions to describe our example problem
remove #include "problem.hh"

```cpp
#include "properties.hh"
```

</details>

### The main function
We will now discuss the main program flow implemented within the `main` function.
At the beginning of each program using Dune, an instance `Dune::MPIHelper` has to
be created. Moreover, we parse the run-time arguments from the command line and the
input file:

```cpp
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);
```

For convenience we define the type tag for this problem.
The type tags contain all the properties that are needed to run the simulations.

```cpp
    using TypeTag = Properties::TTag::MICPColumnSimpleChemistry;
```

### Step 1: Create the grid
The `GridManager` class creates the grid from information given in the input file.
This can either be a grid file, or in the case of structured grids, one can specify the coordinates
of the corners of the grid and the number of cells to be used to discretize each spatial direction. The latter is the case for this example.

```cpp
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
```

### Step 2: Setting up the problem
We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.

```cpp
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
```

We now instantiate the problem, in which we define the boundary and initial conditions.

```cpp
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
```

get some time loop parameters

```cpp
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
```

We initialize the solution vector

```cpp
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;
```

on the basis of this solution, we initialize the grid variables

```cpp
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
```

We initialize the vtk output module. Each model has a predefined model specific output with relevant parameters for that model.

```cpp
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    GetPropType<TypeTag, Properties::IOFields>::initOutputModule(vtkWriter); //!< Add model specific output fields
    //we add permeability as a specific output
    vtkWriter.addField(problem->getPermeability(), "Permeability");
    //We update the output fields before write
    problem->updateVtkOutput(x);
    vtkWriter.write(0.0);
```

We instantiate the time loop and define an episode index to keep track of the the actual episode number

```cpp
    int episodeIdx = 0;
    //We  set the initial episodeIdx in the problem to zero
    problem->setEpisodeIdx(episodeIdx);
    //We  set the time loop with episodes (check points)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
```

### Step 3 Setting episodes specified from file
In this example we want to specify the injected reactant solutions at certain times.
This is specified in an external file injections_checkpoints.dat.
Within the following code block, the parameter file read and the time for the injection are extracted.
Based on that, the checkpoints for the simulation are set.

```cpp
    //We use a time loop with episodes if an injection parameter file is specified in the input
    if (hasParam("Injection.NumInjections"))
    {
        //We first read the times of the checkpoints from the injection file
        const auto injectionCheckPoints = readFileToContainer<std::vector<double>>("injection_checkpoints.dat");

        // We set the time loop with episodes of various lengths as specified in the injection file
        // set the episode ends /check points:
        timeLoop->setCheckPoint(injectionCheckPoints);

        // We also set the initial episodeIdx in the problem to zero, as in the problem the boundary conditions depend on the episode.
        problem->setEpisodeIdx(episodeIdx);
    }

    // If nothing specified, we do not need to use episodes
    else
    {
        // In this case, we set the time loop with one big episodes
        timeLoop->setCheckPoint(tEnd);
    }
```

### Step 4 solving the instationary problem
We create and initialize the assembler with a time loop for the transient problem.
Within the time loop, we will use this assembler in each time step to assemble the linear system.
Additionally the linear and non-linear solvers are set

```cpp
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    //We set the linear solver
    using LinearSolver = Dumux::ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    //We set the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
```

#### The time loop
We start the time loop and solve a new time step as long as `tEnd` is not reached. In every time step,
the problem is assembled and solved, the solution is updated, and when a checkpoint is reached the solution
is written to a new vtk file. In addition, statistics related to CPU time, the current simulation time
and the time step sizes used is printed to the terminal.

```cpp
    timeLoop->start(); do
    {
        // We set the time and time step size to be used in the problem
        problem->setTime( timeLoop->time() + timeLoop->timeStepSize() );
        problem->setTimeStepSize( timeLoop->timeStepSize() );

        // Solve the linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // Update the old solution with the one computed in this time step and move to the next one
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        // We update the output fields before write
        problem->updateVtkOutput(x);

        // We write vtk output on checkpoints or every 20th timestep
        if (timeLoop->timeStepIndex() % 20 == 0 || timeLoop->isCheckPoint())
        {
            // update the output fields before write
            problem->updateVtkOutput(x);

            // write vtk output
            vtkWriter.write(timeLoop->time());
        }
        // We report statistics of this time step
        timeLoop->reportTimeStep();

        //If episodes/check points are used, we count episodes and update the episode indices in the main file and the problem.
        if (hasParam("Injection.NumInjections") || hasParam("TimeLoop.EpisodeLength"))
        {
            if (timeLoop->isCheckPoint())
            {
                episodeIdx++;
                problem->setEpisodeIdx(episodeIdx);
            }
        }

        // set the new time step size that is suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());
    timeLoop->finalize(leafGridView.comm());
```


```cpp

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
}
```

### Exception handling
In this part of the main file we catch and print possible exceptions that could
occur during the simulation.
<details><summary> Click to show error handler</summary>

```cpp
// errors related to run-time parameters
catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
```

</details>

</details>


| [:arrow_left: Back to the main documentation](../README.md) | [:arrow_left: Go back to part 1](modelconcept.md) | [:arrow_right: Continue with part 3](setup.md) |
|---|---|---:|

