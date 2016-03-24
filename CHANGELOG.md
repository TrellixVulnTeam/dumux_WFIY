Differences Between DuMuX 2.8 and DuMuX 2.9
===================================================

* IMPORTANT NOTES:
    - DuMuX 2.9 is expected to run based on either Dune 2.4 or Dune 3.0. We will
      try to keep the compatibility with Dune 3.0 as long as it is technically
      feasible and our resources allow it. If you want to use Dumux multidomain
      models, you have to stick with the Dune 2.4 core and specific versions of
      other modules, see `test/multidomain/README` for details.

    - The ALUGrid stand-alone library cannot be used any longer. Use the module
      [dune-alugrid](https://gitlab.dune-project.org/extensions/dune-alugrid)
      instead. Both the `releases/2.4` branch and the `master` should work.

    - The above is not true if you like to run a sequential MPFA model on an
      `ALUGrid`. Then, you currently have to use the `master-old` branch of
      dune-alugrid. We will try to fix this as soon as possible. Alternatively,
      `UGGrid` or `YaspGrid` can be chosen as grid managers.

    - Instead of using AlbertaGrid for the tests where dim < dimWorld, we now
      employ
      [dune-foamgrid](https://gitlab.dune-project.org/extensions/dune-foamgrid).
      Dune-foamgrid provides 1d and 2d simplex grids embedded in an arbitrary
      dimension world space. It features element parametrizations, runtime growth,
      runtime-movable vertices. You might still use AlbertaGrid, but it is not
      supported by our GridCreator.

    - If you like/have to use corner-point grids by means of the module
      dune-cornerpoint, you have to use (and partially patch) the 2015.10 release
      of [OPM](http://opm-project.org/?page_id=36). See `patches/README` for
      details.

* IMPROVEMENTS and ENHANCEMENTS:
    - The folder structure has been changed according to
      [FS#250](http://www.dumux.org/flyspray/index.php?do=details&task_id=250).
      This has been a rather massive change affecting more than 1000 files. Close
      to 400 files have been moved and/or renamed.
      We made everything backwards-compatible, the worst thing that should happen
      after switching to Dumux 2.9, will be some warnings when including headers
      from old destinations/names. You can fix the include statements and get rid
      of the warnings by applying the bash script `bin/fix_includes.sh` to your
      source files, for example by executing
      ```
      bash ../dumux/bin/fix_includes.sh file1 [file2 ...]
      ```
      or
      ```
      find . -name '*.[ch][ch]' -exec bash ../dumux/bin/fix_includes.sh {} \;
      ```
      inside the folder that contains your files.

      A patch is available to remove deprecated header files:
      ```
      patch -p1 < patches/dumux-2.9-no-deprecated-headers.patch
      ```

      The benefits are hopefully:
        + A clearer structure in terms of the problems that you want to apply Dumux
          for. Three main application areas on the top level: `porousmediumflow`,
          `freeflow` and `geomechanics`. The different numerical treatments "fully
          implicit" or "sequential" appear as discretization detail after the
          choice of the physical model. That's of course currently rather wishful
          thinking, but nevertheless where we are headed. The folder `implicit` on
          the top level now only contains physics-agnostic classes that can be used
          by any class of an application area.  Please note the change from
          "decoupled" to "sequential" according to the related task
          [FS#252](http://www.dumux.org/flyspray/index.php?do=details&task_id=252).

        + Nicer include statements due to relaxation of the naming conventions for
          the header files. Compare the old
          ```
          #include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnilocaloperator.hh>
          ```
          with the new
          ```
          #include <dumux/multidomain/2cnistokes2p2cni/localoperator.hh> 
          ```
      The structure change is reflected in the `test` folder:
        + The tests from`test/implicit/particular_model` have been moved to
          `test/porousmediumflow/particular_model/implicit`. For example,
          `test/implicit/2p` has been moved to `test/porousmediumflow/2p/implicit`.

        + Analogously, the tests from `test/decoupled/particular_model` have been
          moved to `test/porousmediumflow/particular_model/sequential`.

        + The subfolders `decoupled` and `implicit` of `test` have been removed.

        + If you have cloned the Dumux Git repository and have local changes in the
          folders `test/implicit` or `test/decoupled`, you can expect merge
          conflicts for your next `git pull`. You can either deal with these
          conflicts directly or create a patch, remove the local changes, pull, and
          apply the patch afterwards with some care to respect the changed
          structure.

    - A two-phase multiple-interacting-continua (MINC) model has been added to
      the Dumux model portfolio. See `test/porousmediumflow/2pminc/implicit` for
      details.

    - The multidomain models have been restructured. Duplicated code has been
      reduced; isothermal and non-isothermal models are treated in a more
      consistent manner.

    - It is now possible to specify point sources for implicit models. A point
      source is a source term specified at any point location in e.g. kg/s. Dumux
      will compute the correct control volume the source belongs to for you. Point
      sources can be e.g. solution and/or time-dependent. See tests
      (1p/implicit/pointsources, 2p/implicit/pointsources) for examples.

    - All tests use our standard `GridCreator` now. If it is possible to specify
      the grid entirely in the input-file, the corresponding DGF files have been
      deleted. In particular, a YaspGrid tensor grid can now also be specified via
      the input file only.

    - Several sections on our fluid/material framework have been moved from the
      handbook to the Doxygen documentation.

    - The three-phase constitutive relations from
      `material/fluidmatrixinteractions` have been reworked to be consistent with
      their two-phase analogues. In particular, an `EffToAbsLaw` and
      regularization classes have been implemented.

    - In case of a simulation stop due to too many timestep subdivisions, restart
      files of both the current and the old solution are automatically generated.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
    - All flux variables are now default-constructible. While the non-trivial
      constructors are deprecated, model implementers might be required to make
      their flux variables default-constructible too. In particular, this affects
      you if you develop your own flux variables that
        + inherit from flux variables from dumux-stable, such as the 
          `ImplicitDaryFluxVariables`,
        + and/or are used in a local residual from dumux-stable.
      See the
      [mailing list](https://listserv.uni-stuttgart.de/pipermail/dumux/2016q1/001551.html)
      for details.

    - For the multidomain models, the notation of the boundary condition types
      has changed. This is especially important for all momentum boundary
      conditions. In general:
        + `couplingInflow`  -> `couplingNeumann`
        + `couplingOutflow` -> `couplingDirichlet`
    - But for the momentum balances:
        + `couplingInflow`  -> `couplingDirichlet`
        + `couplingOutflow` -> `couplingNeumann`

    - Due to the change in the three-phase fluid-matrix-interactions, you might
      have to adjust your spatial parameters. You should get a compiler warning
      message that gives you more details.

    - The TypeTags `ImplicitModel` and `ExplicitModel` have been deleted. They
      haven't been used apart from one internal inheritance. See FS#304 for
      details.

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.9: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.
    - The word `Decoupled` in the TypeTags has been replaced by `Sequential`:
        + `DecoupledModel` -> `SequentialModel`
        + `DecoupledOneP` -> `SequentialOneP`
        + `DecoupledTwoP` -> `SequentialTwoP`
        + `DecoupledTwoPTwoC` -> `SequentialTwoPTwoC`
        + `DecoupledTwoPTwoCAdaptive` -> `SequentialTwoPTwoCAdaptive`


* Deprecated CLASSES/FILES, to be removed after 2.9:
    - Self-written parallel linear solvers and corresponding infrastructure,
      according to
      [FS#293](http://www.dumux.org/flyspray/index.php?do=details&task_id=293).
      For parallel runs, use the `AMGBackend` instead. For sequential runs,
      direct replacements are:
        + `BoxBiCGStabILU0Solver` -> `ILU0BiCGSTABBackend`
        + `BoxBiCGStabSORSolver` -> `SORBiCGSTABBackend`
        + `BoxBiCGStabSSORSolver` -> `SSORBiCGSTABBackend`
        + `BoxBiCGStabJacSolver` -> `JacBiCGSTABBackend`
        + `BoxBiCGStabGSSolver` -> `GSBiCGSTABBackend`
        + `BoxCGILU0Solver` -> `ILUnCGBackend`
        + `BoxCGSORSolver` -> `SORCGBackend`
        + `BoxCGSSORSolver` -> `SSORCGBackend`
        + `BoxCGJacSolver` -> `JacCGBackend`
        + `BoxCGGSSolver` -> `GSCGBackend`
        + `IMPETBiCGStabILU0Solver` -> `ILU0BiCGSTABBackend`

    - `CubeGridCreator`, functionality available in default `GridCreator`

    - `SimplexGridCreator`, functionality available in default `GridCreator`

    - `DgfGridCreator`, functionality available in default `GridCreator` (since 2.8)

    - `Decoupled...Indices` -> `Sequential...Indices` (BEWARE: maybe no compiler
    warnings)

* Deprecated MEMBER FUNCTIONS, to be removed after 2.9:

* Deprecated protected MEMBER VARIABLES, to be removed after 2.9: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable
  is used.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMuX 2.8:
    - Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.7 and DuMuX 2.8
===================================================

* IMPORTANT NOTES:
  - DuMuX 2.8 is expected to run based on either Dune 2.3 or Dune 2.4. However,
    no attempt has been made to fix the warnings arising from the deprecation of
    EntityPointer in Dune 2.4. This will be made during the next release cycle.
    Moreover, using the multidomain models based on Dune 2.4 is currently only
    possible by patching dune-multidomaingrid. See test/multidomain/README for
    details.

* DELETED BUILD SYSTEM: The Autotools based build system was removed, use the
  CMake based build system as it is default since Dune 2.4.

* IMPROVEMENTS and ENHANCEMENTS:
  - New fully-implicit porous-media models for two fluid phases that may consist
    of an arbitrary number of components have been added. The basic one is
    associated with the property TwoPNC, see test/implicit/2pnc. A more advanced
    one that incorporates solid-fluid phase changes is indicated by TwoPNCMin,
    see test/implicit/2pncmin.

  - The implicit cell-centered models now can use adaptive grid refinement. To
    make a test problem adaptive, just set the property AdaptiveGrid to true and
    choose corresponding indicators via AdaptionInitializationIndicator and
    AdaptionIndicator, see test/implicit/2p/lensproblem.hh for an example. So
    far, indicators are only provided for the TwoPModel. Indicators for other
    models will be provided in the future, as well as parallelization and box
    discretization.

  - With the CpGridCreator, a grid creator has been introduced that reads from a
    Petrel output / Eclipse input file and generates a CpGrid that is offered by
    the OPM module dune-cornerpoint. The fully-implicit cell-centered models are
    now able to deal with cornerpoint grids. See
    test/implicit/2p/test_cc2pcornerpoint for a test of the functionality. A
    realistic corner-point grid will be provided in dumux-lecture soon. The OPM
    modules need to be patched to be compatible with Dune's CMake based build
    system, see patches/README for details.

  - Zero equation turbulence models (zeroeq) have been added as new models
    to the freeflow folder. Tests for coupling a turbulent free flow using
    zeroeq turbulence models with flow in a porous medium have been added.

  - A new class GridCreator is now the new standard grid creator replacing
    DgfGridCreator. It comprises all functionality from the DgfGridCreator (see
    also immediate interface changes), plus the ability to read gmsh, or to
    build a structured grid (only with Dune 2.4) by merely changing the input
    file.

  - Multidomain problems can now be run by using the general start routine,
    just as most other problems. For this, the constructor of the multidomain
    problems has been changed and the InterfaceMeshCreator has been replaced by
    the InterfaceGridCreator, see below.

  - The Richards model has now an additional flag useHead, which can be used to
    switch between a pressure-saturation and a pressureHead-watercontent
    formulation. The primary variables are either pressure in [Pa] or pressure
    head in [cm], respectively. Default is useHead = false. See
    test/implicit/richards for details.

  - A bug in the diffusion term in the freeflow models has been fixed.

  - A lot of work has been devoted to improving the testing environment, adding
    new tests, restructuring the handbook and improving the documentation.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:

  - The (new standard) GridCreator's public method gridPtr has been removed. For
    using PARAMETERS from DGF files use the GridCreators parameters method via
    GridCreator::parameters(...). See test/implicit/co2 for an example.

  - The use and support for SGrid is dropped. SGrid is deprecated in Dune 2.4.
    Use YaspGrid instead.

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.8: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.

* Deprecated CLASSES/FILES, to be removed after 2.8:
  - The InterfaceMeshCreator has been moved to InterfaceGridCreator and adapted
    to the structure of other grid creators, it can simply be used by specifying
    the GridCreator TypeTag.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.8:
  - The constructor of the multidomain problems has changed. They will now be
    called directly from the start.hh, identical to the other problems. The new
    parameters are the TimeManager and the HostGrid.

  - The method simulate(Scalar dtInitial, Scalar tEnd) from MultiDomainProblem,
    is unused and will be dropped.

  - The GnuplotInterface functions are now called without giving a window
    number. If plots should be plotted in different windows, different
    GnuplotInterface objects are now required. This affects also all other plots
    in the "io" folder.

  - The write() function in plotoverline2d.hh now has an append function, to be
    able to decide whether the previously written file should be kept or not.

* Deprecated protected MEMBER VARIABLES, to be removed after 2.8: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable
  is used.

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMuX 2.7:
  Everything listed as deprecated below has been removed.

Differences Between DuMuX 2.6 and DuMuX 2.7
===================================================

* IMPORTANT NOTES:
  - DuMuX 2.7 should work with DUNE 2.3 as well as 2.4. However, at the time of
    writing, DUNE-multidomain(grid) doesn't work with DUNE 2.4. Therefore, if a
    DuMuX multidomain model should be used, DUNE 2.3 is required. See
    test/multidomain/README for details.

  - The 2.3 branch of dune-alugrid has no CMake support, use dune-alugrid master
    respectivly 2.4. Or you can fall back to Autotools or use legacy ALUGrid
    1.52.

* IMPROVEMENTS and ENHANCEMENTS:
  - Since 2.6, all isothermal implicit porous-media models (except 2pdfm) can be
    easily enhanced with the energy equation. For 2.7, this has been also
    carried out for the models that were only isothermal before, namely, 1p, 3p
    and richards. Tests have been written and are provided in test/implicit. In
    order to keep the number of subfolders bearable, isothermal as well as
    thermal tests are gathered in the model folders "1p", "1p2c", ..., "3p3c",
    "mpnc", "richards" (without the "ni") and the corresponding "ni"-folders
    have been deleted.

  - All implicit porous-media models (except 2pdfm) are now able to run on grids
    with dim < dimWorld. In implicit/1p, four new tests are added that run the
    1p test problem on 1d-3d and 2d-3d Alberta grids with box and cell-centered,
    respectively. Compilation has been tested also for all other models, but no
    runtime testing has been performed.

  - The terminology for the Newton method has been improved according to FS#238.
    In particular, what has been referred to as "relative error" is now termed
    "maximum relative shift", while "absolute error" has been renamed to
    "residual reduction". This is particularly important, if corresponding
    parameters or properties are set, see below.

  - The geomechanics ElTwoPModel runs in parallel now. This is made possible by
    a dedicated solver, the El2PAMGBackend which has to be set for the property
    LinearSolver in the problem file. See test/geomechanics/el2p for details.

  - Before, velocity output for the implicit porous-media models only worked for
    cube grids. This has been generalized to simplices (box and cc) and prisms/
    pyramids (box only).

  - Revised and fixed restart capability for the multidomain models.

  - A gnuplot interface has been added (works only with CMake). With
    this interface it is possible to plot material laws (like in the 2p2cni test),
    or to generate live-updating output (like in test_2cnistokes2p2cni).
    The gnuplot interface reads analytical functions, data file or data arrays.

  - The fuzzycompare script for automatic testing has been improved. Instead of
    printing only the first deviation from the reference solution, it now
    prints the maximum deviation in each field/variable.

* DEPRECATED BUILD SYSTEM: DuMuX 2.7 will be the last release which can be built
    with the Autotools based build system. It is deprecated and will be removed
    for DuMuX 2.8. We encourage the change towards CMake, especially with the
    upcoming DUNE 2.4.
    The warning can be suppressed with --disable-dumux-deprecated-autotools

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - Before, the "heatCapacity" function in the spatial parameters and volume
    variables of the implicit nonisothermal models was a misnomer, since it
    returned an effective quantity, namely,
    heatCapacity*density*(1 - porosity) in [J/(K m^3)].
    Except for mpnc, which resulted in an additional inconsistency.
    Corresponding to the decision documented in FS#216, the function has been
    renamed to "solidHeatCapacity" and returns always the "true" (non-effective)
    heat capacity in [J/(kg K)]. This requires an additional function
    "solidDensity" which returns the mass density of the porous matrix.
    Moreover, the functions "thermalConductivitySolid/Fluid" are renamed to
    "solid/fluidThermalConductivity". The decision to prepend with "solid/fluid"
    rather than to append is motivated by consistency with components and fluid
    systems, where "gas" and "liquid" are always prepended to the corresponding
    function names.
    Therefore, it might be necessary to adapt your thermal solid parameters in
    the spatialparams file such that they offer functions "solidHeatCapacity",
    "solidDensity" and "solidThermalConductivity". See
    test/implicit/2p2c/injectionspatialparams.hh for an example.

  - Due to the change in the Newton terminology (see above), there exist two
    backward-compatibility breakages:
    . If a model re-implements the function "relativeErrorDof", it has to be
      renamed to "relativeShiftAtDof". See dumux/implicit/implicitmodel.hh for
      an example.

    . If a NewtonController re-implements the function "newtonUpdateRelError",
      it has to be renamed to "newtonUpdateShift". See
      dumux/nonlinear/newtoncontroller.hh for an example.

  - The properties "AMGPDELabBackend" and "AMGLocalFemMap" have been unified to
    "AmgTraits".

* Deprecated PROPERTY and PARAMETER NAMES, to be removed after 2.7: BEWARE: The
  compiler will not print any warning if a deprecated property or parameter name
  is used. However, a run-time warning should appear in the summary lines after
  the corresponding run.
  - Corresponding to the improved Newton terminology, the following properties
    (prepended with "Newton") and parameters (in the group "Newton") are
    renamed:
    AbsTolerance -> ResidualReduction
    EnableAbsoluteCriterion -> EnableResidualCriterion
    RelTolerance -> MaxRelativeShift
    EnableRelativeCriterion -> EnableShiftCriterion
    SatisfyAbsAndRel -> SatisfyResidualAndShiftCriterion

* Deprecated CLASSES/FILES, to be removed after 2.7:
  - SeqAMGBackend and ScaledSeqAMGBackend, replaced by AMGBackend.

  - P0LocalFiniteElementMap.

  - CellData2P2Cmultiphysics, replaced by CellData2P2CMultiPhysics.

  - BoxLocalOperator from dumux/multidomain/common/pdelablocaloperator.hh.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.7:
  - The functions "heatCapacity", "densitySolid" (mpnc only) and
    "thermalConductivitySolid/Fluid" in the VolumeVariables of the nonisothermal
    implicit porous-media models: use "solidHeatCapacity", "solidDensity" and
    "solid/fluidThermalConductivity" instead. See also the immediate interface
    changes above.

  - In dumux/implicit/common/implicitmodel.hh and
    dumux/geomechanics/el2p/elp2basemodel.hh:
    "relativeErrorDof" -> "relativeShiftAtDof"

  - In dumux/nonlinear/newtoncontroller.hh:
    "setRelTolerance" -> "setMaxRelativeShift"
    "setAbsTolerance" -> "setResidualReduction"
    "newtonUpdateRelError" -> "newtonUpdateShift"

  - The 1p2c volume variables no longer use the method tortuosity() from
    spatial params class, the value is now calculated within the effective
    diffusivity model. Thus the method is deprecated in the spacial params
    classes FVSpatialParamsOneP and ImplicitSpatialParamsOneP.

* Deprecated protected MEMBER VARIABLES, to be removed after 2.7: BEWARE: Older
  compilers will not print any warning if a deprecated protected member variable 
  is used.
  - In dumux/nonlinear/newtoncontroller.hh:
    "error_" -> "shift_"
    "lastError_" -> "lastShift_"
    "tolerance_" -> "shiftTolerance_"
    "absoluteError_" -> "reduction_"
    "lastAbsoluteError_" -> "lastReduction_"
    "initialAbsoluteError_" -> "initialResidual_"
    "absoluteTolerance_" -> "reductionTolerance_"
    "enableRelativeCriterion_" -> "enableShiftCriterion_"
    "enableAbsoluteCriterion_" -> "enableResidualCriterion_"
    "satisfyAbsAndRel_" -> "satisfyResidualAndShiftCriterion_"

* DELETED classes/files, property names, constants/enums,
  member functions, which have been deprecated in DuMuX 2.6:
  Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.5 and DuMuX 2.6
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - For the non-isothermal porous media models, the energy equation is
    implemented in a more generic way for all models in
    dumux/implicit/nonisothermal. The existing TypeTag names like TwoPNI stay
    the same. If a new non-isothermal model should be used, it is important
    to NOT include anything from the old model-specific implementation like
    from dumux/implicit/2pni, but to include from the model folder without the
    "ni". See test/implicit/2pni for details. In principle, any isothermal
    porous media model can be enhanced with the energy equation. Ideally, only
    the corresponding property files have to be augmented. See
    dumux/implicit/2p/2ppropert*.hh for details. The 1p2c model already has
    been enhanced, the remaining models will follow in 2.7.

  - The AMG backend is based directly on dune-istl now. No PDELab is required
    anymore. The tests so far exhibit an improved robustness. Thanks to Markus
    Blatt for the work.

  - The multidomain models can now be used with the 2.3 release versions of the
    DUNE core modules and dune-multidomaingrid, and the 2.0 release versions
    of dune-pdelab and dune-multidomain. See test/multidomain/README for
    details.

  - In the fully implicit mpnc model, a further specialization allows now to
    describe two-phase flow with two energy equations.

  - The free flow models now include the component enthalpy fluxes transported
    by diffusion processes (h^k D grad x), which was not considered before.

  - UMFPack is a new direct linear solver and can be use as a drop-in
    replacement for SuperLU. Some users claim a speed-up up to a factor of
    seven. We know cases where it was 10% slower, so please measure for your
    problems.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - The 3p3cni model now also uses an effective thermal conductivity model
    (ETCM). The ETCM is easily exchangeable. The default one is
    ThermalConductivitySomerton, which is implemented in
    dumux/material/fluidmatrixinteractions/3p. The ETCM requires that 3p3cni
    spatial parameters provide a function thermalConductivitySolid instead of
    matrixHeatFlux. See test/implicit/3p3cni/columnxylolspatialparams.hh for
    details. Moreover, the employed fluid system has to actually implement the
    function thermalConductivity. See
    dumux/material/fluidsystems/h2oairxylenefluidsystem.hh for details.

  - The non-isothermal flux variables call the effective thermal conductivity
    models (ETCM) in a different way. If you used a self-written ETCM and want
    to use a new non-isothermal model, the ETCM has to be adapted. See
    material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh for
    details.

  - Fully implicit mpnc model: in order to account for the possibility of using
    two energy equations, the boolean property EnableKineticEnergy has been
    changed to the integer property NumEnergyEquations.

* Deprecated way of setting command line parameters, to be removed after 2.6:
  - To set paramaters from the command line, the notation --parameterFile=NAME
    is deprecated. Use from now on -ParameterFile NAME.

* Deprecated CLASSES/FILES, to be removed after 2.6:
  - The old non-isothermal porous media models are deprecated. Technically,
    including a ..niproperties.hh file triggers a deprecation warning.

  - FVPressure2P2CAdaptive, use dimension-specific implementations
    FV2dPressure2P2CAdaptive and FV3dPressure2P2CAdaptive instead.

  - FVTransport2P2CAdaptive, use dimension-specific implementations
    FV2dTransport2P2CAdaptive and FV3dTransport2P2CAdaptive instead.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.6:
  - In the Stokes flux variables, the method eddyViscosity() is deprecated, use
    dynamicEddyViscosity() instead.

  - In the Stokes non-isothermal flux variables, the method eddyConductivity()
    is deprecated, use thermalEddyConductivity() instead.

  - Already in 2.5, the following member functions of MultiDomainModel/Problem
    have been deprecated: subProblemX, subModelX, subIDX, gridViewX with X=1,2.
    They are replaced by sdProblemX, sdModelX, sdIDX, sdGridViewX.

* DELETED classes/files, property names, constants/enums, 
  member functions, which have been deprecated in DuMuX 2.5:
  Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.4 and DuMuX 2.5
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - The three-dimensional implementation of the MPFA L-method is made
    available for the decoupled compositional 2p2c models. It also allows
    for simulation on an adaptive grid.
  - Coupling of 2c2p with stokesnc and 2p2cNI with Stokesncni was
    added. The stokes2c and stokes2cni are now DEPRECATED and will be kicked 
    out by the next release. Instead generalized stokesnc and stokesncni 
    models are introduced. Unlike 2c models the transport equations in 
    the nc models are capapable of using both mass and mole fractions.
    NOTE: For coupling test examples be aware of the harsh version 
    restrictions mentioned in dumux/test/modelcoupling/README.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - Dropped support for PDELab 1.0.1, PDELab 1.1 is required.

* Deprecated CLASSES/FILES, to be removed after 2.5:
  - Stokes2cModel was replaced by StokesNCModel, similar for more
    Stokes2c* classes.

* DELETED classes/files, property names, constants/enums, 
  member functions, which have been deprecated in DuMuX 2.4:
  Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.3 and DuMuX 2.4
===================================================

* IMPORTANT NOTES:
  - If the current trunk version of DUNE (2.3) is used, the co2 and co2ni tests
    require that the DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS flag is set. This
    flag is needed in order to be able to use the function boundaryId(). An
    error is thrown during compilation if the flag is not set. For reference see
    the commented part in optim.opts and debug.opts.

  - All two-component models (1p2c, 2p2c, 2p2cni, co2, co2ni) can now be used 
    with either mole or mass fractions. The property useMoles has to be set in 
    the problem file and the boundary conditions have to be choosen accordingly.
    . 1p2c, 2p2c, 2p2cni use MOLE fractions by default.
    . co2, co2ni use MASS fractions by default.
    . For completeness: 3p3c, 3p3cni and mpnc use MOLE fractions only.

* IMPROVEMENTS and ENHANCEMENTS:
  - Three geomechanics models have been added, thanks to the previous PhD work
    of Melanie Darcis. The models are located in dumux/geomechanics, a test for
    each model is provided in test/geomechanics:
    . elastic: linear elasticity equations that account only for the solid
               mechanics.
    . el1p2c:  poroelasticity equations for one solid and one fluid phase, where
               the fluid phase is admitted to consist of two components.
    . el2p:    poroelasticity equations for one solid phase and two immiscible
               fluid phases.

  - A three-dimensional implementation of the MPFA L-method is added and made
    available for the decoupled 2p models. It also allows for simulation on
    an adaptive grid.

  - The MPNC model now allows to drop the assumptions of local thermal and/or
    chemical equilibrium. Dropping the local chemical equilibrium assumption
    leads to mole fractions in different phases that are not determined via
    equilibrium relations. If the local thermal equilibrium assumption is
    dropped, phases (fluids and solid) are locally allowed to have different
    temperatures. See test/implicit/mpnc/test_boxmpnckinetic for an example.

  - A fully-implicit three-phase immiscible model has been added. The headers
    are contained in dumux/implicit/3p, tests are provided in test/implicit/3p.

  - The handling of Dirichlet boundary conditions for the fully-implicit cell-
    centered models has been improved. Now, the conditions are evaluated at the
    centers of the corresponding boundary faces. In general, they are
    incorporated into the local residual in a weak sense. Only for mixed
    Dirichlet/Neumann conditions and for the MPNC model, they still are
    incorporated in a strong sense.

  - The sequential models can use a sub-timestepping in the transport scheme, if
    the property "SubCFLFactor" is chosen smaller than "CFLFactor", which in
    that case specifies the CFL factor used in the pressure equation.

  - All fully-implicit porous-media models now provide the possibility to write
    the velocities to the VTK output. This can be achieved by setting the
    parameter "AddVelocity" of the group "Vtk" or the corresponding property 
    "VtkAddVelocity." See test/implicit/1p2c for an example.

  - The CMake build system support uses the experimental mechanisms from DUNE.
    Recent versions of Dune 2.3-svn or newer are required to use CMake.

  - Our naming rules have been refined and enforced for member functions and
    variables. Many inconsistencies could be removed, especially the special
    status of the capitalized "S" indicating saturation. See the deprecation
    listings below or FS#180 for details.

  - Misleading names in the fully-implicit models that still contained "b/Box"
    have been replaced by more generic ones. See the deprecation listings below
    or FS#194 for details.

  - The ...FVElementGeometry classes have been cleaned up a bit. See the 
    deprecation listings below or FS#192 for details.

  - Added compiler support for Clang 3.2, Clang 3.3, and GCC 4.8.

* Deprecated CLASSES/FILES, to be removed after 2.4:
  - OnePBoxModel, OnePTwoCBoxModel -> OnePModel, OnePTwoCModel

  - Headers containing capitalized letters are replaced by their non-capitalized
    analogs. 
    In dumux/decoupled: 1p/cellData1p.hh, 1p/fluxData1p.hh, 
    2p/cellData2padaptive.hh, 2p/fluxData2p.hh, 2p/cellData2p.hh,
    2p2c/cellData2p2c.hh, 2p2c/cellData2p2cadaptive.hh, 2p2c/fluxData2p2c.hh,
    2p2c/cellData2p2cmultiphysics.hh.
    In dumux/material/fluidmatrixinteractions/3p: parkerVanGen3p.hh,
    parkerVanGen3pparams.hh.

* Deprecated CONSTANTS/ENUMS, to be removed after 2.4: BEWARE: Older compilers 
  will not print any warning if a deprecated constant/enum is used.
  - saturation indices: S0Idx, SnIdx, SnOrXIdx, SOrX1Idx, SOrX2Idx, SwIdx, 
    SwOrXIdx
    -> s0Idx, snIdx, snOrXIdx, sOrX1Idx, sOrX2Idx, swIdx, swOrXIdx

  - two-phase formulations: pnSn, pnSw, pwSn, pwSw -> pnsn, pnsw, pwsn, pwsw

  - DecoupledTwoPCommonIndices: pressureNW, saturationNW, velocityNW 
    -> pressureNw, saturationNw, velocityNw

  - DecoupledTwoPIndices: pressEqIdx -> pressureEqIdx

  - MPNCIndices: NumPrimaryEnergyVars, NumPrimaryVars
    -> numPrimaryEnergyVars, numPrimaryVars

* Deprecated public MEMBER VARIABLES, to be removed after 2.4: BEWARE: Older
  compilers will not print any warning if a deprecated public member variable 
  is used.
  - ...FVElementGeometry: numEdges, numFaces, numFap, numVertices
    -> numScvf, -, numFap of each subcontrolvolume face, numScv

  - BoxFVElementGeometry: edgeCoord, faceCoord

* Deprecated MACROS, to be removed after 2.4: BEWARE: The compiler will not
  print any warning if a deprecated macro is used.
  - DUMUX_ALWAYS_INLINE whether the according attribute is supported

* Deprecated MEMBER FUNCTIONS, to be removed after 2.4:
  - all problems: bboxMin/Max() -> bBoxMin/Max()

  - ImplicitProblem: boxSDNeumann(), boxSDSource()
    -> solDependentNeumann(), solDependentSource()

  - ImplicitPorousMediaProblem: boxGravity(), boxTemperature() 
    -> gravityAtPos(), temperatureAtPos() (different signatures!)

  - fluid-matrix-interactions: dkrn_dSw(), dkrw_dSw(), dpc_dSw(), pC(), 
    dSw_dpC(), Sgr(), Snr(), SnToSne(), Sw(), Swr(), SwToSwe()
    -> dkrn_dsw(), dkrw_dsw(), dpc_dsw(), pc(), dsw_dpc(), sgr(), snr(), 
    snToSne(), sw(), swr(), swToSwe()

  - LinearMaterial(Params): entryPC(), maxPC(), setEntryPC(), setMaxPC() 
    -> entryPc(), maxPc(), setEntryPc(), setMaxPc()

  - RegularizedVanGenuchten(Params): pCHighSw(), pCLowSw() 
    -> pcHighSw(), pcLowSw()

  - VanGenuchtenParams, ParkerVanGen3PParams: setVgM(), setVgN(), vgM(), vgN()
    -> setVgm(), setVgn(), vgm(), vgn()

  - ParkerVanGen3P(Params): betaGN(), betaGW(), betaNW(), pCAlpha(), pCGN(), 
    pCGW(), pCNW(), setBeta..., setkrRegardsSnr(), Swrx()
    -> betaGn(), betaGw(), betaNw(), pcAlpha(), pcgn(), pcgw(), pcnw(), 
    setBeta..., setKrRegardsSnr(), swrx()

  - MPLinearMaterialParams: Sreg() -> sReg()

  - EvalCflFlux...: getCFLFluxFunction() -> getCflFluxFunction

  - FVMPFAOInteractionVolume: getNTKNu_by_dF(), getNTKNu(), getNTKrKNu_by_dF(), 
    getNTKrKNu()
    -> getNtkNu_df(), getNtkNu(), getNtkrkNu_df(), getNtkrkNu()

  - TwoPDFMVolumeVariables: dSM_dSF() -> dsm_dsf()

  - Stokes...Variables: viscosity() -> dynamicViscosity()

  - IAPWS water: ddgamma_ddpi, ddgamma_ddtau, ddgamma_dtaudpi, dgamma_dpi, 
    dgamma_dtau, dp_dpi, dpi_dp, dtau_dt 
    -> ddGamma_ddPi, ddGamma_ddTau, ddGamma_dTaudPi, dGamma_dPi, dGamma_dTau, 
    dp_dPi, dPi_dp, dTau_dt

* DELETED classes/files, property names, constants/enums, 
  member functions/variables, which have been deprecated in DuMuX 2.3:
  Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.2 and DuMuX 2.3
===================================================

* IMPROVEMENTS and ENHANCEMENTS:
  - A fully implicit two-phase discrete-fracture-matrix model has been added, 
    see test/implicit/2pdfm.

  - Almost all porous media fully implicit models now can either use a
    vertex-centered (box) or a cell-centered spatial discretization. The choice
    of the spatial discretization method is controlled by deriving the problem
    type tag either from BoxModel or CCModel. This allows for a uniform problem
    description, as long as the boundaryTypesAtPos and dirichletAtPos methods
    can be used. By evaluating the compile-time property ImplicitIsBox, it is 
    easily possible to separately handle the different discretizations inside 
    am common method. See the tests in test/implicit for examples. 
    Correspondingly, the directory structure has been adapted. 
    Old:             New:
    dumux/           dumux/
      boxmodels/       implicit/
        common/          common/
        1p/              box/
        1p2c/            cellcentered/
        2p/              1p/
        ...              ...
    test/            test/
      boxmodels/       implicit/
        1p/              1p/
          test_1p          test_box1p
        ...                test_cc1p
                         ...

  - A backend for the ISTL AMG solver has been included, based on the
    corresponding DUNE-PDELab backends. It can be used for the fully 
    implicit and the decoupled models, see test_*1pwithamg in 
    test/implicit/1p and test_impeswithamg in test/decoupled/2p. 
    DUNE-PDELab and possibly DUNE-ISTL have to be patched, see the file 
    README in the patches directory. 

  - The decoupled models have been parallelized, see test_impeswithamg in 
    test/decoupled/2p. They work in parallel only if the AMGBackend is used 
    as linear solver. No dynamic loadbalancing can be done yet. 

  - The MPNC model can use either the most wetting or the most non-wetting phase
    pressure as primary variable. This is controlled via the property 
    "PressureFormulation."

  - The table of available parameters has been improved, see 
    http://www.dumux.org/doxygen-stable/html-2.2/a00838.php

  - Improved handling of the conductive heat fluxes in the non-isothermal implicit  
    two-phase models, see the problem files in test/implicit/2p(2c)ni.

  - Introduced new selection of start/stop messages.

* IMMEDIATE INTERFACE CHANGES not allowing/requiring a deprecation period:
  - The property Salinity used in the BrineCO2FluidSystem 
    has been renamed to ProblemSalinity. 

  - The matrixHeatFlux(...) and boundaryMatrixHeatFlux(...) methods in the spatial
    parameters of nonisothermal implicit twophase models have been removed.
    Instead, the computation of the effective thermal conductivity has been sourced 
    out to the fluidmatrixinteractions in a separate file 
    dumux/material/fluidmatrixinteractions/thermalconductivitysomerton.hh, which
    can be exchanged. The spatial parameters file needs a method 
    thermalConductivitySolid(...), where the thermal conductivity of the solid
    material only is specified. The rest is computed in the respective
    flux variables.

* Deprecated CLASSES/FILES, to be removed after 2.3:
  - The following headers in dumux/boxmodels/ have been deprecated and forward 
    to the corresponding headers in dumux/implicit/box: 
    boxassembler.hh, boxelementvolumevariables.hh, boxlocalresidual.hh,
    boxpropertydefaults.hh, boxelementboundarytypes.hh, boxfvelementgeometry.hh,
    boxproperties.hh, intersectiontovertexbc.hh

  - All headers in the following subdirectories of dumux/boxmodels have been 
    deprecated and forward to the headers in the corresponding subdirectories 
    of dumux/implicit:
    1p, 1p2c, 2p, 2p2c, 2p2cni, 2pdfm, 2pni, 
    3p3c, 3p3cni, co2, co2ni, mpnc, richards

  - Some box-specific classes "Box..." in dumux/boxmodels/common could be 
    completely replaced by unified "Implicit..." classes in 
    dumux/implicit/common: 
    ...DarcyFluxVariables, ...darcyfluxvariables.hh 
    ...ForchheimerFluxVariables, ...forchheimerfluxvariables.hh 
    ...LocalJacobian, ...localjacobian.hh
    ...Model, ...model.hh
    ...PorousMediaProblem, ...porousmediaproblem.hh
    ...Problem, ...problem.hh
    ...VolumeVariables, ...volumevariables.hh

  - The box-specific spatial parameter classes BoxSpatialParams... in 
    dumux/material/boxspatialparams....hh have been deprecated in favor of 
    ImplicitSpatialParams... in dumux/material/implicitspatialparams....hh.

  - The GridCreatorheaders from dumux/common have been moved to dumux/io: 
    cubegridcreator.hh, dgfgridcreator.hh, simplexgridcreator.hh

* Deprecated PROPERTY NAMES, to be removed after 2.3: BEWARE: The compiler will
  not print any warning if a deprecated property name is used.
  - CompositionFromFugacitiesSolver has been renamed to Constraintsolver.

* Deprecated public MEMBER VARIABLES, to be removed after 2.3: BEWARE: The
  compiler will not print any warning if a deprecated public member variable 
  is used.
  - numFAP and numSCV in Box(CC)FVElementGeometry have been renamed to 
    numFap and numScv, respectively.

* Deprecated MEMBER FUNCTIONS, to be removed after 2.3:
  - boundaryMatrixHeatFlux, markVertexRed and relativeErrorVertex 
    from ImplicitSpatialParams, ImplicitAssembler and ImplicitModel, 
    respectively. In favor of using 
    thermalConductivitySolid (see above), markDofRed and relativeErrorDof, 
    respectively. 

* DELETED classes/files, property names, constants/enums, 
  member functions, which have been deprecated in DuMuX 2.2:
  Everything listed as deprecated below has been removed.


Differences Between DuMuX 2.1 and DuMuX 2.2
===================================================

* IMPROVEMENTS and ENHANCEMENTS: 
  - Two new fully implicit models dedicated to simulate compositional
    (non-isothermal) CO2-brine systems have been added, together with 
    corresponding components and a fluid system. See test/boxmodels/co2(ni) 
    for details. These tests also illustrate the usage of element and vertex 
    parameters as well as boundary ids provided by DGF files for setting 
    permeability and porosity as well as boundary conditions.

  - Decoupled Models: An h-adaptive model using an MPFA L-method was added
    that simulates 2p and 2p2c flow on unstructured grids with hanging nodes 
    in two dimensions. See test/decoupled/2p(2c) for details.

  - All fully implicit porous media models are now capable of employing 
    the Forchheimer equation as an alternative to the commonly used 
    Darcy law. See test_forchheimer*p in test/boxmodels/mpnc for details. 

  - The Stokes models are now able to simulate the full Navier-Stokes 
    equations for momentum transport. See test/freeflow/navierstokes 
    for details.

  - The fully implicit models have been (partially) generalized to allow 
    for a cell-centered discretization in addition to the default 
    vertex-centered (box) one. Cell-centered fully implicit 2p and 2p2c 
    models are already available in the developers part of Dumux. Further 
    generalizations and the inclusion in the stable part are planned for 
    Dumux 2.3.

  - Several model-specific features and classes have been unified, like 
    the calculation of the Darcy velocity for the fully implicit flux 
    variables, or the temperature, gravity, and spatial parameter 
    functionalities of the fully implicit problems. Moreover, many
    names have been made more consistent. This includes the naming 
    and grouping of several parameters and corresponding properties, 
    the indexing of phases and components, and the preference of the 
    partial name "params" over "parameters." For details, see also the 
    deprecations listed below.

  - Added compiler support for GCC 4.7 and Clang 3.1.

* IMMEDIATE INTERFACE CHANGES not allowing a deprecation period:
  - From Dune 2.2 on, FieldVector::size is a method rather than an enum value.
    It is mandatory to add the flag --enable-fieldvector-size-is-method to the
    CONFIGURE_FLAGS. An example is given in the opts file dumux/debug.opts.
  - Implicit models: TwoPIndices, TwoPNIIndices, and RichardsIndices 
    additionally get TypeTag as template parameter. If the Indices are not 
    obtained via the property, this has to be adapted.

  - Implicit models: All model-specific computeFlux functions in 
    ...localresidual.hh have to get an additional bool parameter onBoundary, 
    which is by default set to false. If outflow conditions should 
    be properly implemented, also the constructor of the flux variables in 
    ...fluxvariables.hh has to get the additional argument and the 
    class has to be adapted to deal with boundary faces. See FS#117 and #99
    for details.
    
* Deprecated CLASSES/FILES, to be removed after 2.2:
  - Model specific base box problems: The common functionality has been 
    collected in PorousMediaBoxProblem in 
    dumux/boxmodels/common/porousmediaboxproblem.hh. The problem can be derived 
    from PorousMediaBoxProblem, instead of the model specific base problem: 
    OnePBoxProblem, dumux/boxmodels/1p/1pproblem.hh,
    OnePTwoCBoxProblem, dumux/boxmodels/1p2c/1p2cproblem.hh,
    TwoPProblem, dumux/boxmodels/2p/2pproblem.hh,
    TwoPNIProblem, dumux/boxmodels/2pni/2pniproblem.hh,
    TwoPTwoCProblem, dumux/boxmodels/2p2c/2p2cproblem.hh,
    TwoPTwoCNIProblem, dumux/boxmodels/2p2cni/2p2cniproblem.hh,
    ThreePThreeCProblem, dumux/boxmodels/3p3c/3p3cproblem.hh,
    ThreePThreeCNIProblem, dumux/boxmodels/3p3cni/3p3cniproblem.hh,
    MPNCProblem, dumux/boxmodels/mpnc/mpncproblem.hh.
    
  - All "...SpatialParameters" base classes have been replaced by 
    "...SpatialParams" classes:
    BoxSpatialParameters, dumux/material/spatialparameters/boxspatialparameters.hh,
    BoxSpatialParametersOneP, dumux/material/spatialparameters/boxspatialparameters1p.hh,
    FVSpatialParameters, dumux/material/spatialparameters/fvspatialparameters.hh,
    FVSpatialParametersOneP, dumux/material/spatialparameters/fvspatialparameters1p.hh.

  - Due to the unification of flux variables for the fully implicit models, 
    some model-specific flux variables have become obsolete: 
    OnePFluxVariables, dumux/boxmodels/1p/1pfluxvariables.hh,
    TwoPFluxVariables, dumux/boxmodels/2p/2pfluxvariables.hh,
    RichardsFluxVariables, dumux/boxmodels/richards/richardsfluxvariables.hh.

  - Two components have new names and locations in dumux/material/components: 
    SimpleDNAPL, simplednapl.hh -> DNAPL, napl.hh
    Oil, oil.hh -> LNAPL, lnapl.hh

  - Some MPFA-O method files/classes have been moved to a new subdirectory 
    "omethod" in dumux/decoupled/2p/diffusion/fvmpfa: 
    fvmpfaopressure2p.hh, fvmpfaovelocity2p.hh, fvmpfaopressureproperties2p.hh

  - DUMUX_UNUSED is deprecated and will be removed after 2.2. It should be 
    replaced by the upstream version DUNE_UNUSED.
  
  - DUMUX_DEPRECATED_MSG is deprecated and will be removed after 2.2. It should
    be replaced by the upstream version DUNE_DEPRECATED_MSG.

* Deprecated PROPERTY NAMES, to be removed after 2.2: BEWARE: The compiler will
  not print any warning if a deprecated property name is used. 
  - The "SpatialParameters" property has been renamed to "SpatialParams".
  
  - The model specific "...Indices" property has been renamed to "Indices".

* Deprecated CONSTANTS/ENUMS, to be removed after 2.2: BEWARE: The compiler will
  not print any warning if a deprecated constant/enum is used.
  - In the 2p2c/ni and 3p3c/ni models, all indices related to phase and 
    components can be pre/suffixed with "w", "n" and, 
    for three phases, with "g". 
    boxmodels/2p2c/...: "l", "g" pre/suffixes have been replaced by "w", "n".
    boxmodels/3p3c/...: "c", "a" pre/suffixes have been replaced by "n", "g".

* Deprecated MEMBER FUNCTIONS, to be removed after 2.2:
  - Spatial parameters: The spatialParameters member functions of the base 
    problems have been replaced by spatialParams: 
    dumux/boxmodels/common/porousmediaboxproblem.hh, 
    dumux/decoupled/1p/diffusion/diffusionproblem...hh,
    dumux/decoupled/2p/impes/impesproblem2p.hh,
    dumux/decoupled/2p/transport/transportproblem2p.hh.
  
  - Flux variables: Renaming of members 
    "...AtIP" -> "...", 
    "concentration..." -> "massFraction...",
    "molarConc..." -> "moleFraction..."
    The "massFraction..." members have been deprecated, instead 
    "moleFraction..." should be used.
    Affected files:
    dumux/boxmodels/1p2c/1p2cfluxvariables.hh, 
    dumux/boxmodels/2p2c/2p2cfluxvariables.hh,
    dumux/boxmodels/mpnc/.../...fluxvariables.hh,
    dumux/freeflow/stokes.../stokes...fluxvariables.hh.

  - Box models: The primaryVarWeight() functions are no longer used for the 
    evaluation of the relative error. 

  - Element and FVElementGeometry: The elem_() and fvElemGeom_() member function
    of BoxLocalResidual have been replaced by element_() and fvGeometry_(). 
    
  - Primary variables: All "...primaryVar/s" member functions have been replaced
    by "...priVar/s": 
    dumux/boxmodels/common/boxlocalresidual.hh, 
    dumux/boxmodels/common/boxvolumevariables.hh.

  - Start functionality in dumux/common/start.hh: printUsageDGF and 
    printUsageGrid are no longer needed.
  
* DELETED member functions, which have been deprecated in DuMuX 2.1: 
  - dumux/material/spatialparameters/boxspatialparameters1p.hh:
    extrusionFactorScv and extrusionFactorScvf, now part of the volume variables

  - dumux/material/idealgas.hh:
    concentration, replaced by molarDensity

  - dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh:
    pC(const Params &params, Scalar Sw, const Scalar temperature)

  - dumux/common/start.hh:
    startFromDGF, startWithGrid, startWithParameters, all replaced by start

  - dumux/common/spline.hh:
    set(const ScalarArray&, const ScalarArray&, Scalar, Scalar), replaced by setXYArrays,
    set(const PointArray&, Scalar, Scalar), replaced by setArrayOfPoints

  - dumux/common/variablelengthspline_.hh, dumux/common/fixedlengthspline_.hh:
    various set routines, replaced by more descriptive names

  - dumux/io/vtkmultiwriter.hh:
    VtkMultiWriter(const std::string&, std::string), replaced by VtkMultiWriter(const GridView&, ...),
    beginTimeStep, replaced by beginWrite,
    createField, replaced by allocateManagedBuffer,
    addVertexData, replaced by attachVertexData,
    addCellData, replaced by attachCellData,
    endTimeStep, replaced by endWrite

  - dumux/decoupled/2p2c/2p2cproblem.hh:
    IMPETProblem2P2C(const GridView&, bool) replaced by IMPETProblem2P2C(TimeManager&, ...),
    IMPETProblem2P2C(..., SpatialParameters&, ...) replaced by IMPETProblem2P2C(TimeManager&, ...),
    initSat(const GlobalPosition&, const Element&) replaced by initSat(const Element&)
    initConcentration(const GlobalPosition&, const Element&) replaced by initConcentration(const Element&)
    
  - DUMUX_DEPRECATED has been removed.


Notable Differences Between DuMuX 2.0 and DuMuX 2.1
===================================================

- The thermodynamics framework has been overhauled:
  - The programming interfaces for fluid systems, fluid states and
    components has been formalized and cleaned up.
  - Fluid systems now have the option to cache computationally
    expensive parameters if they are needed for several relations.
  - Fluid systems are not charged with the computation of the
    chemical equilibrium anymore.
  - Fluid states are now centralized infrastructure instead of being
    model-specific.
  - Constraint solvers (which simplify solving thermodynamic
    constraints) have been introduced.
- Outflow boundary conditions have been implemented for the
  fully-implicit models 1p2c, 2p2c(ni) and stokes(2cni).
- The problem and spatial parameter base classes also provide optional
  model-independent interfaces. These methods only get the position in
  global coordinates as argument and are named *AtPos()
  (e.g. boundaryTypesAtPos()). This allows an easy transfer of problem
  definitions between implicit and sequential models.
- The following fully-implicit models have been added:
  - 3p3c, 3p3cni: Isothermal and non-isothermal three-phase,
    three-component models for flow and transport in porous media
    based on primary variable switching.
  - MpNc: A model for arbitrary number of phases M > 0, and components
    (N >= M - 1 >= 1) for flow and transport in porous media. This
    model also comes with an energy and a molecular diffusion module.
  - stokes, stokes2c, stokes2cni: Models for the plain Stokes
    equation as well as isothermal and non-isothermal Stokes models
    for two-component fluids.
- The sequentially-coupled models have been overhauled:
  - A common structure for cell centered standard finite volume
    implementations has been introduced.
  - The data structures where overhauled to avoid large clumps of data
    in large-scale simulations: Each cell stores data in its own
    storage object.
  - The too large assemble() methods have been split into submethods
    getStorage(), getFlux() etc. By this, inheritance of classes has
    been improved and code duplication was reduced.
  - The conceptual seperation of the "VariableClass" (central
    infrastructure), data storage, transport model and pressure model
    has been improved.
  - More of infrastructure is now shared with the implicit models
    (e.g. the BoundaryTypes). This results in significant performance
    improvements and makes debugging easier.
- The 2padaptive sequentially coupled model has been added. This model
  implements a grid-adaptive finite volume scheme for immiscible
  two-phase flow in porous media on non-conforming quadrilateral
  grids.
- The dependencies for the external dune-pdelab and boost packages
  have been removed.
- The build system has received major improvements:
  - There is now much better test coverage of build-time dependencies
    on packages for the default autotools-based build system.
  - Experimental support for building DuMuX using CMake has been much
    improved. In the long run, CMake is projected to become the
    default build system.
- All headers can now be included without any preconditions.
- DuMuX now compiles without warnings if the -pedantic flag used for GCC.
- Specifying run-time parameters is now possible. The mechanism allows
  to use parameter files or to specify parameters directly on the
  command line and fallback parameter input files have been added for
  each test application.  As a consequence, applications can be run
  now without specifying any command line arguments.
- The DuMuX property system has been fine-tuned:
  - Encapsulating property names with the PTAG() is no longer required
    for the GET_PROP* macros (but is still allowed).
  - Setting property defaults has been deprecated.
  - All properties defined for a type tag can now be printed. Also,
    their value and the location in the source where they where
    specified is included in the output.
- Using quadruple precision math has been made possible for GCC 4.6 or newer:
  - To use it, the configure option '--enable-quad' needs to be added
    and the type of scalar values needs to be changed to 'quad'. This
    can be done in the problem file using

        SET_TYPE_PROP(YourProblemTypeTag, Scalar, quad);

    It should be noted, that performance is very poor when using
    quadruple precision arithmetic. This feature is primarily meant as
    a debugging tool to quickly check whether there are machine
    precision related convergence problems.
