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
/**
 * \file
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 */
#ifndef DUMUX_1PNI_CONDUCTION_PROBLEM_HH
#define DUMUX_1PNI_CONDUCTION_PROBLEM_HH

#include <math.h>

#if PROBLEM==1
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#elif PROBLEM==2
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#elif PROBLEM==3
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#endif
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include "1pnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePNIConductionProblem;

namespace Properties
{
#if PROBLEM==1
NEW_TYPE_TAG(OnePNIConductionProblem, INHERITS_FROM(OnePNIMimetic));
#elif PROBLEM==2
NEW_TYPE_TAG(OnePNIConductionProblem, INHERITS_FROM(CCTpfaModel, OnePNI));
#elif PROBLEM==3
NEW_TYPE_TAG(OnePNIConductionProblem, INHERITS_FROM(CCMpfaModel, OnePNI));
#endif

// Set the grid type
SET_TYPE_PROP(OnePNIConductionProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePNIConductionProblem, Problem,
              OnePNIConductionProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(OnePNIConductionProblem, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                        H2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
// Set the spatial parameters
SET_TYPE_PROP(OnePNIConductionProblem,
              SpatialParams,
              OnePNISpatialParams<TypeTag>);

SET_TYPE_PROP(OnePNIConductionProblem, LinearSolver, SuperLUBackend<TypeTag> );

SET_BOOL_PROP(OnePNIConductionProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePNIConductionProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePNIConductionProblem, EnableGlobalVolumeVariablesCache, true);

}


/*!
 * \ingroup OnePNIModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the OnePModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side there is a Dirichlet boundary condition with an increased temperature and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution for a diffusion process:
  \f[
     T =T_{high} + (T_{init} - T_{high})erf \left(0.5\sqrt{\frac{x^2 S_{total}}{t \lambda_{eff}}}\right)
 \f]
 *
 * This problem uses the \ref OnePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1pniconduction -ParameterFile ./test_box1pniconduction.input</tt> or <br>
 * <tt>./test_cc1pniconduction -ParameterFile ./test_cc1pniconduction.input</tt>
 */
template <class TypeTag>
class OnePNIConductionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using VtkOutputModule = typename GET_PROP_TYPE(TypeTag, VtkOutputModule);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dimWorld : 0 };

    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx,
#if PROBLEM==1
        facePressureIdx = Indices::facePressureIdx,
        faceTemperatureIdx = Indices::faceTemperatureIdx
#else
        facePressureIdx = Indices::pressureIdx,
        faceTemperatureIdx = Indices::temperatureIdx,
#endif
    };

    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    OnePNIConductionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);

        temperatureHigh_ = 300.0;

        pIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureIn);
        pOut_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureOut);
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeFinished() ||
            this->timeManager().willBeFinished();
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 5 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& temperatureExact = outputModule.createScalarField("temperatureExact", dofCodim);

        const auto someElement = *(elements(this->gridView()).begin());
        const auto someElemSol = this->model().elementSolution(someElement, this->model().curSol());
        const auto someInitSol = initial_(someElement.geometry().center());

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storage = densityW*heatCapacityW*porosity + densityS*heatCapacityS*(1 - porosity);
        const auto effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars, this->spatialParams(),
                                                                                                         someElement, someFvGeometry, someScv);
        using std::max;
        Scalar time = max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto globalIdx = scv.dofIndex();
                const auto& globalPos = scv.dofPosition();

                using std::erf; using std::sqrt;
                temperatureExact[globalIdx] = temperatureHigh_ + (someInitSol[temperatureIdx] - temperatureHigh_)
                                              *erf(0.5*sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));
            }
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        if(globalPos[1] < eps_ || globalPos[1] > this->bBoxMax()[1] - eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(initial_(globalPos));
        if (globalPos[1] < eps_)
        {
            priVars[pressureIdx] = pIn_;
            priVars[temperatureIdx] = temperatureHigh_;
            priVars[facePressureIdx] = pIn_;
            priVars[faceTemperatureIdx] = temperatureHigh_;
        }

        return priVars;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = pOut_; // initial condition for the pressure
        priVars[temperatureIdx] = 290.0;
        priVars[facePressureIdx] = pOut_;
        priVars[faceTemperatureIdx] = 290.0;
        return priVars;
    }

    Scalar pIn_;
    Scalar pOut_;
    Scalar temperatureHigh_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    int outputInterval_;
};

} //end namespace Dumux

#endif // DUMUX_1PNI_CONDUCTION_PROBLEM_HH
