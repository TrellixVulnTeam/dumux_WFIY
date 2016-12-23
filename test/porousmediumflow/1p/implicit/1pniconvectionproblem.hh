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
 * \brief Test for the OnePModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 */
#ifndef DUMUX_1PNI_CONVECTION_PROBLEM_HH
#define DUMUX_1PNI_CONVECTION_PROBLEM_HH

#include <math.h>

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1pnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePNIConvectionProblem;

namespace Properties
{

NEW_TYPE_TAG(OnePNIConvectionProblem, INHERITS_FROM(OnePNI));
NEW_TYPE_TAG(OnePNIConvectionBoxProblem, INHERITS_FROM(BoxModel, OnePNIConvectionProblem));
NEW_TYPE_TAG(OnePNIConvectionCCProblem, INHERITS_FROM(CCTpfaModel, OnePNIConvectionProblem));
NEW_TYPE_TAG(OnePNIConvectionCCMpfaProblem, INHERITS_FROM(CCMpfaModel, OnePNIConvectionProblem));

// Set the grid type
SET_TYPE_PROP(OnePNIConvectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(OnePNIConvectionProblem, Problem, OnePNIConvectionProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(OnePNIConvectionProblem, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                        H2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Set the spatial parameters
SET_TYPE_PROP(OnePNIConvectionProblem, SpatialParams, OnePNISpatialParams<TypeTag>);

} // end namespace Properties


/*!
 * \ingroup OnePNIModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the OnePModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side water is injected at a constant rate and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution where a retarded front velocity is calculated as follows:
  \f[
     v_{Front}=\frac{q S_{water}}{\phi S_{total}}
 \f]
 *
 * The result of the analytical solution is written into the vtu files.
 * This problem uses the \ref OnePModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1pniconvection -ParameterFile ./test_box1pniconvection.input</tt> or <br>
 * <tt>./test_cc1pniconvection -ParameterFile ./test_cc1pniconvection.input</tt>
 */
template <class TypeTag>
class OnePNIConvectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using IapwsH2O = H2O<Scalar>;

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
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };


    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    OnePNIConvectionProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);
        darcyVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, DarcyVelocity);

        temperatureHigh_ = 291.0;
        temperatureLow_ = 290.0;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeOver() ||
            this->timeManager().willBeFinished();
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
    {
        //Here we calculate the analytical solution
        auto numDofs = this->model().numDofs();

        //create required scalar fields
        auto& temperatureExact = *(this->resultWriter().allocateManagedBuffer(numDofs));
        auto& temperature = *(this->resultWriter().allocateManagedBuffer(numDofs));

        const auto someElement = *(elements(this->gridView()).begin());
        const auto someElemSol = this->model().elementSolution(someElement, this->model().curSol());

        auto someFvGeometry = localView(this->model().globalFvGeometry());
        someFvGeometry.bindElement(someElement);
        const auto someScv = *(scvs(someFvGeometry).begin());

        VolumeVariables volVars;
        volVars.update(someElemSol, *this, someElement, someScv);

        const auto porosity = this->spatialParams().porosity(someElement, someScv, someElemSol);
        const auto densityW = volVars.density();
        const auto heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
        const auto storageW =  densityW*heatCapacityW*porosity;
        const auto densityS = this->spatialParams().solidDensity(someElement, someScv, someElemSol);
        const auto heatCapacityS = this->spatialParams().solidHeatCapacity(someElement, someScv, someElemSol);
        const auto storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);

        std::cout << "storage: " << storageTotal << std::endl;

        const auto time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);
        const auto retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;

        std::cout << "retarded velocity: " << retardedFrontVelocity << std::endl;

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto globalIdx = scv.dofIndex();
                const auto& globalPos = scv.dofPosition();
                temperatureExact[globalIdx] = globalPos[0] < retardedFrontVelocity*time ? temperatureHigh_ : temperatureLow_;
                temperature[globalIdx] = this->model().curSol()[globalIdx][temperatureIdx];
            }
        }
        this->resultWriter().attachDofData(temperature, "temperature", isBox);
        this->resultWriter().attachDofData(temperatureExact, "temperatureExact", isBox);
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

        if(globalPos[0] > this->bBoxMax()[0] - eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolvars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        const auto globalPos = scvf.ipGlobal();
        const auto& volVars = elemVolvars[scvf.insideScvIdx()];

        if(globalPos[0] < eps_)
        {
            values[conti0EqIdx] = -darcyVelocity_*volVars.density();
            values[energyEqIdx] = values[conti0EqIdx]*IapwsH2O::liquidEnthalpy(temperatureHigh_, volVars.pressure());
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0);
    }

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
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[temperatureIdx] = temperatureLow_;
        return priVars;
    }

    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    int outputInterval_;
};

} //end namespace Dumux

#endif // DUMUX_1PNI_CONVECTION_PROBLEM_HH
