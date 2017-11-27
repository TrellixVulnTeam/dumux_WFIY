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
/*!
 * \file
 *
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_INFILTRATION_THREEP_PROBLEM_HH
#define DUMUX_INFILTRATION_THREEP_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/3p/implicit/model.hh>
#include <dumux/material/fluidsystems/h2oairmesitylene.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include "infiltration3pspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class InfiltrationThreePProblem;

namespace Properties
{
NEW_TYPE_TAG(InfiltrationThreePProblem, INHERITS_FROM(ThreeP, InfiltrationThreePSpatialParams));
NEW_TYPE_TAG(InfiltrationThreePBoxProblem, INHERITS_FROM(BoxModel, InfiltrationThreePProblem));
NEW_TYPE_TAG(InfiltrationThreePCCProblem, INHERITS_FROM(CCTpfaModel, InfiltrationThreePProblem));

// Set the grid type
SET_TYPE_PROP(InfiltrationThreePProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationThreePProblem, Problem, InfiltrationThreePProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(InfiltrationThreePProblem,
              FluidSystem,
              FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

}// end namespace Properties

/*!
 * \ingroup ThreePModel
 * \ingroup ImplicitTestProblems
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone.
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file naplinfiltrationexercise.input
 *
 * To run the simulation execute the following line in shell:
 * <tt>./naplinfiltrationexercise -parameterFile naplinfiltrationexercise.input</tt>
 *  */
template <class TypeTag >
class InfiltrationThreePProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        pressureIdx = Indices::pressureIdx,
        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InfiltrationThreePProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(282.15, 284.15, 3, 8e4, 3e5, 200);

        name_ = getParam<std::string>("Problem.Name");

        this->spatialParams().plotMaterialLaw();
        time_ = 0.0;
    }

    /*!
     * \name Problem parameters
     */
    // \{
    void setTime(Scalar time)
    { time_ = time; }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_;
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
        BoundaryTypes values;

        if(globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_
           || globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
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
        PrimaryVariables values(0.0);
        initial_(values, globalPos);
        return values;
    }
    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NeumannFluxes neumannAtPos(const GlobalPosition &globalPos) const
    {
        NeumannFluxes values(0.0);

        // negative values for injection
        if (time_ < 2592000.0 - eps_)
        {
            if ((globalPos[0] < 175.0 + eps_) && (globalPos[0] > 155.0 - eps_)
                 && (globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_))
            {
                // mol fluxes, convert with M(Mesit.)=0,120 kg/mol --> 1.2e-4  kg/(sm)
                values[Indices::contiNEqIdx] = -0.001;
            }
        }

        return values;
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
        PrimaryVariables values(0.0);
        initial_(values, globalPos);
        return values;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a uniform temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }



private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar sw, swr=0.12, sgr=0.03;

        if(y > (-1.E-3*x+5) - eps_)
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            sw = invertPcgw_(pc,
                             this->spatialParams().materialLawParamsAtPos(globalPos));
            if (sw < swr) sw = swr;
            if (sw > 1.-sgr) sw = 1.-sgr;

            values[pressureIdx] = 1e5 ;
            values[swIdx] = sw;
            values[snIdx] = 0.;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[swIdx] = 1.-sgr;
            values[snIdx] = 0.;
        }
    }

    // small solver inverting the pc curve
    static Scalar invertPcgw_(Scalar pcIn, const MaterialLawParams &pcParams)
    {
        Scalar lower,upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar sw, pcgw;
        lower=0.0; upper=1.0;
        for (k=1; k<=25; k++)
        {
            sw = 0.5*(upper+lower);
            pcgw = MaterialLaw::pcgw(pcParams, sw);
            Scalar delta = pcgw-pcIn;
            if (delta<0.) delta*=-1.;
            if (delta<bisLimit)
            {
                return(sw);
            }
            if (k==maxIt) {
                return(sw);
            }
            if (pcgw>pcIn) lower=sw;
            else upper=sw;
        }
        return(sw);
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    Scalar time_;
};

} // end namespace Dumux

#endif
