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
 * \brief test problem for the sequential 2p2c model
 */
#ifndef DUMUX_TEST_2P2C_PROBLEM_HH
#define DUMUX_TEST_2P2C_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dumux/common/cubegridcreator.hh>

#include <dumux/decoupled/2p2c/2p2cproblem.hh>
#include <dumux/decoupled/2p2c/fvpressure2p2c.hh>
#include <dumux/decoupled/2p2c/fvtransport2p2c.hh>

// fluid properties
//#include <dumux/material/old_fluidsystems/simple_h2o_n2_system.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include "test_dec2p2c_spatialparams.hh"
#include <dumux/linear/impetbicgstabilu0solver.hh>
namespace Dumux
{

template<class TypeTag>
class TestDecTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
NEW_TYPE_TAG(TestDecTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC, Test2P2CSpatialParams));

// set the GridCreator property
SET_TYPE_PROP(TestDecTwoPTwoCProblem, GridCreator, CubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(TestDecTwoPTwoCProblem, Grid)
{
        typedef Dune::YaspGrid<3> type;
//    typedef Dune::SGrid<3, 3> type;
};

// Set the problem property
SET_PROP(TestDecTwoPTwoCProblem, Problem)
{
    typedef Dumux::TestDecTwoPTwoCProblem<TypeTag> type;
};

// Set the model properties
SET_PROP(TestDecTwoPTwoCProblem, TransportModel)
{
    typedef Dumux::FVTransport2P2C<TypeTag> type;
};

SET_PROP(TestDecTwoPTwoCProblem, PressureModel)
{
    typedef Dumux::FVPressure2P2C<TypeTag> type;
};

SET_INT_PROP(TestDecTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureNW);

// Select fluid system
SET_PROP(TestDecTwoPTwoCProblem, FluidSystem)
{
    typedef Dumux::H2OAirFluidSystem<TypeTag> type;
};
// Select fluid system
SET_BOOL_PROP(TestDecTwoPTwoCProblem, EnableComplicatedFluidSystem, true);


// Select water formulation
SET_PROP(TestDecTwoPTwoCProblem, Components) : public GET_PROP(TypeTag, DefaultComponents)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//    typedef Dumux::TabulatedComponent<Scalar, typename Dumux::H2O<Scalar> > H20;
        typedef Dumux::H2O<Scalar> H2O;
};

//SET_TYPE_PROP(TestDecTwoPTwoCProblem, LinearSolver, IMPETBiCGStabILU0Solver<TypeTag> );

// Enable gravity
SET_BOOL_PROP(TestDecTwoPTwoCProblem, ProblemEnableGravity, true);
SET_BOOL_PROP(TestDecTwoPTwoCProblem, EnableCapillarity, true);
SET_INT_PROP(TestDecTwoPTwoCProblem,
        BoundaryMobility,
        GET_PROP_TYPE(TypeTag, Indices)::satDependent);
SET_SCALAR_PROP(TestDecTwoPTwoCProblem, CFLFactor, 0.8);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p2c model
 *
 * The domain is box shaped (3D). All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Nitrogen)
 * is injected over a vertical well in the center of the domain.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_dec2p2c -parameterFile ./test_dec2p2c.input</tt>
 */
template<class TypeTag>
class TestDecTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
TestDecTwoPTwoCProblem(TimeManager &timeManager, const GridView &gridView) :
ParentType(timeManager, gridView), eps_(1e-6), depthBOR_(1000.0)
{
//    this->setOutputInterval(20);
    // initialize the tables of the fluid system
//    FluidSystem::init(/*tempMin=*/280,
//            /*tempMax=*/290,
//            /*numTemp=*/10,
//            /*pMin=*/190000,
//            /*pMax=*/280000,
//            /*numP=*/400);
}

/*!
 * \name Problem parameters
 */
// \{

//! The problem name.
/*! This is used as a prefix for files generated by the simulation.
*/
const char *name() const
{
    return "test_dec2p2c";
}
//!  Returns true if a restart file should be written.
/* The default behaviour is to write no restart file.
 */
bool shouldWriteRestartFile() const
{
    return false;
}

//! Returns the temperature within the domain.
/*! This problem assumes a temperature of 10 degrees Celsius.
 * \param globalPos The global Position
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
//! Returns the reference pressure.
 /*This pressure is used in order to calculate the material properties
 * at the beginning of the initialization routine. It should lie within
 * a reasonable pressure range for the current problem.
 * \param globalPos The global Position
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
//! Type of boundary condition.
/*! Defines the type the boundary condition for the conservation equation,
 *  e.g. Dirichlet or Neumann. The Pressure equation is acessed via
 *  Indices::pressureEqIdx, while the transport (or mass conservation -)
 *  equations are reached with Indices::contiWEqIdx and Indices::contiNEqIdx
 *
 * \param bcTypes The boundary types for the conservation equations
 * \param globalPos The global Position
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > this->bboxMax()[0]-1E-6 || globalPos[0] < 1e-6)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

//! Flag for the type of Dirichlet conditions
/*! The Dirichlet BCs can be specified by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 *
 * \param bcFormulation The boundary formulation for the conservation equations.
 * \param intersection The intersection on the boundary.
 */
const void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
//! Values for dirichlet boundary condition \f$ [Pa] \f$ for pressure and \f$ \frac{mass}{totalmass} \f$ or \f$ S_{\alpha} \f$ for transport.
/*! In case of a dirichlet BC, values for all primary variables have to be set. In the sequential 2p2c model, a pressure
 * is required for the pressure equation and component fractions for the transport equations. Although one BC for the two
 * transport equations can be deduced by the other, it is seperately defined for consistency reasons with other models.
 * Depending on the boundary Formulation, either saturation or total mass fractions can be defined.
 *
 * \param bcValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
 */
void dirichletAtPos(PrimaryVariables &bcValues ,const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    // Dirichlet for pressure equation
    bcValues[Indices::pressureEqIdx] = (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);

    // Dirichlet values for transport equations
    bcValues[Indices::contiWEqIdx] = 1.;
    bcValues[Indices::contiNEqIdx] = 1.- bcValues[Indices::contiWEqIdx];

}

//! Value for neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$.
/*! In case of a neumann boundary condition, the flux of matter
 *  is returned. Both pressure and transport module regard the neumann-bc values
 *  with the tranport (mass conservation -) equation indices. So the first entry
 *  of the vector is superflous.
 *  An influx into the domain has negative sign.
 *
 * \param neumannValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    neumannValues[Indices::contiNEqIdx] =0.;
    neumannValues[Indices::contiWEqIdx] =0.;
//    if (globalPos[1] < 15 && globalPos[1]> 5)
//    {
//        neumannValues[Indices::contiNEqIdx] = -0.015;
//    }
}
//! Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$
/*! Evaluate the source term for all phases within a given
 *  volume. The method returns the mass generated (positive) or
 *  annihilated (negative) per volume unit.
 *  Both pressure and transport module regard the neumann-bc values
 *  with the tranport (mass conservation -) equation indices. So the first entry
 *  of the vector is superflous.
 *
 * \param sourceValues Vector holding values for all equations (pressure and transport).
 * \param globalPos The global Position
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    sourceValues[Indices::contiWEqIdx]=0.;
    sourceValues[Indices::contiNEqIdx]=0.;
    if (fabs(globalPos[0] - 4.8) < 0.5 && fabs(globalPos[1] - 4.8) < 0.5)
        sourceValues[Indices::contiNEqIdx] = 0.0001;
}
//! Flag for the type of initial conditions
/*! The problem can be initialized by a given concentration (mass of
 * a component per total mass inside the control volume) or by means
 * of a saturation.
 */
const void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
//! Concentration initial condition (dimensionless)
/*! The problem is initialized with the following concentration.
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1;
}
private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

const Scalar eps_;
const Scalar depthBOR_;
};
} //end namespace

#endif
