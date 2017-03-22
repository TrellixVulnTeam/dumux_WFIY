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
 * \brief problem for the sequential tutorial
 */
#ifndef DUMUX_EX2TUTORIALPROBLEM_SEQUENTIAL_HH // guardian macro /*@\label{tutorial-sequential:guardian1}@*/
#define DUMUX_EX2TUTORIALPROBLEM_SEQUENTIAL_HH // guardian macro /*@\label{tutorial-sequential:guardian2}@*/

// dumux 2p-sequential environment
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/impes/problem.hh> /*@\label{tutorial-sequential:parent-problem}@*/

// assign parameters dependent on space (e.g. spatial parameters)
#include "ex2tutorialspatialparams_sequential.hh" /*@\label{tutorial-sequential:spatialparameters}@*/

// include cfl-criterion after coats: more suitable if the problem is not advection dominated
#include<dumux/porousmediumflow/2p/sequential/transport/cellcentered/evalcflfluxcoats.hh>

// the components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/lnapl.hh>

namespace Dumux
{

template<class TypeTag>
class Ex2TutorialProblemSequential;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
// create a new type tag for the problem
NEW_TYPE_TAG(Ex2TutorialProblemSequential, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP,
                                                        Ex2TutorialSpatialParamsSequential)); /*@\label{tutorial-sequential:create-type-tag}@*/

// Set the problem property
SET_PROP(Ex2TutorialProblemSequential, Problem) /*@\label{tutorial-sequential:set-problem}@*/
{
    typedef Ex2TutorialProblemSequential<TypeTag> type;
};

// Set the grid type
SET_TYPE_PROP(Ex2TutorialProblemSequential, Grid, Dune::YaspGrid<2>); /*@\label{tutorial-sequential:set-grid-type}@*/

// Set the wetting phase
SET_PROP(Ex2TutorialProblemSequential, WettingPhase) /*@\label{tutorial-sequential:2p-system-start}@*/
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, H2O<Scalar> > type; /*@\label{tutorial-sequential:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(Ex2TutorialProblemSequential, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, LNAPL<Scalar> > type; /*@\label{tutorial-sequential:nonwettingPhase}@*/
}; /*@\label{tutorial-sequential:2p-system-end}@*/

SET_TYPE_PROP(Ex2TutorialProblemSequential, EvalCflFluxFunction, EvalCflFluxCoats<TypeTag>); /*@\label{tutorial-sequential:cflflux}@*/
SET_SCALAR_PROP(Ex2TutorialProblemSequential, ImpetCFLFactor, 0.95); /*@\label{tutorial-sequential:cflfactor}@*/

// Disable gravity
SET_BOOL_PROP(Ex2TutorialProblemSequential, ProblemEnableGravity, false); /*@\label{tutorial-sequential:gravity}@*/
} /*@\label{tutorial-sequential:propertysystem-end}@*/

/*! \ingroup SequentialProblems
 * @brief Problem class for the sequential tutorial
*/
template<class TypeTag>
class Ex2TutorialProblemSequential: public IMPESProblem2P<TypeTag> /*@\label{tutorial-sequential:def-problem}@*/
{
    typedef IMPESProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    enum
    {
        dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pwIdx = Indices::pwIdx,
        swIdx = Indices::swIdx,
        pressEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    Ex2TutorialProblemSequential(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-6)/*@\label{tutorial-sequential:constructor-problem}@*/
    {
        //write only every 10th time step to output file
        this->setOutputInterval(10);/*@\label{tutorial-sequential:outputinterval}@*/
    }

    //! The problem name.
    /*! This is used as a prefix for files generated by the simulation.
    */
    std::string name() const    /*@\label{tutorial-sequential:name}@*/
    {
        return "tutorial_sequential";
    }

    //!  Returns true if a restart file should be written.
    /* The default behaviour is to write no restart file.
     */
    bool shouldWriteRestartFile() const /*@\label{tutorial-sequential:restart}@*/
    {
        return false;
    }

    //! Returns the temperature within the domain at position globalPos.
    /*! This problem assumes a temperature of 10 degrees Celsius.
     *
     *  \param element The finite volume element
     *
     * Alternatively, the function temperatureAtPos(const GlobalPosition& globalPos) could be
     * defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    Scalar temperature(const Element& element) const /*@\label{tutorial-sequential:temperature}@*/
    {
        return 273.15 + 10; // -> 10°C
    }

    //! Returns a constant pressure to enter material laws at position globalPos.
    /* For incrompressible simulations, a constant pressure is necessary
     * to enter the material laws to gain a constant density etc. In the compressible
     * case, the pressure is used for the initialization of material laws.
     *
     * \param element The finite volume element
     *
     * Alternatively, the function referencePressureAtPos(const GlobalPosition& globalPos) could be
     * defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    Scalar referencePressure(const Element& element) const /*@\label{tutorial-sequential:refPressure}@*/
    {
        return 2e5;
    }

    //! Source of mass \f$ [\frac{kg}{m^3 \cdot s}] \f$ of a finite volume.
    /*! Evaluate the source term for all phases within a given
     *  volume.
     *
     *  \param values Includes sources for the two phases
     *  \param element The finite volume element
     *
     *  The method returns the mass generated (positive) or
     *  annihilated (negative) per volume unit.
     *
     * Alternatively, the function sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos)
     * could be defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    void source(PrimaryVariables &values, const Element& element) const /*@\label{tutorial-sequential:source}@*/
    {
        values = 0;
    }

    //! Type of boundary conditions at position globalPos.
    /*! Defines the type the boundary condition for the pressure equation,
     *  either pressure (dirichlet) or flux (neumann),
     *  and for the transport equation,
     *  either saturation (dirichlet) or flux (neumann).
     *
     *  \param bcTypes Includes the types of boundary conditions
     *  \param globalPos The position of the center of the finite volume
     *
     *  Alternatively, the function boundaryTypes(PrimaryVariables &values, const Intersection&
     *  intersection) could be defined, where intersection is the boundary intersection.
     */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const /*@\label{tutorial-sequential:bctype}@*/
    {
            if (globalPos[0] < this->bBoxMin()[0] + eps_)
            {
               bcTypes.setAllDirichlet(); // alternative if the same BC is used for both types of equations
            }
            // all other boundaries
            else
            {
               bcTypes.setAllNeumann(); // alternative if the same BC is used for both types of equations
            }
    }
    //! Value for dirichlet boundary condition at position globalPos.
    /*! In case of a dirichlet BC for the pressure equation the pressure \f$ [Pa] \f$, and for
     *  the transport equation the saturation [-] have to be defined on boundaries.
     *
     *  \param values Values of primary variables at the boundary
     *  \param intersection The boundary intersection
     *
     *  Alternatively, the function dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    void dirichlet(PrimaryVariables &values, const Intersection& intersection) const /*@\label{tutorial-sequential:dirichlet}@*/
    {
        values[pwIdx] = 2e5;
        values[swIdx] = 0.0;
    }
    //! Value for neumann boundary condition \f$ [\frac{kg}{m^3 \cdot s}] \f$ at position globalPos.
    /*! In case of a neumann boundary condition, the flux of matter
     *  is returned as a vector.
     *
     *  \param values Boundary flux values for the different phases
     *  \param globalPos The position of the center of the finite volume
     *
     *  Alternatively, the function neumann(PrimaryVariables &values, const Intersection& intersection) could be defined,
     *  where intersection is the boundary intersection.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const /*@\label{tutorial-sequential:neumann}@*/
    {
        values = 0;
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            values[wPhaseIdx] = 2e-4;
        }
    }
    //! Initial condition at position globalPos.
    /*! Only initial values for saturation have to be given!
     *
     *  \param values Values of primary variables
     *  \param element The finite volume element
     *
     *  Alternatively, the function initialAtPos(PrimaryVariables &values, const GlobalPosition& globalPos)
     *  could be defined, where globalPos is the vector including the global coordinates of the finite volume.
     */
    void initial(PrimaryVariables &values,
            const Element &element) const /*@\label{tutorial-sequential:initial}@*/
    {
        values[swIdx] = 1.0;
    }

private:
    const Scalar eps_;
};
} //end namespace

#endif
