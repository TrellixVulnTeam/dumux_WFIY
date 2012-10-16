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
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 */
#ifndef DUMUX_KUEVETTE3P3CNIPROBLEM_HH
#define DUMUX_KUEVETTE3P3CNIPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>

#include <dumux/boxmodels/3p3cni/3p3cnimodel.hh>
#include <dumux/boxmodels/common/porousmediaboxproblem.hh>

#include "kuevettespatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class KuevetteProblem;

namespace Properties
{
NEW_TYPE_TAG(KuevetteProblem, INHERITS_FROM(BoxThreePThreeCNI, KuevetteSpatialParams));

// Set the grid type
SET_PROP(KuevetteProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};


// Set the problem property
SET_PROP(KuevetteProblem, Problem)
{
    typedef Dumux::KuevetteProblem<TypeTag> type;
};

// Set the fluid system
SET_TYPE_PROP(KuevetteProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(KuevetteProblem, ProblemEnableGravity, true);

// Use central differences (backward -1, forward +1)
SET_INT_PROP(KuevetteProblem, ImplicitNumericDifferenceMethod, 0);

// Write newton convergence
//SET_BOOL_PROP(KuevetteProblem, NewtonWriteConvergence, true);

// Set the maximum time step
SET_SCALAR_PROP(KuevetteProblem, TimeManagerMaxTimeStepSize, 60.);

// set newton relative tolerance
SET_SCALAR_PROP(KuevetteProblem, NewtonRelTolerance, 1e-6);
}


/*!
 * \ingroup ThreePThreeCNIBoxModel
 * \ingroup BoxTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 *
 * The domain is a quasi-two-dimensional container (kuevette). Its dimensions
 * are 1.5 m x 0.74 m. The top and bottom boundaries are closed, the right
 * boundary is a Dirichlet boundary allowing fluids to escape. From the left,
 * an injection of a hot water-air mixture is applied (Neumann boundary condition
 * for the mass components and the enthalpy), aimed at remediating an initial 
 * NAPL (Non-Aquoeus Phase Liquid) contamination in the heterogeneous domain.
 * The contamination is initially placed partly into the coarse sand
 * and partly into a fine sand lense.
 *
 * This simulation can be varied through assigning different boundary conditions
 * at the left boundary as described in Class (2001):
 * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in
 * NAPL-kontaminierten por"osen Medien, Dissertation, Eigenverlag des Instituts
 * f"ur Wasserbau
 *
 * This problem uses the \ref ThreePThreeCNIModel
 *
 * To see the basic effect and the differences to scenarios with pure steam or
 * pure air injection, it is sufficient to simulated for about 2-3 hours (10000 s).
 * Complete remediation of the domain requires much longer (about 10 days simulated time).
 * To adjust the simulation time it is necessary to edit the file test_3p3cni.input
 *
 * To run the simulation execute:
 *
 * <tt>./test_3p3cni -parameterFile test_3p3cni.input</tt>
 *  */
template <class TypeTag >
class KuevetteProblem : public PorousMediaBoxProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    typedef PorousMediaBoxProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        // Phase State
        threePhases = Indices::threePhases,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    KuevetteProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView), eps_(1e-6)
    {
        FluidSystem::init();
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
    const char *name() const
    { return "kuevette3p3cni"; }

    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0;
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
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if(globalPos[0] > 1.5 - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex for which the boundary type is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeomtry The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeomtry,
                 const Intersection &is,
                 const int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values = 0;

        // negative values for injection
        if (globalPos[0] < eps_)
        {
            values[Indices::contiWEqIdx] = -0.1435; // 0.3435 [mol/(s m)] in total
            values[Indices::contiAEqIdx] = -0.2;
            values[Indices::contiCEqIdx] =  0.0;
            values[Indices::energyEqIdx] = -6929.;
        }
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
     * \param element The finite element
     * \param fvGeomtry The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeomtry,
                 const int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);

    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             const int &globalIdx,
                             const GlobalPosition &globalPos) const
    {
        if((globalPos[0] >= 0.20) && (globalPos[0] <= 0.80) && (globalPos[1] >= 0.4) && (globalPos[1] <= 0.65))
        {
            return threePhases;
        }
        else return wgPhaseOnly;
    }


private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 1e5 ;
        values[switch1Idx] = 0.12;
        values[switch2Idx] = 1.e-6;
        values[temperatureIdx] = 293.0;

        if((globalPos[0] >= 0.20) && (globalPos[0] <= 0.80) && (globalPos[1] >= 0.4) && (globalPos[1] <= 0.65))
        {
            values[switch2Idx] = 0.07;
        }
    }

    const Scalar eps_;
};
} //end namespace

#endif
