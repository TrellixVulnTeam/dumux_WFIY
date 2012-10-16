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
 * \brief A test problem for the one-phase box model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_1PTEST_PROBLEM_HH
#define DUMUX_1PTEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dumux/boxmodels/1p/1pmodel.hh>
#include <dumux/boxmodels/common/porousmediaboxproblem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1ptestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePTestProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTestProblem, INHERITS_FROM(BoxOneP));

SET_PROP(OnePTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_PROP(OnePTestProblem, Grid)
{
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
  //typedef Dune::UGGrid<2> type;
  //typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> type;
};

// Set the problem property
SET_PROP(OnePTestProblem, Problem)
{ typedef Dumux::OnePTestProblem<TypeTag> type; };

// Set the spatial parameters
SET_PROP(OnePTestProblem, SpatialParams)
{ typedef Dumux::OnePTestSpatialParams<TypeTag> type; };

// Linear solver settings
SET_TYPE_PROP(OnePTestProblem, LinearSolver, Dumux::BoxBiCGStabILU0Solver<TypeTag> );
SET_INT_PROP(OnePTestProblem, LinearSolverVerbosity, 0);
SET_INT_PROP(OnePTestProblem, PreconditionerIterations, 1);
SET_SCALAR_PROP(OnePTestProblem, PreconditionerRelaxation, 1.0);

// Enable gravity
SET_BOOL_PROP(OnePTestProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxTestProblems
 * \brief  Test problem for the one-phase box model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p -parameterFile test_1p.input</tt>
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::SGrid<2,2> type;</tt> to
 * <tt>typedef Dune::SGrid<3,3> type;</tt> in the problem file
 * and use <tt>1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class OnePTestProblem : public PorousMediaBoxProblem<TypeTag>
{
    typedef PorousMediaBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
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
    { return "1ptest"; }

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


    void sourceAtPos(PrimaryVariables &priVars,
                const GlobalPosition &globalPos) const
    {
        priVars = 0;
    }
    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specify which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        double eps = 1.0e-3;
        if (globalPos[dim-1] < eps || globalPos[dim-1] > this->bboxMax()[dim-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &priVars, const Vertex &vertex) const
    {
        double eps = 1.0e-3;
        const GlobalPosition globalPos = vertex.geometry().center();

        if (globalPos[dim-1] < eps) {
            priVars[pressureIdx] = 2.0e+5;
        }
        else if (globalPos[dim-1] > this->bboxMax()[dim-1] - eps) {
            priVars[pressureIdx] = 1.0e+5;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    using ParentType::neumann;
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        //  const GlobalPosition &globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        priVars[conti0EqIdx] = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        //const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        priVars[pressureIdx] = 1.0e+5;// + 9.81*1.23*(20-globalPos[dim-1]);
    }

    // \}
};
} //end namespace

#endif
