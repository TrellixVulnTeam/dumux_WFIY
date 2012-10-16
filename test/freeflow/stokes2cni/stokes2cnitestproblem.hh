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
 * @file
 * @brief  Definition of a simple Stokes problem
 */
#ifndef DUMUX_STOKES2CNITESTPROBLEM_HH
#define DUMUX_STOKES2CNITESTPROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/freeflow/stokes2cni/stokes2cnimodel.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cniTestProblem;

//////////
// Specify the properties for the stokes2cni problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cniTestProblem, INHERITS_FROM(BoxStokes2cni));

// Set the grid type
SET_TYPE_PROP(Stokes2cniTestProblem, Grid, Dune::SGrid<2,2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cniTestProblem, Problem, Stokes2cniTestProblem<TypeTag>);

//! Select the fluid system
SET_PROP(BoxStokes2cni, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2OAir<Scalar> type;
};

//! Scalar is set to type long double for higher accuracy
//SET_TYPE_PROP(BoxStokes, Scalar, long double);

//! a stabilization factor. Set to zero for no stabilization
SET_SCALAR_PROP(BoxStokes2cni, StokesStabilizationAlpha, -1.0);

// Enable gravity
SET_BOOL_PROP(Stokes2cniTestProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxTestProblems
 * \brief Stokes2cni problem with air flowing
 *        from the bottom to the top, blowing away a warm and dry square.
 *
 * The domain is sized 1m times 1m. An air flow from the bottom boundary blows a warm and dry square,
 * which is initially in the center of the domain out of the upper boundary.
 * The boundary conditions for the momentum balances
 * are all set to Dirichlet. The mass balance has outflow boundary conditions, which are
 * replaced in the localresidual by the sum of the two momentum balances equations in case of Dirichlet bcs for the momentum balance.
 * On the upper boundary a Dirichlet condition is set for the mass balance to fix the pressure.
 * Gravity is on in this example.
 *
 * This problem uses the \ref BoxStokes2cniModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes2cni  -parameterFile ./test_stokes2cni.input</tt>
 */
template <class TypeTag>
class Stokes2cniTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cniIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };
    enum { // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx, //!< Index of the transport equation
        energyEqIdx =    Indices::energyEqIdx     //!< Index of the energy equation
    };
    enum { // indices for primary variables
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    Stokes2cniTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

        // initialize the tables of the fluid system
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
    { return "stokes2cni"; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc BoxProblem::boundaryTypes()
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        values.setAllDirichlet();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        if (onUpperBoundary_(globalPos) &&
                !onLeftBoundary_(globalPos) && !onRightBoundary_(globalPos))
            values.setAllOutflow();

        // set pressure at the upper boundary
        if (onUpperBoundary_(globalPos) &&
                !onLeftBoundary_(globalPos) && !onRightBoundary_(globalPos))
            values.setDirichlet(massBalanceIdx);
    }

    //! \copydoc BoxProblem::dirichlet()
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    //! \copydoc BoxProblem::neumann()
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0.0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! \copydoc BoxProblem::source()
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
    }

    //! \copydoc BoxProblem::initial()
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }
   // \}

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        const Scalar v1 = 0.5;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = v1*(globalPos[0] - this->bboxMin()[0])*(this->bboxMax()[0] - globalPos[0])
                                   / (0.25*(this->bboxMax()[0] - this->bboxMin()[0])*(this->bboxMax()[0] - this->bboxMin()[0]));
        values[pressureIdx] = 1e5 + 1.189*this->gravity()[1]*globalPos[1];
        values[massOrMoleFracIdx] = 1e-4;
        values[temperatureIdx] = 283.15;
        if(globalPos[0]<0.75 && globalPos[0]>0.25 &&
                globalPos[1]<0.75 && globalPos[1]>0.25)
        {
            values[massOrMoleFracIdx] = 0.9e-4;
            values[temperatureIdx] = 284.15;
        }
    }
    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    Scalar eps_;
};
} //end namespace

#endif
