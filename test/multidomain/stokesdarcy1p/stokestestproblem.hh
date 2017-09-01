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
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dumux/implicit/staggered/properties.hh>
#include <dumux/freeflow/staggered/model.hh>
#include <dumux/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the Stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(StaggeredModel, NavierStokes));

// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<3, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 3> >);
#else
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);
#endif

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, StokesTestProblem<TypeTag>);

SET_PROP(StokesTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::Constant<TypeTag, Scalar> >;
};

SET_BOOL_PROP(StokesTestProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(StokesTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(StokesTestProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(StokesTestProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(StokesTestProblem, EnableInertiaTerms, false);

// Set the grid parameter group
SET_STRING_PROP(StokesTestProblem, GridParameterGroup, "StokesGrid");

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);

}

/*!
 * \ingroup StaggeredStokesModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes flow problem with water (H2O) flowing from the left to the right.
 *
 * The domain is sized 2m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow boundary conditions.
 *
 * This problem uses the \ref StaggeredStokesModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes</tt>
 */
template <class TypeTag>
class StokesTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;

    using SubControlVolumeFace =  typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    using BoundaryValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);
    using InitialValues = typename GET_PROP_TYPE(TypeTag, BoundaryValues);

public:
    StokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
{
        eps_ = 1e-6;
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        vIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Velocity);

        bBoxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[0];
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[0];
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions0)[1];
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, StokesGrid, Positions1)[1];
}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     *
     */
    const std::string name() const
    {
        return name_+"_stokes";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10.0; // -> 10C
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(momentumBalanceIdx);
        }
        else
            values.setOutflow(massBalanceIdx);

        return values;
    }

    /*!
     * \brief Return Dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryValues values;

        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = 0.0;

        if (onLeftBoundary_(globalPos))
        {
            values[velocityXIdx] = vIn_;
        }
        return values;
    }

    /*!
     * \brief Return Neumann boundary values at a given position
     *
     * \param globalPos The global position
     */
    BoundaryValues neumannAtPos(const GlobalPosition &globalPos) const
    {
        return BoundaryValues(0.0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);
        values[pressureIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
     * \brief Return the coupling boundary
     *
     * \param globalPos The global position
     */
    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return onLowerBoundary_(globalPos); }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    Scalar eps_;
    Scalar vIn_;
    std::string name_;

    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace

#endif
