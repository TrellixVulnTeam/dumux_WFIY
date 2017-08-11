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
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_1P_HETEROGENEOUS_MATRIX_PROBLEM_HH
#define DUMUX_1P_HETEROGENEOUS_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "heterogeneousmatrixspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePMpfaMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePMatrixProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePCCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, OnePMatrixProblem));

SET_TYPE_PROP(OnePMatrixProblem, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), Unit<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Set the problem property
SET_TYPE_PROP(OnePMatrixProblem, Problem, OnePMpfaMatrixProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(OnePMatrixProblem, SpatialParams, OnePMatrixSpatialParams<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(OnePMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(OnePMatrixProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(OnePCCMpfaMatrixProblem, SolutionDependentAdvection, false);

SET_BOOL_PROP(OnePCCMpfaMatrixProblem, EnableGlobalFVGeometryCache, true);

// change mpfa method
//SET_PROP(OnePMatrixProblem, MpfaMethod) { static const MpfaMethods value = MpfaMethods::lMethod; };
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model
 */

template <class TypeTag>
class OnePMpfaMatrixProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

    enum
    {
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    OnePMpfaMatrixProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_matrix";
        eps_ = 1e-6;
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }

    /*!
     * \brief Return the sources within the domain.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[0] < this->bBoxMin()[0] + eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Specifies if a given intersection is on an interior boundary
     */
    bool isInteriorBoundary(const Element& element, const Intersection& is) const
    { return couplingManager().isInteriorBoundary(element, is); }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < this->bBoxMin()[0] + eps_)
            return PrimaryVariables(0.0);
        else
            return PrimaryVariables(1.0);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1; // initial condition for the pressure
        return priVars;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    CouplingManager& couplingManager()
    { return *couplingManager_; }

private:
    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
