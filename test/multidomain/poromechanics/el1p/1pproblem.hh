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
 * \ingroup MultiDomain
 * \ingroup OnePTests
 * \ingroup PoroElastic
 * \brief Definition of the spatial parameters for the single-phase flow
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
#ifndef DUMUX_1P_SUB_PROBLEM_HH
#define DUMUX_1P_SUB_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "1pspatialparams.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class OnePSubProblem;

namespace Properties {

NEW_TYPE_TAG(OnePSubTypeTag, INHERITS_FROM(CCTpfaModel, OneP));

// The fluid phase consists of one constant component
SET_TYPE_PROP(OnePSubTypeTag,
              FluidSystem,
              Dumux::FluidSystems::OnePLiquid< typename GET_PROP_TYPE(TypeTag, Scalar),
                                               Dumux::Components::Constant<0, typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Set the grid type
SET_TYPE_PROP(OnePSubTypeTag, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(OnePSubTypeTag, Problem, OnePSubProblem<TypeTag> );
// Set the spatial parameters
SET_PROP(OnePSubTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using type = OnePSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
};
} // end namespace Properties

/*!
 * \ingroup MultiDomain
 * \ingroup OnePTests
 * \ingroup PoroElastic
 *
 * \brief The single-phase sub problem in the el1p coupled problem.
 */
template <class TypeTag>
class OnePSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy pressure index for convenience
    enum { pressureIdx = GET_PROP_TYPE(TypeTag, ModelTraits)::Indices::pressureIdx };

    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    OnePSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    {}

    //! Return the temperature within the domain in [K].
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    //! Evaluate the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0e5); }

    //! Evaluate source terms
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar source = getParam<Scalar>("Problem.InjectionRate");
        if (globalPos[0] > 0.4 && globalPos[0] < 0.6 && globalPos[1] < 0.6 && globalPos[1] > 0.4)
            return NumEqVector(source);
        return NumEqVector(0.0);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
};

} //end namespace Dumux

#endif
