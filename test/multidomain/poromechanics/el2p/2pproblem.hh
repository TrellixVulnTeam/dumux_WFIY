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
 * \ingroup TwoPTests
 * \ingroup PoroElastic
 * \brief Definition of the spatial parameters for the two-phase flow
 *        sub-problem in the coupled poro-mechanical elp problem.
 */
#ifndef DUMUX_2P_SUB_PROBLEM_HH
#define DUMUX_2P_SUB_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>

#include "2pspatialparams.hh"
#include "el2pco2tables.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class TwoPSubProblem;

namespace Properties {

NEW_TYPE_TAG(TwoPSubTypeTag, INHERITS_FROM(CCTpfaModel, TwoP));

// Set the fluid system for TwoPSubProblem
SET_PROP(TwoPSubTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::BrineCO2<Scalar, El2P::CO2Tables>;
};

// Set the grid type
SET_TYPE_PROP(TwoPSubTypeTag, Grid, Dune::YaspGrid<3>);
// Set the problem property
SET_TYPE_PROP(TwoPSubTypeTag, Problem, TwoPSubProblem<TypeTag> );
// Set the spatial parameters
SET_PROP(TwoPSubTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using type = TwoPSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
};
} // end namespace Properties

/*!
 * \ingroup MultiDomain
 * \ingroup TwoPTests
 * \ingroup PoroElastic
 *
 * \brief The two-phase sub problem in the el2p coupled problem.
 */
template <class TypeTag>
class TwoPSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy pressure index for convenience
    enum {
          pressureIdx = GET_PROP_TYPE(TypeTag, ModelTraits)::Indices::pressureIdx,
          saturationNIdx = GET_PROP_TYPE(TypeTag, ModelTraits)::Indices::saturationIdx,
          waterPhaseIdx = FluidSystem::phase0Idx,
          gasPhaseIdx = FluidSystem::phase1Idx,
          dimWorld = GridView::dimensionworld
    };

    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    TwoPSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
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
    {
      PrimaryVariables values;

      values[pressureIdx] = 1.5e7;
      values[saturationNIdx] = 0.0;
      return values;
    }

    //! Evaluate source terms
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        static const Scalar sourceG = getParam<Scalar>("Problem.InjectionRateGas");
        static const Scalar sourceW = getParam<Scalar>("Problem.InjectionRateWater");
        if(globalPos[0] > 250 + eps_ && globalPos[0] < 750 - eps_
           && globalPos[1] > 250 + eps_ && globalPos[1] < 750 - eps_
           && globalPos[dimWorld-1] > 250 + eps_ && globalPos[dimWorld-1] < 750 - eps_)
        {
            values[gasPhaseIdx] = sourceG;
            values[waterPhaseIdx] = sourceW;
        }
        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[dimWorld-1] < eps_)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
};

} //end namespace Dumux

#endif
