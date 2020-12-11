// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief The properties for the one-phase pore-network model with non-creeping flow.
 */
#ifndef DUMUX_PNM1P_NONCREEPING_PROPERTIES_HH
#define DUMUX_PNM1P_NONCREEPING_PROPERTIES_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/properties.hh>

// Pore network model
#include <dumux/porenetwork/1p/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"
#include "problem.hh"

//////////
// Specify the properties
//////////
namespace Dumux::Properties {

namespace TTag {
#if ISOTHERMAL
struct PNMOnePNonCreepingProblem { using InheritsFrom = std::tuple<PNMOneP>; };
#else
struct PNMOnePNonCreepingProblem { using InheritsFrom = std::tuple<PNMOnePNI>; };
#endif
}

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePNonCreepingProblem>
{ using type = PNMOnePNonCreepingProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePNonCreepingProblem>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> > ;
};

// the grid
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePNonCreepingProblem>
{ using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMOnePNonCreepingProblem>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::NonCreepingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMOnePNonCreepingProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = Dumux::PoreNetwork::TransmissibilityPatzekSilin<Scalar, false/*considerPoreBodyResistance*/>;
public:
    using type = Dumux::PoreNetwork::NonCreepingFlow<Scalar, TransmissibilityLaw>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMOnePNonCreepingProblem> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

} //end namespace Dumux::Properties

#endif
