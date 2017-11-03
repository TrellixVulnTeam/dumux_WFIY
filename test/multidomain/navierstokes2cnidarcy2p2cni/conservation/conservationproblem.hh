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
#ifndef DUMUX_CONSERVATION_PROBLEM_HH
#define DUMUX_CONSERVATION_PROBLEM_HH

#include "conservationdarcysubproblem.hh"
#include "conservationstokessubproblem.hh"
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/multidomain/problem.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/map.hh>

#include <dumux/multidomain/staggered-ccfv/properties.hh>
#include <dumux/multidomain/staggered-ccfv/model.hh>
#include <dumux/multidomain/staggered-ccfv/assembler.hh>
#include <dumux/multidomain/staggered-ccfv/newtoncontroller.hh>
#include <dumux/multidomain/staggered-ccfv/subproblemlocaljacobian.hh>

//#include <appl/staggeredgrid/multidomain/navierstokes2ctdarcy2p2ct/problem.hh> TODO - ok?

namespace Dumux
{
template <class TypeTag>
class ConservationProblem;

namespace Properties
{
NEW_TYPE_TAG(ConservationProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(ConservationProblem, Problem, Dumux::ConservationProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(ConservationProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(ConservationProblem, DarcyProblemTypeTag, TTAG(DarcySubProblem));
SET_TYPE_PROP(ConservationProblem, StokesProblemTypeTag, TTAG(StokesSubProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(DarcySubProblem, GlobalProblemTypeTag, TTAG(ConservationProblem));
SET_TYPE_PROP(StokesSubProblem, GlobalProblemTypeTag, TTAG(ConservationProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(DarcySubProblem, ParameterTree) : GET_PROP(TTAG(ConservationProblem), ParameterTree) {};
SET_PROP(StokesSubProblem, ParameterTree) : GET_PROP(TTAG(ConservationProblem), ParameterTree) {};

// Set the fluid system to use complex relations (last argument)
SET_TYPE_PROP(ConservationProblem, FluidSystem, FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);
SET_TYPE_PROP(DarcySubProblem, FluidSystem, typename GET_PROP_TYPE(TTAG(ConservationProblem), FluidSystem));
SET_TYPE_PROP(StokesSubProblem, FluidSystem, typename GET_PROP_TYPE(TTAG(ConservationProblem), FluidSystem));

//// Define whether mole(true) or mass (false) fractions are used
//SET_BOOL_PROP(ConservationProblem, UseMoles, false);
//SET_BOOL_PROP(DarcySubProblem, UseMoles, GET_PROP_VALUE(TTAG(ConservationProblem), UseMoles));
//SET_BOOL_PROP(StokesSubProblem, UseMoles, GET_PROP_VALUE(TTAG(ConservationProblem), UseMoles));

NEW_PROP_TAG(DarcyToStokesMapValue); // TODO: make specialized map value class
SET_TYPE_PROP(ConservationProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);
}

template <class TypeTag>
class ConservationProblem : public MultiDomainProblem<TypeTag>
{
    using ParentType = MultiDomainProblem<TypeTag>;

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    // obtain types from the sub problem type tags
    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

public:
    /*!
     * \brief Base class for the multi domain problem
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The GridView
     */
    ConservationProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/343.15, /*numTemp=*/140,
                          /*pMin=*/5e4, /*pMax=*/1.5e5, /*numP=*/100);
    }

    bool shouldWriteOutput() const //define output
    {
        return true;
    }

private:
};

} // end namespace Dumux

#endif // DUMUX_CONSERVATION_PROBLEM_HH
