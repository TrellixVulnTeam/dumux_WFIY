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
 * \ingroup FreeflowNCModel
 *
 * \brief A single-phase, multi-component low-Re k-epsilon model
 *
 * \copydoc Dumux::FreeflowNCModel
 */

#ifndef DUMUX_LOWREKEPSILON_NC_MODEL_HH
#define DUMUX_LOWREKEPSILON_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>

#include "vtkoutputfields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component low-Re k-epsilon model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal low-Re k-epsilon model
NEW_TYPE_TAG(LowReKEpsilonNC, INHERITS_FROM(NavierStokesNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the low-Reynolds k-epsilon multi-component model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam useM Use molar or mass balances
 * \tparam replaceCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 */
template<int dimension, int nComp, bool useMoles, int replaceCompEqIdx>
struct LowReKEpsilonNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp+2; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! the indices
    using Indices = LowReKEpsilonIndices<dimension, nComp>;

    //! return the names of the primary variables in cells
    template <class FluidSystem>
    static std::string primaryVariableNameCell(int pvIdx, int state = 0)
    {
        using ParentType = NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>;
        if (pvIdx < nComp)
            return ParentType::template primaryVariableNameCell<FluidSystem>(pvIdx, state);
        else if (pvIdx == nComp)
            return "k";
        else
            return "epsilon";
    }
};

//!< states some specifics of the isothermal multi-component low-Reynolds k-epsilon model
SET_PROP(LowReKEpsilonNC, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = LowReKEpsilonNCModelTraits<dimension, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilonNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(LowReKEpsilonNC, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = LowReKEpsilonResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(LowReKEpsilonNC, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = LowReKEpsilonFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(LowReKEpsilonNC, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SinglePhaseVtkOutputFields = LowReKEpsilonVtkOutputFields<FVGridGeometry>;
public:
    using type = FreeflowNCVtkOutputFields<SinglePhaseVtkOutputFields, ModelTraits, FVGridGeometry, FluidSystem>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component low-Re k-epsilon model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal low-Re k-epsilon models
NEW_TYPE_TAG(LowReKEpsilonNCNI, INHERITS_FROM(NavierStokesNCNI));

//! The model traits of the non-isothermal model
SET_PROP(LowReKEpsilonNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    using IsothermalModelTraits = LowReKEpsilonNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilonNCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(LowReKEpsilonNCNI, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = LowReKEpsilonResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(LowReKEpsilonNCNI, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = LowReKEpsilonFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(LowReKEpsilonNCNI, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BaseVtkOutputFields = LowReKEpsilonVtkOutputFields<FVGridGeometry>;
    using NonIsothermalFields = FreeflowNonIsothermalVtkOutputFields<BaseVtkOutputFields, ModelTraits>;
public:
    using type = FreeflowNCVtkOutputFields<NonIsothermalFields, ModelTraits, FVGridGeometry, FluidSystem>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
