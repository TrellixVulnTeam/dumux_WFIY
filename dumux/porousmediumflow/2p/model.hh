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
 * \ingroup TwoPModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPFormulation::pwsn</tt> or <tt>TwoPFormulation::pnsw</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/fv.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/vtkoutputfields.hh>

#include "formulation.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "saturationreconstruction.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPModel
 * \brief Specifies a number properties of two-phase models.
 */
template<TwoPFormulation formulation>
struct TwoPModelTraits
{
    using Indices = TwoPIndices;

    static constexpr TwoPFormulation priVarFormulation()
    { return formulation; }

    static constexpr int numEq() { return 2; }
    static constexpr int numPhases() { return 2; }
    static constexpr int numComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }

    static std::string primaryVariableName(int pvIdx)
    {
        if (priVarFormulation() == TwoPFormulation::p0s1)
            return pvIdx == 0 ? "pw" : "Sn";
        else
            return pvIdx == 0 ? "pn" : "Sw";
    }
};

/*!
 * \ingroup TwoPModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam SSY The solid system type
 * \tparam SST The solid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 * \tparam SR The class used for reconstruction of
 *            non-wetting phase saturations in scvs
 */
template<class PV, class FSY, class FST,class SSY, class SST, class PT, class MT, class SR>
struct TwoPVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using SaturationReconstruction = SR;
};

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal two-phase model
NEW_TYPE_TAG(TwoP, INHERITS_FROM(PorousMediumFlow));
//! The type tag for the non-isothermal two-phase model
NEW_TYPE_TAG(TwoPNI, INHERITS_FROM(TwoP));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////
 //!< Set the default formulation to pwsn
SET_PROP(TwoP, Formulation)
{ static constexpr auto value = TwoPFormulation::p0s1; };

SET_TYPE_PROP(TwoP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);         //!< Use the immiscible local residual operator for the 2p model

//! The model traits class
SET_TYPE_PROP(TwoP, ModelTraits, TwoPModelTraits<GET_PROP_VALUE(TypeTag, Formulation)>);

//! Set the vtk output fields specific to the twop model
SET_TYPE_PROP(TwoP, VtkOutputFields, TwoPVtkOutputFields);

//! Set the volume variables property
SET_PROP(TwoP, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    static constexpr auto DM = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod;
    static constexpr bool enableIS = GET_PROP_VALUE(TypeTag, EnableBoxInterfaceSolver);
    // class used for scv-wise reconstruction of non-wetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;

    using Traits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;
public:
    using type = TwoPVolumeVariables<Traits>;
};

//! The two-phase model uses the immiscible fluid state
SET_PROP(TwoP, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

////////////////////////////////////////////////////////
// properties for the non-isothermal two-phase model
////////////////////////////////////////////////////////

//! The non-isothermal model traits class
SET_TYPE_PROP(TwoPNI, ModelTraits, PorousMediumFlowNIModelTraits<TwoPModelTraits<GET_PROP_VALUE(TypeTag, Formulation)>>);

//! Set the vtk output fields specific to the non-isothermal twop model
SET_TYPE_PROP(TwoPNI, VtkOutputFields, EnergyVtkOutputFields<TwoPVtkOutputFields>);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
