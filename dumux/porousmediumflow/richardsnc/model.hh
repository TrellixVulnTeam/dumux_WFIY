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
 * \ingroup RichardsNCModel
 * \brief Base class for all models which use the Richards,
 *        n-component fully implicit model.
 *
 * In the unsaturated zone, Richards' equation
 *\f{eqnarray*}
 && \frac{\partial (\sum_w \varrho_w X_w^\kappa \phi S_w )}
 {\partial t}
 - \sum_w  \text{div} \left\{ \varrho_w X_w^\kappa
 \frac{k_{rw}}{\mu_w} \mbox{\bf K}
 (\text{grad}\, p_w - \varrho_{w}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_w \text{div} \left\{{\bf D_{w, pm}^\kappa} \varrho_{w} \text{grad}\, X^\kappa_{w} \right\}
 - \sum_w q_w^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 w \in \{w, g\}
 \f}
 * is frequently used to
 * approximate the water distribution above the groundwater level.
 *
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the non-wetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitly large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus choses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_n - p_w)
 \f]
 * holds, where \f$p_n\f$ is a given reference pressure. Nota bene,
 * that the last step is assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */

#ifndef DUMUX_RICHARDSNC_MODEL_HH
#define DUMUX_RICHARDSNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/vtkoutputfields.hh>

#include <dumux/porousmediumflow/richards/model.hh>

#include "volumevariables.hh"
#include "indices.hh"
#include "vtkoutputfields.hh"

namespace Dumux {

/*!
 * \ingroup RichardsNCModel
 * \brief Specifies a number properties of the Richards n-components model.
 *
 * \tparam nComp the number of components to be considered.
 * \tparam useMol whether to use mass or mole balances
 */
template<int nComp, bool useMol, int repCompEqIdx = nComp>
struct RichardsNCModelTraits
{
    using Indices = RichardsNCIndices;

    static constexpr int numEq() { return nComp; }
    static constexpr int numPhases() { return 1; }
    static constexpr int numComponents() { return nComp; }
    static constexpr int replaceCompEqIdx() { return repCompEqIdx; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }

    static constexpr bool useMoles() { return useMol; }

    template <class FluidSystem, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx, int state = 0)
    {
        const std::string xString = useMoles() ? "x" : "X";
        if (pvIdx == 0)
            return "p_" + FluidSystem::phaseName(0);
        else
            return xString + "^" + FluidSystem::componentName(pvIdx)
                   + "_" + FluidSystem::phaseName(0);
    }
};

namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
NEW_TYPE_TAG(RichardsNC, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(RichardsNCNI, INHERITS_FROM(RichardsNC));
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the model traits class
SET_PROP(RichardsNC, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = RichardsNCModelTraits<FluidSystem::numComponents, GET_PROP_VALUE(TypeTag, UseMoles), GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx)>;
};

//! Define that per default mole fractions are used in the balance equations
SET_BOOL_PROP(RichardsNC, UseMoles, true);

//! Use the dedicated local residual
SET_TYPE_PROP(RichardsNC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! We set the replaceCompIdx to 0, i.e. the first equation is substituted with
//! the total mass balance, i.e. the phase balance
SET_INT_PROP(RichardsNC, ReplaceCompEqIdx, 0);

//! Set the volume variables property
SET_PROP(RichardsNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SSY = typename GET_PROP_TYPE(TypeTag, SolidSystem);
    using SST = typename GET_PROP_TYPE(TypeTag, SolidState);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using PT = typename GET_PROP_TYPE(TypeTag, SpatialParams)::PermeabilityType;

    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

    using Traits = RichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;
public:
    using type = RichardsNCVolumeVariables<Traits>;
};

//! The default richardsnc model computes no diffusion in the air phase
//! Turning this on leads to the extended Richards equation (see e.g. Vanderborght et al. 2017)
SET_BOOL_PROP(RichardsNC, EnableWaterDiffusionInAir, false);

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the liquid phase fluid system with simple H2O.
 */
SET_PROP(RichardsNC, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::SimpleH2O<Scalar>, Components::Constant<1, Scalar>>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(RichardsNC, FluidState)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the vtk output fields specific to this model
SET_TYPE_PROP(RichardsNC, VtkOutputFields, RichardsNCVtkOutputFields);

//! The model after Millington (1961) is used for the effective diffusivity
SET_TYPE_PROP(RichardsNC, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! average is used as default model to compute the effective thermal heat conductivity
SET_TYPE_PROP(RichardsNCNI, ThermalConductivityModel, ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Richards n-components model
//////////////////////////////////////////////////////////////////

//! set non-isothermal model traits
SET_PROP(RichardsNCNI, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using IsothermalTraits = RichardsNCModelTraits<FluidSystem::numComponents, GET_PROP_VALUE(TypeTag, UseMoles), GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx)>;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
