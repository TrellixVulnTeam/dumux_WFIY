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
 * \brief Helper class to be used to obtain lambda functions for the tensors
 *        involved in the laws that describe the different kind of fluxes that
 *        occur in DuMuX models (i.e. advective, diffusive and heat conduction fluxes).
 *        The local systems appearing in Mpfa methods have to be solved subject to
 *        the different tensors. This class returns lambda expressions to be used in the
 *        local systems. The specialization for other discretization methods allows
 *        compatibility with the TPFA scheme, which could be used for one or more of the tensors.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_TENSOR_LAMBDA_FACTORY_HH
#define DUMUX_DISCRETIZATION_MPFA_TENSOR_LAMBDA_FACTORY_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/tensorlambdafactory.hh>

namespace Dumux
{

//! forward declaration of properties
namespace Properties
{
NEW_PROP_TAG(EffectiveDiffusivityModel);
NEW_PROP_TAG(ThermalConductivityModel);
};

/*!
 * \ingroup MpfaModel
 * \brief Helper class to be used to obtain lambda functions for the tensors
 *        involved in the laws that describe the different kind of fluxes that
 *        occur in DuMuX models (i.e. advective, diffusive and heat conduction fluxes).
 *        The local systems appearing in Mpfa methods have to be solved subject to
 *        the different tensors. This class returns lambda expressions to be used in the
 *        local systems. The specialization for discretization methods other than mpfa allows
 *        compatibility with the TPFA scheme, which could be used for one or more of the tensors.
 *        The interfaces of the lambdas are chosen such that all involved tensors can be extracted
 *        with the given arguments.
 */
template<class TypeTag, DiscretizationMethods Method>
class TensorLambdaFactory
{
public:

    //! We return zero scalars here in the functions below.
    //! We have to return something as the local systems expect a type
    //! to perform actions on. We return 0.0 as this call should never happen
    //! for a tensor which is not treated by an mpfa method anyway.

    //! lambda for the law describing the advective term
    static auto getAdvectionLambda()
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return 0.0; };
    }

    //! lambda for the diffusion coefficient/tensor
    static auto getDiffusionLambda(unsigned int phaseIdx, unsigned int compIdx)
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return 0.0; };
    }

    //! lambda for the fourier coefficient
    static auto getHeatConductionLambda()
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return 0.0; };
    }
};

//! Specialization for mpfa schemes
template<class TypeTag>
class TensorLambdaFactory<TypeTag, DiscretizationMethods::CCMpfa>
{
public:

    //! return advection tensor/scalar lambda here
    static auto getAdvectionLambda()
    {
        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return volVars.permeability(); };
    }

    //! return diffusion tensor/scalar lambda here (single-phasic case)
    template<class T = TypeTag>
    static auto getDiffusionLambda(unsigned int phaseIdx,
                                   unsigned int compIdx,
                                   typename std::enable_if<(GET_PROP_VALUE(T, NumPhases) == 1), void>::type* dummy = nullptr)
    {
        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
        return [phaseIdx, compIdx] (const auto& problem,
                                    const auto& element,
                                    const auto& volVars,
                                    const auto& fvGeometry,
                                    const auto& scv)
               { return EffDiffModel::effectiveDiffusivity(volVars.porosity(),
                                                           volVars.saturation(phaseIdx),
                                                           volVars.diffusionCoefficient(phaseIdx, compIdx)); };
    }

    //! return diffusion tensor/scalar lambda here (multi-phasic case)
    template<class T = TypeTag>
    static auto getDiffusionLambda(unsigned int phaseIdx,
                                   unsigned int compIdx,
                                   typename std::enable_if<(GET_PROP_VALUE(T, NumPhases) > 1), void>::type* dummy = nullptr)
    {
        return [phaseIdx, compIdx] (const auto& problem,
                                    const auto& element,
                                    const auto& volVars,
                                    const auto& fvGeometry,
                                    const auto& scv)
               { return volVars.diffusionCoefficient(phaseIdx, compIdx); };
    }

    //! return fourier coefficient/tensor lambda here
    static auto getHeatConductionLambda()
    {
        using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

        return [] (const auto& problem,
                   const auto& element,
                   const auto& volVars,
                   const auto& fvGeometry,
                   const auto& scv)
               { return ThermalConductivityModel::effectiveThermalConductivity(volVars,
                                                                               problem.spatialParams(),
                                                                               element,
                                                                               fvGeometry,
                                                                               scv); };
    }
};

} // end namespace

#endif
