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
 * \ingroup MPNCModel
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 */
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_HH

#include <dune/istl/bvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup MPNCModel
 * \brief MpNc specific details needed to approximately calculate the local
 *        defect in the fully implicit scheme.
 *
 * This class is used to fill the gaps in ImplicitLocalResidual for the
 * MpNc flow.
 */
template<class TypeTag>
class MPNCLocalResidual : public CompositionalLocalResidual<TypeTag>
{
    using ParentType = CompositionalLocalResidual<TypeTag>;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {phase0NcpIdx = Indices::phase0NcpIdx};

public:
    using ParentType::ParentType;

    using typename ParentType::ElementResidualVector;

    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache,
                                            const ElementBoundaryTypes &bcTypes) const
    {
        ElementResidualVector residual = ParentType::evalFluxAndSource(element, fvGeometry, elemVolVars, elemFluxVarsCache, bcTypes);

        for (auto&& scv : scvs(fvGeometry))
        {
            //here we need to set the constraints of the mpnc model into the residual
            const auto localScvIdx = scv.indexInElement();
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                residual[localScvIdx][phase0NcpIdx + phaseIdx] = elemVolVars[scv].phaseNcp(phaseIdx);
        }

        return residual;
    }
};

} // end namespace

#endif
