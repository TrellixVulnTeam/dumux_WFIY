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
  * \ingroup OneEqModel
  * \copydoc Dumux::OneEqResidualImpl
  */
#ifndef DUMUX_STAGGERED_ONEEQ_LOCAL_RESIDUAL_HH
#define DUMUX_STAGGERED_ONEEQ_LOCAL_RESIDUAL_HH

#include <dune/common/hybridutilities.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/localresidual.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class BaseLocalResidual, DiscretizationMethod discMethod>
class OneEqResidualImpl;

/*!
 * \ingroup OneEqModel
 * \brief Element-wise calculation of the residual for one-equation turbulence models
 *        using the staggered discretization
 */
template<class TypeTag, class BaseLocalResidual>
class OneEqResidualImpl<TypeTag, BaseLocalResidual, DiscretizationMethod::staggered>
: public BaseLocalResidual
{
    using ParentType = BaseLocalResidual;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static constexpr int viscosityTildeEqIdx = Indices::viscosityTildeEqIdx - ModelTraits::dim();

public:
    using ParentType::ParentType;

    //! Evaluate fluxes entering or leaving the cell center control volume.
    CellCenterPrimaryVariables computeStorageForCellCenter(const Problem& problem,
                                                           const SubControlVolume& scv,
                                                           const VolumeVariables& volVars) const
    {
        CellCenterPrimaryVariables storage = ParentType::computeStorageForCellCenter(problem, scv, volVars);
        storage[viscosityTildeEqIdx] = volVars.viscosityTilde();
        return storage;
    }

    CellCenterPrimaryVariables computeSourceForCellCenter(const Problem& problem,
                                                          const Element &element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const ElementVolumeVariables& elemVolVars,
                                                          const ElementFaceVariables& elemFaceVars,
                                                          const SubControlVolume &scv) const
    {
        CellCenterPrimaryVariables source = ParentType::computeSourceForCellCenter(problem, element, fvGeometry,
                                                                                   elemVolVars, elemFaceVars, scv);

        const auto& volVars = elemVolVars[scv];

        source[viscosityTildeEqIdx] += volVars.cb1() * (1.0 - volVars.ft2())
                                       * volVars.stressTensorScalarProductTilde()
                                       * volVars.viscosityTilde();

        source[viscosityTildeEqIdx] -= (volVars.cw1() * volVars.fW()
                                        - volVars.cb1() * volVars.ft2() / problem.karmanConstant() / problem.karmanConstant())
                                       * volVars.viscosityTilde() * volVars.viscosityTilde()
                                       / volVars.wallDistance() / volVars.wallDistance();

        for (unsigned int dimIdx = 0; dimIdx < ModelTraits::dim(); ++dimIdx)
        {
            source[viscosityTildeEqIdx] += volVars.cb2() / volVars.sigma()
                                           * volVars.storedViscosityTildeGradient()[dimIdx]
                                           * volVars.storedViscosityTildeGradient()[dimIdx];
        }

        return source;
    }
};
}

#endif
