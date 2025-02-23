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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredFreeFlowVelocityOutput
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_VELOCITYOUTPUT_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_VELOCITYOUTPUT_HH

#include <dumux/io/velocityoutput.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityreconstruction.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Velocity output for staggered free-flow models
 */
template<class GridVariables>
class NavierStokesVelocityOutput : public VelocityOutput<GridVariables>
{
    using ParentType = VelocityOutput<GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FieldType = typename ParentType::FieldType;

public:
    using VelocityVector = typename ParentType::VelocityVector;

    NavierStokesVelocityOutput(const std::string& paramGroup = "")
    {
        enableOutput_ = getParamFromGroup<bool>(paramGroup, "Vtk.AddVelocity", true);
    }

    //! Returns whether to enable the velocity output or not
    bool enableOutput() const override { return enableOutput_; }

    //! returns the phase name of a given phase index
    std::string phaseName(int phaseIdx) const override { return FluidSystem::phaseName(phaseIdx); }

    //! returns the number of phases
    int numFluidPhases() const override { return VolumeVariables::numFluidPhases(); }

    //! returns the field type
    FieldType fieldType() const override { return FieldType::element; }

    //! Calculate the velocities for the scvs in the element
    //! We assume the local containers to be bound to the complete stencil
    void calculateVelocity(VelocityVector& velocity,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVarsCache& elemFluxVarsCache,
                           int phaseIdx) const override
    {
        using CouplingManager = std::decay_t<decltype(elemVolVars.gridVolVars().problem().couplingManager())>;
        using MomGG = std::decay_t<decltype(std::declval<CouplingManager>().problem(CouplingManager::freeFlowMomentumIndex).gridGeometry())>;
        if constexpr (MomGG::discMethod == DiscretizationMethods::fcstaggered)
            calculateVelocityForStaggeredGrid_(velocity, element, fvGeometry, elemVolVars);
        else if constexpr (
            MomGG::discMethod == DiscretizationMethods::fcdiamond
            || MomGG::discMethod == DiscretizationMethods::pq1bubble
            || MomGG::discMethod == DiscretizationMethods::box
        )
            calculateVelocityForCVFESchemes_(velocity, element, fvGeometry, elemVolVars);
    }

private:
    void calculateVelocityForStaggeredGrid_(VelocityVector& velocity,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto getFaceVelocity = [&](const FVElementGeometry& fvG, const auto& scvf)
        {
            return elemVolVars.gridVolVars().problem().faceVelocity(element, fvGeometry, scvf);
        };

        velocity[eIdx] = StaggeredVelocityReconstruction::cellCenterVelocity(getFaceVelocity, fvGeometry);
    }

    void calculateVelocityForCVFESchemes_(VelocityVector& velocity,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        velocity[eIdx] = elemVolVars.gridVolVars().problem().elementVelocity(fvGeometry);
    }


    bool enableOutput_;
};

} // end namespace Dumux

#endif
