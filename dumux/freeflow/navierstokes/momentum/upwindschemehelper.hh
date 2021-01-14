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
 * \copydoc Dumux::MomentumUpwindSchemeHelper
 */
#ifndef DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH
#define DUMUX_MOMENTUM_UPWIND_SCHEME_HELPER_HH

#include <array>
#include <optional>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/method.hh>
#include <dumux/freeflow/staggeredupwindmethods.hh>
#include "velocitygradients.hh"

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The upwinding variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, int upwindSchemeOrder>
class FaceCenteredStaggeredUpwindHelper
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;

    using UpwindScheme = StaggeredUpwindMethods<Scalar, upwindSchemeOrder>;
    static constexpr bool useHigherOrder = upwindSchemeOrder > 1;
    static_assert(upwindSchemeOrder <= 2, "Not implemented: Order higher than 2!");

public:
    FaceCenteredStaggeredUpwindHelper(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const Problem& problem,
                                      const SubControlVolumeFace& scvf,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementBoundaryTypes& elemBcTypes,
                                      const UpwindScheme& upwindScheme)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , problem_(problem)
    , scvf_(scvf)
    , elemVolVars_(elemVolVars)
    , elemBcTypes_(elemBcTypes)
    , upwindScheme_(upwindScheme)
    {}

    /*!
     * \brief Returns the momentum in the frontal directon.
     *
     *        Checks if the model has higher order methods enabled and if the scvf in
     *        question is far enough from the boundary such that higher order methods can be employed.
     *        Then the corresponding set of momenta are collected and the prescribed
     *        upwinding method is used to calculate the momentum.
     */
    Scalar computeUpwindFrontalMomentum(const bool selfIsUpstream) const
    {
        const auto density = problem_.density(element_, fvGeometry_, scvf_);

        // for higher order schemes do higher order upwind reconstruction
        if constexpr (useHigherOrder)
        {
            // only second order is implemented so far
            if (canDoFrontalSecondOrder_(selfIsUpstream))
            {
                const auto distances = getFrontalDistances_(selfIsUpstream);
                const auto upwindMomenta = getFrontalSecondOrderUpwindMomenta_(density, selfIsUpstream);
                return upwindScheme_.tvd(upwindMomenta, distances, selfIsUpstream, upwindScheme_.tvdApproach());
            }
        }

        // otherwise apply first order upwind scheme
        const auto upwindMomenta = getFrontalFirstOrderUpwindMomenta_(density, selfIsUpstream);
        return upwindScheme_.upwind(upwindMomenta[0], upwindMomenta[1]);
    }
private:
    /*!
     * \brief Returns whether or not the face in question is far enough from the wall to handle higher order methods.
     *
     *        Evaluates which face is upstream.
     *        If the face is upstream, and the scvf has a forward neighbor, higher order methods are possible.
     *        If the face is downstream, and the scvf has a backwards neighbor, higher order methods are possible.
     *        Otherwise, higher order methods are not possible.
     */
    bool canDoFrontalSecondOrder_(bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        // Depending on selfIsUpstream I have to check if I have a forward or a backward neighbor to retrieve
        return selfIsUpstream ? fvGeometry_.hasForwardNeighbor(scvf_) : fvGeometry_.hasBackwardNeighbor(scvf_);
    }

    /*!
     * \brief Returns an array of the three momenta needed for second order upwinding methods.
     */
    std::array<Scalar, 3> getFrontalSecondOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");
        const auto momentumSelf = elemVolVars_[scvf_.insideScvIdx()].velocity() * density;
        const auto momentumOpposite = elemVolVars_[scvf_.outsideScvIdx()].velocity() * density;
        if (selfIsUpstream)
        {
            const auto momentumForward =  elemVolVars_[fvGeometry_.forwardScvIdx(scvf_)].velocity() * density;
            return {momentumOpposite, momentumSelf, momentumForward};
        }
        else
        {
            const auto momentumBackward = elemVolVars_[fvGeometry_.backwardScvIdx(scvf_)].velocity() * density;
            return {momentumSelf, momentumOpposite, momentumBackward};
        }
    }

    /*!
     * \brief Returns an array of the two momenta needed for first order upwinding method
     */
    std::array<Scalar, 2> getFrontalFirstOrderUpwindMomenta_(const Scalar density, bool selfIsUpstream) const
    {
        const auto momentumSelf = elemVolVars_[scvf_.insideScvIdx()].velocity() * density;
        const auto momentumOpposite = elemVolVars_[scvf_.outsideScvIdx()].velocity() * density;
        if (selfIsUpstream)
            return {momentumOpposite, momentumSelf};
        else
            return {momentumSelf, momentumOpposite};
    }

    /*!
     * \brief Returns an array of distances needed for non-uniform higher order upwind schemes
     * Depending on selfIsUpstream the downstream and the (up)upstream distances are saved.
     * distances {upstream to downstream distance, up-upstream to upstream distance, downstream staggered cell size}
     */
    std::array<Scalar, 3> getFrontalDistances_(const bool selfIsUpstream) const
    {
        static_assert(useHigherOrder, "Should only be reached if higher order methods are enabled");

        if (selfIsUpstream)
        {
            std::array<Scalar, 3> distances;
            distances[0] = fvGeometry_.selfToOppositeDistance(scvf_);
            distances[1] = fvGeometry_.selfToForwardDistance(scvf_);
            distances[2] = 0.5 * (fvGeometry_.selfToOppositeDistance(scvf_) + fvGeometry_.oppositeToBackwardDistance(scvf_));
            return distances;
        }
        else
        {
            std::array<Scalar, 3> distances;
            distances[0] = fvGeometry_.selfToOppositeDistance(scvf_);
            distances[1] = fvGeometry_.oppositeToBackwardDistance(scvf_);
            distances[2] = 0.5 * (fvGeometry_.selfToOppositeDistance(scvf_) + fvGeometry_.selfToForwardDistance(scvf_));
            return distances;
        }
    }

    const Element& element_;
    const FVElementGeometry& fvGeometry_;
    const Problem& problem_;
    const SubControlVolumeFace& scvf_;
    const ElementVolumeVariables& elemVolVars_;
    const ElementBoundaryTypes& elemBcTypes_;
    const UpwindScheme& upwindScheme_;
};

} // end namespace Dumux

#endif
