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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index set within an interaction volume of the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH
#define DUMUX_DISCRETIZATION_MPFA_O_INTERACTIONVOLUME_INDEXSET_HH

#include <dune/common/reservedvector.hh>

#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux
{
/*!
 * \ingroup CCMpfaDiscretization
 * \brief The interaction volume index set class for the mpfa-o scheme.
 *
 * \tparam DualGridNodalIndexSet The type used for the nodal index set in the dual grid.
 */
template< class DualGridNodalIndexSet >
class CCMpfaOInteractionVolumeIndexSet
{
public:
    //! Export the type used for the nodal grid index sets
    using NodalIndexSet = DualGridNodalIndexSet;

    //! Export the types used for local/grid indices
    using LocalIndexType = typename DualGridNodalIndexSet::LocalIndexType;
    using GridIndexType = typename DualGridNodalIndexSet::GridIndexType;

    // Export the types used for local/grid stencils
    using LocalStencilType = typename DualGridNodalIndexSet::LocalStencilType;
    using GridStencilType = typename DualGridNodalIndexSet::GridStencilType;
    using GridScvfStencilType = typename DualGridNodalIndexSet::GridScvfStencilType;

    //! Export the type used for the neighbor scv index sets of the scvfs
    using ScvfNeighborLocalIndexSet = typename DualGridNodalIndexSet::ScvfNeighborLocalIndexSet;

    //! The constructor
    template< class FlipScvfIndexSet >
    CCMpfaOInteractionVolumeIndexSet(const NodalIndexSet& nodalIndexSet, const FlipScvfIndexSet& flipScvfIndexSet)
    : nodalIndexSet_(nodalIndexSet)
    {
        const auto numNodalScvfs = nodalIndexSet.numScvfs();

        // keeps track of which nodal scvfs have been handled already
        std::vector<bool> isHandled(numNodalScvfs, false);

        // go over faces in nodal index set, check if iv-local face has been
        // inserted already for this scvf and if not, insert index mapping
        numFaces_ = 0;
        for (LocalIndexType i = 0; i < numNodalScvfs; ++i)
        {
            // check if the nodal scvf still has to be handled
            if (isHandled[i])
                continue;

            // for scvfs touching the boundary there are no "outside" scvfs
            if (nodalIndexSet.scvfIsOnBoundary(i))
            {
                scvfNeighborScvLocalIndices_.push_back({nodalIndexSet.insideScvIdxLocal(i)});
                nodeToIvScvf_[i] = ivToNodeScvf_.size();
                isHandled[i] = true;
                ivToNodeScvf_.push_back(i);
                numFaces_++;
                continue;
            }

            // insert a new iv-local face
            const auto curIvLocalScvfIdx = ivToNodeScvf_.size();
            nodeToIvScvf_[i] = curIvLocalScvfIdx;
            isHandled[i] = true;

            // construct local index sets
            const auto& flipScvfIndices = flipScvfIndexSet[nodalIndexSet.scvfIdxGlobal(i)];
            const auto numFlipIndices = flipScvfIndices.size();

            ScvfNeighborLocalIndexSet neighborsLocal;
            neighborsLocal.resize(numFlipIndices + 1);
            neighborsLocal[0] = nodalIndexSet.insideScvIdxLocal(i);

            // mappings for all flip scvf
            for (unsigned int j = 0; j < numFlipIndices; ++j)
            {
                const auto outsideScvfIdx = flipScvfIndices[j];
                for (unsigned int nodeLocalIdx = 0; nodeLocalIdx < nodalIndexSet.numScvfs(); ++nodeLocalIdx)
                {
                    if (nodalIndexSet.scvfIdxGlobal(nodeLocalIdx) == outsideScvfIdx)
                    {
                        neighborsLocal[j+1] = nodalIndexSet.insideScvIdxLocal(nodeLocalIdx);
                        nodeToIvScvf_[nodeLocalIdx] = curIvLocalScvfIdx;
                        isHandled[nodeLocalIdx] = true;
                        break; // go to next outside scvf
                    }
                }
            }

            scvfNeighborScvLocalIndices_.push_back(neighborsLocal);
            ivToNodeScvf_.push_back(i);
            numFaces_++;
        }
    }

    //! returns the corresponding nodal index set
    const DualGridNodalIndexSet& nodalIndexSet() const { return nodalIndexSet_; }

    //! returns the global scv indices connected to this dual grid node
    const GridStencilType& globalScvIndices() const { return nodalIndexSet_.globalScvIndices(); }

    //! returns the global scvf indices connected to this dual grid node
    const GridScvfStencilType& globalScvfIndices() const { return nodalIndexSet_.globalScvfIndices(); }

    //! returns the number of faces in the interaction volume
    std::size_t numFaces() const { return numFaces_; }

    //! returns the number of scvs in the interaction volume
    std::size_t numScvs() const { return nodalIndexSet_.numScvs(); }

    //! returns a global scvf idx for a given iv-local scvf index
    GridIndexType scvfIdxGlobal(LocalIndexType ivLocalScvfIdx) const
    {
        assert(ivLocalScvfIdx < numFaces());
        return nodalIndexSet_.scvfIdxGlobal( ivToNodeScvf_[ivLocalScvfIdx] );
    }

    //! returns the iv-local scvf idx of the i-th scvf embedded in a local scv
    LocalIndexType scvfIdxLocal(LocalIndexType scvIdxLocal, unsigned int i) const
    {
        assert(nodalIndexSet_.scvfIdxLocal(scvIdxLocal, i) < nodeToIvScvf_.size());
        return nodeToIvScvf_[ nodalIndexSet_.scvfIdxLocal(scvIdxLocal, i) ];
    }

    //! returns the local indices of the neighboring scvs of an scvf
    const ScvfNeighborLocalIndexSet& neighboringLocalScvIndices(LocalIndexType ivLocalScvfIdx) const
    {
        assert(ivLocalScvfIdx < numFaces());
        return scvfNeighborScvLocalIndices_[ivLocalScvfIdx];
    }

private:
    //! returns the local scv index to a given global scv index
    unsigned int findLocalScvIdx_(GridIndexType globalScvIdx) const
    {
        auto it = std::find( nodalIndexSet_.globalScvIndices().begin(), nodalIndexSet_.globalScvIndices().end(), globalScvIdx );
        assert(it != nodalIndexSet_.globalScvIndices().end() && "Global scv index not found in local container!");
        return std::distance(nodalIndexSet_.globalScvIndices().begin(), it);
    }

    const DualGridNodalIndexSet& nodalIndexSet_;

    std::size_t numFaces_;
    Dune::ReservedVector< LocalIndexType, NodalIndexSet::maxNumScvfsAtNode > ivToNodeScvf_;
    Dune::ReservedVector< LocalIndexType, NodalIndexSet::maxNumScvfsAtNode > nodeToIvScvf_;

    // maps to each scvf a list of neighbouring scv indices
    // ordering: 0 - inside scv idx; 1..n - outside scv indices
    Dune::ReservedVector< ScvfNeighborLocalIndexSet, NodalIndexSet::maxNumScvfsAtNode > scvfNeighborScvLocalIndices_;
};

} // end namespace Dumux

#endif
