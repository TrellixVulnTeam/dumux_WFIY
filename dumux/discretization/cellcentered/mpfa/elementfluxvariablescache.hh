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
 * \brief The local object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_ELEMENT_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>

#include "fluxvariablescachefiller.hh"

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the local flux variables cache.
 *        Prepares the cache on all the faces in the stencil.
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class CCMpfaElementFluxVariablesCache;

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    CCMpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf) {}

    // Specialization for the global caching being enabled - do nothing here
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*globalFluxVarsCachePtr_)[scvf.index()]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCMpfaElementFluxVariablesCache<TypeTag, false>
{
    // the flux variables cache filler needs to be friend to fill the matrices
    friend CCMpfaFluxVariablesCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using DataHandle = typename BoundaryInteractionVolume::Traits::DataHandle;

    static constexpr int dim = GridView::dimension;

public:
    CCMpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // This function has to be called prior to flux calculations on the element.
    // Prepares the transmissibilities of the scv faces in an element. The FvGeometry is assumed to be bound.
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        // TODO
        DUNE_THROW(Dune::NotImplemented, "Local element binding of the flux variables cache in mpfa schemes");
    }

    // This function is called by the CCLocalResidual before flux calculations during assembly.
    // Prepares the transmissibilities of the scv faces in the stencil. The FvGeometries are assumed to be bound.
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        // clear data
        clear_();

        // reference to the problem for convenience
        const auto& problem = globalFluxVarsCache().problem_();

        // the global index the cache will be bound to
        globalI_ = problem.elementMapper().index(element);

        // the assembly map of the given element
        const auto& assemblyMapI = problem.model().localJacobian().assemblyMap()[globalI_];

        // reserve memory for scvf index container
        unsigned int numNeighborScvfs = 0;
        for (auto&& dataJ : assemblyMapI)
            numNeighborScvfs += dataJ.scvfsJ.size();
        globalScvfIndices_.resize(fvGeometry.numScvf() + numNeighborScvfs);

        // set the scvf indices in scvf index container
        unsigned int i = 0;
        for (const auto& scvf : scvfs(fvGeometry))
            globalScvfIndices_[i++] = scvf.index();
        for (const auto& dataJ : assemblyMapI)
            for (auto scvfIdx : dataJ.scvfsJ)
                globalScvfIndices_[i++] = scvfIdx;

        // reserve memory estimate for interaction volumes and corresponding data
        const auto ivEstimate = getNoInteractionVolumesEstimate(element, assemblyMapI);
        const auto maxBoundaryIv = element.subEntities(dim);
        interactionVolumes_.reserve(ivEstimate);
        boundaryInteractionVolumes_.reserve(maxBoundaryIv);
        ivDataHandles_.reserve(ivEstimate);
        boundaryIvDataHandles_.reserve(maxBoundaryIv);

        // helper class to fill flux variables caches
        FluxVariablesCacheFiller filler(problem);

        // resize the cache container
        fluxVarsCache_.resize(globalScvfIndices_.size());

        // go through the caches and fill them
        i = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            auto& scvfCache = fluxVarsCache_[i++];
            if (!scvfCache.isUpdated())
                filler.fill(*this, scvfCache, element, fvGeometry, elemVolVars, scvf);
        }

        for (const auto& dataJ : assemblyMapI)
        {
            const auto elementJ = problem.model().globalFvGeometry().element(dataJ.globalJ);
            for (const auto scvfIdx : dataJ.scvfsJ)
            {
                auto& scvfCache = fluxVarsCache_[i++];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx));
            }
        }
    }

    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        // TODO
        DUNE_THROW(Dune::NotImplemented, "Local element binding of the flux variables cache in mpfa schemes");
    }

    // This function updates the transmissibilities (e.g. after the solution has been deflected)
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        // make sure this is not called for a different element than a previous bind()
        assert(globalI_ == globalFluxVarsCache().problem_().elementMapper().index(element)
               && "The element flux variables cache can only be updated for the element it was bound to");

        static const bool isSolDependent = FluxVariablesCacheFiller::isSolutionDependent();
        if (isSolDependent)
        {
            const auto& problem = globalFluxVarsCache().problem_();
            const auto& assemblyMapI = problem.model().localJacobian().assemblyMap()[globalI_];

            // helper class to fill flux variables caches
            FluxVariablesCacheFiller filler(problem);

            // set all the caches to "outdated"
            for (auto& cache : fluxVarsCache_)
                cache.setUpdateStatus(false);

            // go through the caches maybe update them
            unsigned int i = 0;
            for (const auto& scvf : scvfs(fvGeometry))
            {
                auto& scvfCache = fluxVarsCache_[i++];
                if (!scvfCache.isUpdated())
                    filler.update(*this, scvfCache, element, fvGeometry, elemVolVars, scvf);
            }

            for (const auto& dataJ : assemblyMapI)
            {
                const auto elementJ = problem.model().globalFvGeometry().element(dataJ.globalJ);
                for (const auto scvfIdx : dataJ.scvfsJ)
                {
                    auto& scvfCache = fluxVarsCache_[i++];
                    if (!scvfCache.isUpdated())
                        filler.update(*this, scvfCache, elementJ, fvGeometry, elemVolVars, fvGeometry.scvf(scvfIdx));
                }
            }
        }
    }

    // access operators in the case of no caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    const FluxVariablesCache& operator [](const IndexType scvfIdx) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    FluxVariablesCache& operator [](const IndexType scvfIdx)
    { return fluxVarsCache_[getLocalScvfIdx_(scvfIdx)]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;

    void clear_()
    {
        fluxVarsCache_.clear();
        globalScvfIndices_.clear();
        interactionVolumes_.clear();
        boundaryInteractionVolumes_.clear();
        ivDataHandles_.clear();
        boundaryIvDataHandles_.clear();
    }

    // get number of interaction volumes that are going to be required
    template<class AssemblyMap>
    std::size_t getNoInteractionVolumesEstimate(const Element& element, const AssemblyMap& assemblyMap)
    {
        if (GET_PROP_VALUE(TypeTag, MpfaMethod) == MpfaMethods::oMethod || GET_PROP_VALUE(TypeTag, MpfaMethod) == MpfaMethods::oMethodFps)
            return element.subEntities(dim);
        else if (GET_PROP_VALUE(TypeTag, MpfaMethod) == MpfaMethods::lMethod)
        {
            std::size_t numInsideScvfs = MpfaHelper::getNumLocalScvfs(element.geometry().type());
            std::size_t numOutsideScvf = 0;
            for (const auto& dataJ : assemblyMap) numOutsideScvf += dataJ.scvfsJ.size();
            return numOutsideScvf - numInsideScvfs;
        }
    }

    // get index of an scvf in the local container
    unsigned int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    // the global index of the actual bound element
    IndexType globalI_;

    // the local flux vars caches and the index set
    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;

    // store the interaction volumes and handles
    std::vector<InteractionVolume> interactionVolumes_;
    std::vector<BoundaryInteractionVolume> boundaryInteractionVolumes_;
    std::vector<DataHandle> ivDataHandles_;
    std::vector<DataHandle> boundaryIvDataHandles_;
};

} // end namespace

#endif
