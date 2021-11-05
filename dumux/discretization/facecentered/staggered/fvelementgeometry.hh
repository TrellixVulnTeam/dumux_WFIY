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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <utility>

#include <dune/common/rangeutilities.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>
#include <dumux/discretization/facecentered/staggered/consistentlyorientedgrid.hh>

namespace Dumux {

template<class GG, bool cachingEnabled>
class FaceCenteredStaggeredFVElementGeometry;

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for face-centered staggered models
 *        Specialization for grid caching enabled
 */
template<class GG>
class FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/true>
{
    using ThisType = FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    static constexpr auto numScvsPerElement = GG::StaticInformation::numScvsPerElement;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = numScvsPerElement;

    FaceCenteredStaggeredFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometry_(&gridGeometry)
    {}

    //! Get a sub control volume  with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    { return gridGeometry().scv(scvIdx); }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return gridGeometry().scvf(scvfIdx); }

    //! Return a the lateral sub control volume face which is orthogonal to the given sub control volume face
    const SubControlVolumeFace& lateralOrthogonalScvf(const SubControlVolumeFace& scvf) const
    {
        assert(scvf.isLateral());
        const auto otherGlobalIdx = scvfIndices_()[GridGeometry::GeometryHelper::lateralOrthogonalScvfLocalIndex(scvf.localIndex())];
        return gridGeometry().scvf(otherGlobalIdx);

    }

    //! Return the frontal sub control volume face on a the boundary for a given sub control volume
    const SubControlVolumeFace& frontalScvfOnBoundary(const SubControlVolume& scv) const
    {
        assert(scv.boundary());

        // frontal boundary faces are always stored after the lateral faces
        auto scvfIter = scvfs(*this, scv).begin();
        const auto end = scvfs(*this, scv).end();
        while (!(scvfIter->isFrontal() && scvfIter->boundary()) && (scvfIter != end))
            ++scvfIter;

        assert(scvfIter->isFrontal());
        assert(scvfIter->boundary());
        return *scvfIter;
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto
    scvs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    { return fvGeometry.gridGeometry().scvs(fvGeometry); }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    {
        using IndexContainerType = std::decay_t<decltype(fvGeometry.scvfIndices_())>;
        using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, IndexContainerType, ThisType>;
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(fvGeometry.scvfIndices_().begin(), fvGeometry),
                                                 ScvfIterator(fvGeometry.scvfIndices_().end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element belonging to the given sub control volume.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry, scv))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry, const SubControlVolume& scv)
    {
        using IndexContainerType = std::decay_t<decltype(fvGeometry.scvfIndices_())>;
        using ScvfIterator = Dumux::SkippingScvfIterator<SubControlVolumeFace, IndexContainerType, ThisType>;
        auto begin = ScvfIterator::makeBegin(fvGeometry.scvfIndices_(), fvGeometry, scv.index());
        auto end = ScvfIterator::makeEnd(fvGeometry.scvfIndices_(), fvGeometry, scv.index());
        return Dune::IteratorRange<ScvfIterator>(begin, end);
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    { return numScvsPerElement; }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    { return scvfIndices_().size(); }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(eIdx_); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    FaceCenteredStaggeredFVElementGeometry bind(const Element& element) &&
    {
        this->bind(element);
        return std::move(*this);
    }

    //! Binding of an element, called by the local jacobian to prepare element assembly
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    FaceCenteredStaggeredFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Bind only element-local
    void bindElement(const Element& element) &
    {
        elementPtr_ = &element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    {
        assert(gridGeometry_);
        return *gridGeometry_;
    }

    std::size_t elementIndex() const
    { return eIdx_; }

    //! The bound element
    const Element& element() const
    { return *elementPtr_; }

private:

    const auto& scvfIndices_() const
    {
        return gridGeometry().scvfIndicesOfElement(eIdx_);
    }

    const Element* elementPtr_;
    GridIndexType eIdx_;
    const GridGeometry* gridGeometry_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for face-centered staggered models
 *        Specialization for grid caching disabled
 */
template<class GG>
class FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/false>
{
    using ThisType = FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/false>;

    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    //TODO include assert that checks for quad geometry

    using StaticInfo = typename GG::StaticInformation;
    static constexpr auto dim = StaticInfo::dim;
    static constexpr auto numScvsPerElement = StaticInfo::numScvsPerElement;
    static constexpr auto numLateralScvfsPerScv = StaticInfo::numLateralScvfsPerScv;
    static constexpr auto numLateralScvfsPerElement = StaticInfo::numLateralScvfsPerElement;
    static constexpr auto minNumScvfsPerElement = StaticInfo::minNumScvfsPerElement;
    static constexpr auto maxNumScvfsPerElement = StaticInfo::maxNumScvfsPerElement;

    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = numScvsPerElement;

    FaceCenteredStaggeredFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometry_(&gridGeometry)
    , geometryHelper_(gridGeometry.gridView())
    {}

    //! Get a sub control volume face with a local scv index
    const SubControlVolumeFace& scvf(const GridIndexType scvfIdx) const
    { return scvfs_[findLocalIndex_(scvfIdx, scvfIndices_())]; }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    {
        using Iter = typename decltype(fvGeometry.scvfs_)::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element belonging to the given sub control volume.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry, scv))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry, const SubControlVolume& scv)
    {
        using IndexContainerType = std::decay_t<decltype(fvGeometry.scvfIndices_())>;
        using ScvfIterator = Dumux::SkippingScvfIterator<SubControlVolumeFace, IndexContainerType, ThisType>;
        auto begin = ScvfIterator::makeBegin(fvGeometry.scvfIndices_(), fvGeometry, scv.index());
        auto end = ScvfIterator::makeEnd(fvGeometry.scvfIndices_(), fvGeometry, scv.index());
        return Dune::IteratorRange<ScvfIterator>(begin, end);
    }

    //! Get a half sub control volume with a global scv index
    const SubControlVolume& scv(const GridIndexType scvIdx) const
    {
        if (scvIdx >= scvIndicesOfElement_.front() && scvIdx <= scvIndicesOfElement_.back())
            return scvs_[findLocalIndex_(scvIdx, scvIndicesOfElement_)];
        else
            return neighborScvs_[findLocalIndex_(scvIdx, neighborScvIndices_)];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto
    scvs(const FaceCenteredStaggeredFVElementGeometry& g)
    {
        using IteratorType = typename std::array<SubControlVolume, 1>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvs_.begin(), g.scvs_.end());
    }

    //! Return a the lateral sub control volume face which is orthogonal to the given sub control volume face
    const SubControlVolumeFace& lateralOrthogonalScvf(const SubControlVolumeFace& scvf) const
    {
        assert(scvf.isLateral());
        const auto otherLocalIdx = GridGeometry::GeometryHelper::lateralOrthogonalScvfLocalIndex(scvf.localIndex());
        return scvfs_[otherLocalIdx];
    }

    //! Return the frontal sub control volume face on a the boundary for a given sub control volume
    const SubControlVolumeFace& frontalScvfOnBoundary(const SubControlVolume& scv) const
    {
        assert(scv.boundary());

        // frontal boundary faces are always stored after the lateral faces
        auto scvfIter = scvfs(*this, scv).begin();
        const auto end = scvfs(*this, scv).end();
        while (!(scvfIter->isFrontal() && scvfIter->boundary()) && (scvfIter != end))
            ++scvfIter;

        assert(scvfIter->isFrontal());
        assert(scvfIter->boundary());
        return *scvfIter;
    }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(eIdx_); }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    { return numScvsPerElement; }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    { return scvfIndices_().size(); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    FaceCenteredStaggeredFVElementGeometry bind(const Element& element) &&
    {
        this->bind_(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    { this->bind_(element); }

    /*!
    * \brief bind the local view (r-value overload)
    * This overload is called when an instance of this class is a temporary in the usage context
    * This allows a usage like this: `const auto view = localView(...).bind(element);`
    */
    FaceCenteredStaggeredFVElementGeometry bindElement(const Element& element) &&
    {
        typename GG::LocalIntersectionMapper localIsMapper;
        localIsMapper.update(gridGeometry().gridView(), element);
        this->bindElement_(element, localIsMapper);
        return std::move(*this);
    }

    void bindElement(const Element& element) &
    {
        typename GG::LocalIntersectionMapper localIsMapper;
        localIsMapper.update(gridGeometry().gridView(), element);
        this->bindElement_(element, localIsMapper);
    }

    GridIndexType elementIndex() const
    { return eIdx_; }

    //! The bound element
    const Element& element() const
    { return *elementPtr_; }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    {
        assert(gridGeometry_);
        return *gridGeometry_;
    }

private:
    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind_(const Element& element)
    {
        typename GG::LocalIntersectionMapper localIsMapper;
        localIsMapper.update(gridGeometry().gridView(), element);

        bindElement_(element, localIsMapper);
        neighborScvIndices_.clear();
        neighborScvs_.clear();

        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (intersection.neighbor())
            {
                const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
                const auto localOppositeScvIdx = geometryHelper_.localOppositeIdx(localScvIdx);
                const auto& neighborElement = intersection.outside();
                const auto neighborElementIdx = gridGeometry().elementMapper().index(neighborElement);

                typename GG::LocalIntersectionMapper localNeighborIsMapper;
                localNeighborIsMapper.update(gridGeometry().gridView(), neighborElement);

                for (const auto& neighborIntersection : intersections(gridGeometry().gridView(), neighborElement))
                {
                    const auto localNeighborScvIdx = localNeighborIsMapper.realToRefIdx(neighborIntersection.indexInInside());
                    const bool addScvForFirstOrder = localNeighborScvIdx != localScvIdx && localNeighborScvIdx != localOppositeScvIdx;
                    const bool addScvForSecondOrder = GG::upwindSchemeOrder > 1 && localNeighborScvIdx == localScvIdx;
                    if (addScvForFirstOrder || addScvForSecondOrder)
                    {
                        const auto dofIndex = gridGeometry().intersectionMapper().globalIntersectionIndex(neighborElement, neighborIntersection.indexInInside());
                        const auto globalScvIndex = gridGeometry().globalScvIndex(neighborElementIdx, localNeighborScvIdx);
                        neighborScvs_.push_back(SubControlVolume(
                            neighborElement.geometry(),
                            neighborIntersection.geometry(),
                            globalScvIndex,
                            localNeighborScvIdx,
                            dofIndex,
                            Dumux::normalAxis(neighborIntersection.centerUnitOuterNormal()),
                            neighborElementIdx,
                            onDomainBoundary_(neighborIntersection)
                        ));

                        neighborScvIndices_.push_back(globalScvIndex);
                    }

                    if constexpr (GG::upwindSchemeOrder > 1)
                    {
                        // treat element right next to neighbor element
                        if (localNeighborScvIdx == localScvIdx)
                        {
                            const auto& secondNeighborElement = neighborIntersection.outside();
                            typename GG::LocalIntersectionMapper localSecondNeighborIsMapper;
                            localSecondNeighborIsMapper.update(gridGeometry().gridView(), secondNeighborElement);
                            for (const auto& secondNeighborIntersection : intersections(gridGeometry().gridView(), secondNeighborElement))
                            {
                                const auto localSecondNeighborScvIdx = localSecondNeighborIsMapper.realToRefIdx(secondNeighborIntersection.indexInInside());
                                if (localSecondNeighborScvIdx != localOppositeScvIdx)
                                {
                                    const auto dofIndex = gridGeometry().intersectionMapper().globalIntersectionIndex(secondNeighborElement, secondNeighborIntersection.indexInInside());
                                    const auto globalScvIndex = gridGeometry().globalScvIndex(secondNeighborElement, localSecondNeighborScvIdx);
                                    neighborScvs_.push_back(SubControlVolume(
                                        secondNeighborElement.geometry(),
                                        neighborIntersection.geometry(),
                                        globalScvIndex,
                                        localSecondNeighborScvIdx,
                                        dofIndex,
                                        Dumux::normalAxis(secondNeighborIntersection.centerUnitOuterNormal()),
                                        neighborElementIdx,
                                        onDomainBoundary_(secondNeighborIntersection)
                                    ));

                                    neighborScvIndices_.push_back(globalScvIndex);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //! Bind only element-local
    void bindElement_(const Element& element, const typename GG::LocalIntersectionMapper& localIsMapper)
    {
        elementPtr_ = &element;
        eIdx_ = gridGeometry().elementMapper().index(element);
        std::iota(scvIndicesOfElement_.begin(), scvIndicesOfElement_.end(), eIdx_*numScvsPerElement);
        scvfs_.clear();
        scvfs_.resize(minNumScvfsPerElement);

        const auto& elementGeometry = element.geometry();

        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
            auto localScvfIdx = localScvIdx*(1 + numLateralScvfsPerScv);
            const auto& intersectionGeometry = intersection.geometry();
            const auto dofIndex = gridGeometry().intersectionMapper().globalIntersectionIndex(element, intersection.indexInInside());

            scvs_[localScvIdx] = SubControlVolume(
                elementGeometry,
                intersectionGeometry,
                scvIndicesOfElement_[localScvIdx],
                localScvIdx,
                dofIndex,
                Dumux::normalAxis(intersection.centerUnitOuterNormal()),
                eIdx_,
                onDomainBoundary_(intersection)
            );

            // the frontal sub control volume face at the element center
            const auto localOppositeScvIdx = geometryHelper_.localOppositeIdx(localScvIdx);
            scvfs_[localScvfIdx] = SubControlVolumeFace(
                elementGeometry,
                intersectionGeometry,
                std::array{scvIndicesOfElement_[localScvIdx], scvIndicesOfElement_[localOppositeScvIdx]},
                localScvfIdx,
                scvfIndices_()[localScvfIdx],
                intersection.centerUnitOuterNormal(),
                SubControlVolumeFace::FaceType::frontal,
                SubControlVolumeFace::BoundaryType::interior
            );
            ++localScvfIdx;

            // the lateral sub control volume faces
            for (const auto lateralFacetIndex : Dune::transformedRangeView(geometryHelper_.localLaterFaceIndices(localScvIdx),
                                                                           [&](auto&& idx) { return localIsMapper.refToRealIdx(idx) ;})
                )
            {
                const auto& lateralIntersection = geometryHelper_.intersection(lateralFacetIndex, element);
                const auto& lateralFacet =  geometryHelper_.facet(lateralFacetIndex, element);
                const auto& lateralFacetGeometry = lateralFacet.geometry();

                // helper lambda to get the lateral scvf's global inside and outside scv indices
                const auto globalScvIndicesForLateralFace = [&]
                {
                    const auto globalOutsideScvIdx = [&]
                    {
                        if (lateralIntersection.neighbor())
                        {
                            const auto parallelElemIdx = gridGeometry().elementMapper().index(lateralIntersection.outside());
                            // todo: could be done easier?
                            std::array<GridIndexType, numScvsPerElement> globalScvIndicesOfNeighborElement;
                            std::iota(globalScvIndicesOfNeighborElement.begin(), globalScvIndicesOfNeighborElement.end(), parallelElemIdx*numScvsPerElement);
                            return globalScvIndicesOfNeighborElement[localScvIdx];
                        }
                        else
                            return gridGeometry().outsideVolVarIndex(scvfIndices_()[localScvfIdx]);
                    }();

                    return std::array{scvIndicesOfElement_[localScvIdx], globalOutsideScvIdx};
                }();

                const auto boundaryType = [&]
                {
                    if (onProcessorBoundary_(lateralIntersection))
                        return SubControlVolumeFace::BoundaryType::processorBoundary;
                    else if (onDomainBoundary_(lateralIntersection))
                        return SubControlVolumeFace::BoundaryType::physicalBoundary;
                    else
                        return SubControlVolumeFace::BoundaryType::interior;
                }();

                scvfs_[localScvfIdx] = SubControlVolumeFace(
                    elementGeometry,
                    intersectionGeometry,
                    lateralFacetGeometry,
                    globalScvIndicesForLateralFace, // TODO higher order
                    localScvfIdx,
                    scvfIndices_()[localScvfIdx],
                    lateralIntersection.centerUnitOuterNormal(),
                    SubControlVolumeFace::FaceType::lateral,
                    boundaryType
                );
                ++localScvfIdx;
            }
        }

        // do a second loop over all intersections to add frontal boundary faces
        auto localScvfIdx = minNumScvfsPerElement;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            // the frontal sub control volume face at a domain boundary (coincides with element face)
            if (onDomainBoundary_(intersection))
            {
                const auto localScvIdx = localIsMapper.realToRefIdx(intersection.indexInInside());
                // the frontal sub control volume face at the boundary
                scvfs_.push_back(SubControlVolumeFace(
                    element.geometry(),
                    intersection.geometry(),
                    std::array{scvIndicesOfElement_[localScvIdx], scvIndicesOfElement_[localScvIdx]},
                    localScvfIdx,
                    scvfIndices_()[localScvfIdx],
                    intersection.centerUnitOuterNormal(),
                    SubControlVolumeFace::FaceType::frontal,
                    SubControlVolumeFace::BoundaryType::physicalBoundary
                ));
                 ++localScvfIdx;
            }
        }

        if constexpr (!ConsistentlyOrientedGrid<typename GridView::Grid>{})
        {
            static const bool makeConsistentlyOriented = getParam<bool>("Grid.MakeConsistentlyOriented", true);
            if (makeConsistentlyOriented && scvfs_.size() > minNumScvfsPerElement)
            {
                // make sure frontal boundary scvfs are sorted correctly
                std::sort(scvfs_.begin() + minNumScvfsPerElement, scvfs_.end(),
                    [](const auto& scvfLeft, const auto& scvfRight)
                    { return scvfLeft.insideScvIdx() < scvfRight.insideScvIdx(); }
                );
            }
        }
    }

    const auto& scvfIndices_() const
    { return gridGeometry().scvfIndicesOfElement(eIdx_); }

    template<class Entry, class Container>
    const LocalIndexType findLocalIndex_(const Entry& entry,
                                         const Container& container) const
    {
        auto it = std::find(container.begin(), container.end(), entry);
        assert(it != container.end() && "Could not find the entry! Make sure to properly bind this class!");
        return std::distance(container.begin(), it);
    }

    bool onDomainBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && intersection.boundary();
    }

    bool onProcessorBoundary_(const typename GridView::Intersection& intersection) const
    {
        return !intersection.neighbor() && !intersection.boundary();
    }

    Dune::ReservedVector<SubControlVolumeFace, maxNumScvfsPerElement> scvfs_;

    Dune::ReservedVector<SubControlVolume, numLateralScvfsPerElement> neighborScvs_;
    Dune::ReservedVector<GridIndexType, numLateralScvfsPerElement> neighborScvIndices_;

    std::array<SubControlVolume, numScvsPerElement> scvs_;
    std::array<GridIndexType, numScvsPerElement> scvIndicesOfElement_;

    const GridGeometry* gridGeometry_;
    const Element* elementPtr_;
    GridIndexType eIdx_;
    typename GridGeometry::GeometryHelper geometryHelper_;
};

} // end namespace Dumux

#endif
