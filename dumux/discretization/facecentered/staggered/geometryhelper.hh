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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredFVGridGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_GEOMETRY_HELPER_HH

#include <dune/common/float_cmp.hh>
#include <dumux/common/indextraits.hh>
#include <dune/common/reservedvector.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/method.hh>

namespace Dune {

template<int dim, class Coordinates>
class YaspGrid;

}

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Helper used for navigating through staggered elements and sub Entities
 *
 * \image html localScvfToElementMappingStaggered.svg
 */
template<class GridView, class Implementation>
class FaceCenteredStaggeredGeometryHelperBase
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    using Facet = typename GridView::template Codim<1>::Entity;

public:
    static constexpr auto dim = GridView::Grid::dimension;
    static constexpr auto numElementFaces = dim * 2;
    static constexpr auto numLateralFacesPerScv = 2 * (dim - 1);

    FaceCenteredStaggeredGeometryHelperBase(const GridView& gridView) : gridView_(gridView) {}

    //! Returns the local index of the opposing face.
    static constexpr SmallLocalIndexType localOppositeIdx(const SmallLocalIndexType ownLocalFaceIndex)
    {
        return (ownLocalFaceIndex % 2) ? (ownLocalFaceIndex - 1) : (ownLocalFaceIndex + 1);
    }

     //! Returns the local indices of the faces lateral to the own one.
    static constexpr auto localLaterFaceIndices(const SmallLocalIndexType ownLocalFaceIndex)
    {
        constexpr auto table = []
        {
            using Table = std::array<std::array<SmallLocalIndexType, numLateralFacesPerScv>, numElementFaces>;
            if constexpr (dim == 1)
                return Table{};
            else if constexpr (dim == 2)
                return Table {{ {2,3}, {2,3}, {0,1}, {0,1} }};
            else
                return Table {{ {2,3,4,5}, {2,3,4,5}, {0,1,4,5}, {0,1,4,5}, {0,1,2,3}, {0,1,2,3} }};
        }();

        return table[ownLocalFaceIndex];
    }

    //! Returns an element's facet based on the local facet index.
    Facet getFacet(const SmallLocalIndexType localFacetIdx, const Element& element) const
    {
        return element.template subEntity <1> (localFacetIdx);
    }

    //! Returns an element's intersection based on the local facet index.
    auto getIntersection(const SmallLocalIndexType localFacetIdx, const Element& element) const
    {
        SmallLocalIndexType counter = 0;
        for (const auto& intersection : intersections(gridView(), element))
        {
            if (counter == localFacetIdx)
                return intersection;
            ++counter;
        }
        DUNE_THROW(Dune::InvalidStateException, "localFacetIdx " << localFacetIdx << " out of range");
    }

    //! Map two local facet indices to a local common entity (vertex in 2D or edge in 3D) index
    SmallLocalIndexType getLocalCommonEntityIndex(SmallLocalIndexType localFacetIdx0, SmallLocalIndexType localFacetIdx1) const
    {
        constexpr auto table = []
        {
            using Table = std::array<std::array<SmallLocalIndexType, 2>, dim == 2 ? 4 : 12>;
            if constexpr (dim == 2)
                return Table {{ {0,2}, {1,2}, {0,3}, {1,3} }};
            else
                return Table {{ {0,2}, {1,2}, {0,3}, {1,3}, {0,4}, {1,4}, {2,4}, {3,4}, {0,5}, {1,5}, {2,5}, {3,5} }};
        }();

        using std::swap;
        if (localFacetIdx0 > localFacetIdx1)
            swap(localFacetIdx0, localFacetIdx1);

        const auto idx = std::find(table.begin(), table.end(), std::array<SmallLocalIndexType, 2>{localFacetIdx0, localFacetIdx1});

        if (idx != std::end(table))
            return std::distance(table.begin(), idx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Faces with local indices " << localFacetIdx0 <<  " and "  << localFacetIdx1 << " do not intersect");
    }

    const GridView& gridView() const
    { return gridView_; }

protected:

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:

    const GridView gridView_;
};


template<class GridView, class Grid>
class FaceCenteredStaggeredGeometryHelper
    : public FaceCenteredStaggeredGeometryHelperBase<GridView, FaceCenteredStaggeredGeometryHelper<GridView, Grid>>
{
    using ParentType = FaceCenteredStaggeredGeometryHelperBase<GridView, FaceCenteredStaggeredGeometryHelper<GridView, Grid>>;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
public:

    using ParentType::ParentType;

    //! Update the map for getting the corresponding local face indices in another element.
    void update(const Element& ownElement, const Element& otherElement)
    {
        SmallLocalIndexType ownIdx = 0;
        SmallLocalIndexType otherIdx = 0;
        for (const auto ownIs : intersections(this->gridView(), ownElement))
        {
            // helper lambda to make sure the inner loop stops once the other index is found
            auto getOtherIdx = [&] (const auto& ownOuterNormal)
            {
                for (const auto& otherIs : intersections(this->gridView(), otherElement))
                {
                    const auto& otherOuterNormal = otherIs.centerUnitOuterNormal();
                    if (Dune::FloatCmp::eq<typename GridView::ctype>(ownOuterNormal*otherOuterNormal, 1.0))
                    {
                        map_[ownIdx++] = otherIdx++;
                        return;
                    }
                }
            };

            getOtherIdx(ownIs.centerUnitOuterNormal());
        }
        if(ownIdx != map_.size())
            DUNE_THROW(Dune::InvalidStateException, "Index could not be mapped");
    }

    //! Return the local index of the other element's facet with the same position as the own element's facet.
    SmallLocalIndexType localFaceIndexInOtherElement(const SmallLocalIndexType localFaceIndexInOwnElement) const
    {
        return map_[localFaceIndexInOwnElement];
    }

private:
    std::array<SmallLocalIndexType, ParentType::numElementFaces> map_;
};

template<class GridView>
class FaceCenteredStaggeredGeometryHelper<GridView, Dune::YaspGrid<GridView::Grid::dimension, typename GridView::Grid::ctype>>
    : public FaceCenteredStaggeredGeometryHelperBase<GridView, FaceCenteredStaggeredGeometryHelper<GridView, typename GridView::Grid>>
{
    using ParentType = FaceCenteredStaggeredGeometryHelperBase<GridView, FaceCenteredStaggeredGeometryHelper<GridView, typename GridView::Grid>>;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
public:

    using ParentType::ParentType;

    //! Update the map for getting the corresponding local face indices in another element.
    //! Nothing needs to be done here.
    void update(const Element& ownElement, const Element& otherElement)
    {}

    //! Return the local index of the other element's facet with the same position as the own element's facet.
    //! For Yasp grids, this is just the same index.
    SmallLocalIndexType localFaceIndexInOtherElement(const SmallLocalIndexType localFaceIndexInOwnElement) const
    {
        return localFaceIndexInOwnElement;
    }
};

} // end namespace Dumux
#endif
