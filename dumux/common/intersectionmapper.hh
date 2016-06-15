// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
#ifndef DUMUX_INTERSECTIONITERATOR_HH
#define DUMUX_INTERSECTIONITERATOR_HH

#include <vector>
#include <unordered_map>

#include <dune/common/deprecated.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

/*!
 * \file
 * \brief defines an intersection mapper for mapping of global DOFs assigned
 *        to faces which also works for adaptive grids.
 */

namespace Dumux
{

/*!
 * \brief defines an intersection mapper for mapping of global DOFs assigned
 *        to faces which also works for adaptive grids.
 */
template<class GridView>
class IntersectionMapper
{
    typedef typename GridView::Grid Grid;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView, Dune::MCMGElementLayout> ElementMapper;

public:
    IntersectionMapper(const GridView& gridview)
    : gridView_(gridview), elementMapper_(gridView_), size_(gridView_.size(1)),
      intersectionMapGlobal_(gridView_.size(0)), intersectionMapLocal_(gridView_.size(0))
    {
//        const auto element = *gridView_.template begin<0>();
//
//        int fIdx = 0;
//        for (const auto& intersection : intersections(gridView_, element))
//        {
//            int idxInInside = intersection.indexInInside();
//
//            standardLocalIdxMap_[idxInInside] = fIdx;
//
//            fIdx++;
//        }
    }

    const ElementMapper& elementMapper() const
    {
        return elementMapper_;
    }

    DUNE_DEPRECATED_MSG("Map is deprecated in 2.9, use index instead.")
    int map(const Element& element) const
    {
        return elementMapper_.index(element);
    }

    /*!
     * \brief Map element to array index.
     *
     * \param element Reference to element which should be mapped.
     * \return An index in the range 0 ... Max number of elements in set - 1.
     */
    int index(const Element& element) const
    {
        return elementMapper_.index(element);
    }

    DUNE_DEPRECATED_MSG("Map is deprecated in 2.9, use subIndex instead.")
    int map(int elemIdx, int fIdx)
    {
        return intersectionMapGlobal_[elemIdx][fIdx];
    }

    /*!
     * \brief Map interface fIdx'th interface of element index to array index.
     *
     * \param elemIdx Index of element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(int elemIdx, int fIdx)
    {
        return intersectionMapGlobal_[elemIdx][fIdx];
    }

    DUNE_DEPRECATED_MSG("Map is deprecated in 2.9, use subIndex instead.")
    int map(int elemIdx, int fIdx) const
    {
        return (intersectionMapGlobal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }

    /*!
     * \brief Map interface fIdx'th interface of element index to array index.
     *
     * \param elemIdx Index of element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(int elemIdx, int fIdx) const
    {
        return (intersectionMapGlobal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }

    DUNE_DEPRECATED_MSG("Map is deprecated in 2.9, use subIndex instead.")
    int map(const Element& element, int fIdx)
    {
        return intersectionMapGlobal_[index(element)][fIdx];
    }

    /*!
     * \brief Map interface fIdx'th interface of element to array index.
     *
     * \param element Reference to element.
     * \param fIdx Index of interface to map.
     * \return An index in the range 0 ... Max number of intersections in set - 1.
     */
    int subIndex(const Element& element, int fIdx)
    {
        return intersectionMapGlobal_[index(element)][fIdx];
    }

    DUNE_DEPRECATED_MSG("Map is deprecated in 2.9, use subIndex instead.")
    int map(const Element& element, int fIdx) const
    {
        return intersectionMapGlobal_[index(element)].find(fIdx)->second;//use find() for const function!
    }

    int subIndex(const Element& element, int fIdx) const
    {
        return intersectionMapGlobal_[index(element)].find(fIdx)->second;//use find() for const function!
    }

    int maplocal(int elemIdx, int fIdx)
    {
        return intersectionMapLocal_[elemIdx][fIdx];
    }

    int maplocal(int elemIdx, int fIdx) const
    {
        return (intersectionMapLocal_[elemIdx].find(fIdx))->second;//use find() for const function!
    }


    int maplocal(const Element& element, int fIdx)
    {
        return intersectionMapLocal_[index(element)][fIdx];
    }

    int maplocal(const Element& element, int fIdx) const
    {
        return (intersectionMapLocal_[index(element)].find(fIdx))->second;//use find() for const function!
    }

    // return number intersections
    unsigned int size () const
    {
        return size_;
    }

    unsigned int size (int elemIdx) const
    {
        return intersectionMapLocal_[elemIdx].size();
    }

    unsigned int size (const Element& element) const
    {
        return intersectionMapLocal_[index(element)].size();
    }

    void update()
    {
        elementMapper_.update();

        intersectionMapGlobal_.clear();
        intersectionMapGlobal_.resize(elementMapper_.size());
        intersectionMapLocal_.clear();
        intersectionMapLocal_.resize(elementMapper_.size());

        for (const auto& element : elements(gridView_))
        {
            int eIdxGlobal = index(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : intersections(gridView_, element))
            {
                int indexInInside = intersection.indexInInside();
                intersectionMapLocal_[eIdxGlobal][fIdx] = indexInInside;

                fIdx++;
            }
        }

        int globalIntersectionIdx = 0;
        for (const auto& element : elements(gridView_))
        {
            int eIdxGlobal = index(element);

            int fIdx = 0;
            // run through all intersections with neighbors
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.neighbor())
                {
                    auto neighbor = intersection.outside();
                    int globalIdxNeighbor = index(neighbor);

                    if (element.level() > neighbor.level() || (element.level() == neighbor.level() && eIdxGlobal < globalIdxNeighbor))
                    {

                        int faceIdxNeighbor = 0;
//                        if (size(globalIdxNeighbor) == 2 * dim)
//                        {
//                            faceIdxNeighbor = standardLocalIdxMap_[intersection.indexInOutside()];
//                        }
//                        else
//                        {
                            for (const auto& intersectionNeighbor : intersections(gridView_, neighbor))
                            {
                                if (intersectionNeighbor.neighbor())
                                {
                                    if (intersectionNeighbor.outside() == element)
                                    {
                                        break;
                                    }
                                }
                                faceIdxNeighbor++;
                            }
//                        }

                        intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                        intersectionMapGlobal_[globalIdxNeighbor][faceIdxNeighbor] = globalIntersectionIdx;
                        globalIntersectionIdx ++;
                    }
                }
                else
                {
                    intersectionMapGlobal_[eIdxGlobal][fIdx] = globalIntersectionIdx;
                    globalIntersectionIdx ++;
                }
                fIdx++;
            }
        }
        size_ = globalIntersectionIdx;
    }

protected:
    const GridView gridView_;
    ElementMapper elementMapper_;
    unsigned int size_;
    std::vector<std::unordered_map<int, int> > intersectionMapGlobal_;
    std::vector<std::unordered_map<int, int> > intersectionMapLocal_;
    //std::unordered_map<int, int> standardLocalIdxMap_;
};

}

#endif
