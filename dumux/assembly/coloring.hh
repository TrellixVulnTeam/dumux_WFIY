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
 * \ingroup Assembly
 * \brief Coloring schemes for shared-memory-parallel assembly
 */
#ifndef DUMUX_ASSEMBLY_COLORING_HH
#define DUMUX_ASSEMBLY_COLORING_HH

#include <vector>
#include <deque>
#include <iostream>
#include <tuple>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
#include <dumux/discretization/method.hh>

#ifndef DOXYGEN // hide from doxygen
namespace Dumux::Detail {

//! Compute a map from dof indices to element indices (helper data for coloring algorithm)
template <class GridGeometry>
std::vector<std::vector<std::size_t>>
computeDofToElementMap(const GridGeometry& gg)
{
    std::vector<std::vector<std::size_t>> dofToElements;

    if constexpr (GridGeometry::discMethod == DiscretizationMethod::cctpfa)
    {
        dofToElements.resize(gg.gridView().size(0));
        const auto& eMapper = gg.elementMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = eMapper.index(element);
            for (const auto& intersection : intersections(gg.gridView(), element))
                if (intersection.neighbor())
                    dofToElements[eMapper.index(intersection.outside())].push_back(eIdx);
        }
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        static constexpr int dim = GridGeometry::GridView::dimension;
        dofToElements.resize(gg.gridView().size(dim));
        const auto& vMapper = gg.vertexMapper();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            for (int i = 0; i < element.subEntities(dim); i++)
                dofToElements[vMapper.subIndex(element, i, dim)].push_back(eIdx);
        }
    }

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method");

    return dofToElements;
}

/*!
 * \brief Compute the colors of neighboring nodes in the dependency graph
 *
 * Neighboring nodes are those elements that manipulate the
 * same data structures (e.g. system matrix, volvars, flux cache) in the same places
 *
 * \param gridGeometry the grid geometry
 * \param element the element we want to color
 * \param colors a vector of current colors for each element (not assigned: -1)
 * \param dofToElement a map from dof indices to element indices
 * \param neighborColors a vector to add the colors of neighbor nodes to
 */
template<class GridGeometry, class DofToElementMap>
void addNeighborColors(const GridGeometry& gg,
                       const typename GridGeometry::LocalView::Element& element,
                       const std::vector<int>& colors,
                       const DofToElementMap& dofToElement,
                       std::vector<int>& neighborColors)
{
    if constexpr (GridGeometry::discMethod == DiscretizationMethod::cctpfa)
    {
        // we modify neighbor elements during the assembly
        // check who else modifies these neighbor elements
        const auto& eMapper = gg.elementMapper();
        for (const auto& intersection : intersections(gg.gridView(), element))
        {
            if (intersection.neighbor())
            {
                // direct face neighbors
                const auto nIdx = eMapper.index(intersection.outside());
                neighborColors.push_back(colors[nIdx]);

                // neighbor-neighbors
                for (const auto nnIdx : dofToElement[eMapper.index(intersection.outside())])
                    neighborColors.push_back(colors[nnIdx]);
            }
        }
    }

    else if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        // we modify the vertex dofs of our element during the assembly
        // check who else modifies these vertex dofs
        const auto& vMapper = gg.vertexMapper();
        static constexpr int dim = GridGeometry::GridView::dimension;
        // direct vertex neighbors
        for (int i = 0; i < element.subEntities(dim); i++)
            for (auto eIdx : dofToElement[vMapper.subIndex(element, i, dim)])
                neighborColors.push_back(colors[eIdx]);
    }

    else
        DUNE_THROW(Dune::NotImplemented,
            "Missing coloring scheme implementation for this discretization method");
}

/*!
 * \brief Find the smallest color (integer >= 0) _not_ present in the given list of colors
 * \param colors list of colors which are already taken
 * \param notAssigned container to store which colors are not yet taken (is resized as required)
 */
int smallestAvailableColor(const std::vector<int>& colors,
                           std::vector<bool>& colorUsed)
{
    const int numColors = colors.size();
    colorUsed.assign(numColors, false);

    // The worst case for e.g. numColors=3 is colors={0, 1, 2}
    // in which case we return 3 as smallest available color
    // That means, we only track candidates in the (half-open) interval [0, numColors)
    // Mark candidate colors which are present in colors
    for (int i = 0; i < numColors; i++)
        if (colors[i] >= 0 && colors[i] < numColors)
            colorUsed[colors[i]] = true;

    // return smallest color not in colors
    for (int i = 0; i < numColors; i++)
        if (!colorUsed[i])
            return i;

    return numColors;
}

} // end namespace Dumux::Detail
#endif // DOXYGEN

namespace Dumux {

/*!
 * \brief Compute iterable lists of element seeds partitioned by color
 *
 * Splits up the elements of a grid view into partitions such that
 * all elements in one partition do not modify global data structures
 * at the same place during assembly. This is used to allow for
 * lock-free thread-parallel (shared memory) assembly routines.
 *
 * Implements a simply greedy graph coloring algorithm:
 * For each node (element), assign the smallest available color
 * not used by any of the neighboring nodes (element with conflicting memory access)
 * The greedy algorithm doesn't necessarily return the smallest
 * possible number of colors (that's a hard problem) but is fast
 *
 * \param gridGeometry the grid geometry
 * \param verbosity the verbosity level
 */
template<class GridGeometry>
auto coloredElementSets(const GridGeometry& gg, int verbosity = 1)
{
    Dune::Timer timer;

    using ElementSeed = typename GridGeometry::GridView::Grid::template Codim<0>::EntitySeed;;
    std::deque<std::vector<ElementSeed>> elementSets;
    std::vector<int> colors(gg.gridView().size(0), -1);

    // pre-reserve some memory for helper arrays to avoid reallocation
    std::vector<int> neighborColors; neighborColors.reserve(30);
    std::vector<bool> colorUsed; colorUsed.reserve(30);

    // dof to element map to speed up neighbor search
    const auto dofToElement = Detail::computeDofToElementMap(gg);

    for (const auto& element : elements(gg.gridView()))
    {
        // compute neighbor colors based on discretization-dependent stencil
        neighborColors.clear();
        Detail::addNeighborColors(gg, element, colors, dofToElement, neighborColors);

        // find smallest color (positive integer) not in neighborColors
        const auto color = Detail::smallestAvailableColor(neighborColors, colorUsed);

        // assign color to element
        colors[gg.elementMapper().index(element)] = color;

        // add element to the set of elements with the same color
        if (color < elementSets.size())
            elementSets[color].push_back(element.seed());
        else
            elementSets.push_back(std::vector<ElementSeed>{ element.seed() });
    }

    if (verbosity > 0)
        std::cout << Fmt::format("Colored {} elements with {} colors in {} seconds.\n",
                                 gg.gridView().size(0), elementSets.size(), timer.elapsed());

    return std::make_tuple(elementSets, colors);
}

//! Traits specifying if a given discretization tag supports coloring
template<DiscretizationMethod discMethod>
struct SupportsColoring : public std::false_type {};

template<> struct SupportsColoring<DiscretizationMethod::cctpfa> : public std::true_type {};
template<> struct SupportsColoring<DiscretizationMethod::box> : public std::true_type {};

} // end namespace Dumux

#endif
