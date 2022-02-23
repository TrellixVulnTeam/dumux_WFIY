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
 * \ingroup Typetraits
 * \copydoc Dumux::ConsistentlyOrientedGrid
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONSISTENTLY_ORIENTED_GRID_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_CONSISTENTLY_ORIENTED_GRID_HH

#include <type_traits>

// forward declare
namespace Dune {
template<int dim, class Coordinates>
class YaspGrid;

template <int dim, class HostGrid, bool mapIndexStorage>
class SubGrid;

}

namespace Dumux {

/*!
 * \brief Helper type to determine whether a grid is guaranteed to be oriented consistently.
 *        This means that the intersection indices always correspond to the ones of a reference element
 *        or, in other words, the elements are never rotated.
 */
template<class T>
struct ConsistentlyOrientedGrid : public std::false_type {};

template<int dim, class Coords>
struct ConsistentlyOrientedGrid<Dune::YaspGrid<dim, Coords>> : public std::true_type {};

template<int dim, class Coords, bool mapIndexStorage>
struct ConsistentlyOrientedGrid<Dune::SubGrid<dim, Dune::YaspGrid<dim, Coords>, mapIndexStorage>> : public std::true_type {};


} // end namespace Dumux

#endif
