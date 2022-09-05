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
 * \ingroup Common
 * \brief A Python-like enumerate function
 */

#ifndef DUMUX_COMMON_ENUMERATE_HH
#define DUMUX_COMMON_ENUMERATE_HH

#include <tuple>
#include <ranges>

namespace Dumux {

/*!
 * \brief A Python-like enumerate function
 * \param inputRange Range to be enumerated
 * Usage example: for (const auto& [i, item] : enumerate(list))
 */
template<std::ranges::range Range>
constexpr auto enumerate(Range&& inputRange)
{
    if constexpr (std::is_reference_v<std::ranges::range_reference_t<Range>>)
    {
        using Ref = std::ranges::range_reference_t<Range>;
        return std::views::transform(inputRange, [i=0] (Ref r) mutable {
            return std::tie(i, r);
        });
    }
    else
    {
        using Value = std::ranges::range_value_t<Range>;
        static_assert(std::is_move_constructible_v<Value>);
        return std::views::transform(inputRange, [i=0] (Value v) mutable {
            return std::make_tuple(i, std::move(v));
        });
    }
}

} // end namespace Dumux

#endif
