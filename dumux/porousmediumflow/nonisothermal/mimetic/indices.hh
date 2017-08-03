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
 * \brief Defines the indices used by the non-isothermal two-phase two-component model
 */
#ifndef DUMUX_MIMETIC_ENERGY_INDICES_HH
#define DUMUX_MIMETIC_ENERGY_INDICES_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup NIModel
 * \ingroup ImplicitIndices
 * \brief Indices for the non-isothermal two-phase two-component model
 *
 * \tparam formulation The formulation, either pwsn or pnsw.
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
class MimeticEnergyIndices : public GET_PROP_TYPE(TypeTag, IsothermalIndices)
{
public:
    static const int numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
    static const int temperatureIdx = PVOffset + numEqCellCenter - 1; //! The index for temperature in primary variable vectors.
    static const int energyEqIdx = PVOffset + numEqCellCenter - 1; //! The index for energy in equation vectors.

    static const int numEqFace = GET_PROP_VALUE(TypeTag, NumEqFace);
    static constexpr int faceFluxBalanceEnergyIdx = energyEqIdx + numEqFace;
    static constexpr int faceTemperatureIdx = faceFluxBalanceEnergyIdx;

};

} // end namespace Dumux

#endif
