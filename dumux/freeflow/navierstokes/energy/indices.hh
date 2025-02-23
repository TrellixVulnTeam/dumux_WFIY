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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::NavierStokesEnergyIndices
 */
#ifndef DUMUX_FREEFLOW_NAVIER_STOKES_ENERGY_INDICES_HH
#define DUMUX_FREEFLOW_NAVIER_STOKES_ENERGY_INDICES_HH

namespace Dumux {

/*!
 * \ingroup FreeflowNIModel
 * \brief Indices for the non-isothermal Navier-Stokes model.
 *
 * \tparam IsothermalIndices The isothermal indices class
 * \tparam numEq the number of equations of the non-isothermal model
 */
template <class IsothermalIndices, int numEq>
class NavierStokesEnergyIndices : public IsothermalIndices
{
public:
    static constexpr int energyEqIdx = numEq - 1;
    static constexpr int temperatureIdx = numEq - 1;
};

} // end namespace Dumux

#endif
