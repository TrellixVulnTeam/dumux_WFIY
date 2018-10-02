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
 * \ingroup ThreePModel
 * \brief Adds vtk output fields specific to the three-phase model
 */
#ifndef DUMUX_THREEP_VTK_OUTPUT_FIELDS_HH
#define DUMUX_THREEP_VTK_OUTPUT_FIELDS_HH

namespace Dumux {

/*!
 * \ingroup ThreePModel
 * \brief Adds vtk output fields specific to the three-phase model
 */
class ThreePVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized vtk output fields
        for (int i = 0; i < VolumeVariables::numPhases(); ++i)
        {
            vtk.addVolumeVariable([i](const auto& v){ return v.saturation(i); }, "S_"+ FluidSystem::phaseName(i));
            vtk.addVolumeVariable([i](const auto& v){ return v.pressure(i); }, "p_"+ FluidSystem::phaseName(i));
            vtk.addVolumeVariable([i](const auto& v){ return v.density(i); }, "rho_"+ FluidSystem::phaseName(i));
        }

        vtk.addVolumeVariable( [](const auto& v){ return v.porosity(); },"porosity");
        vtk.addVolumeVariable( [](const auto& v){ return v.permeability(); },"permeability");
    }
};

} // end namespace Dumux

#endif
