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
 * \ingroup TracerModel
 * \brief Adds vtk output fields specific to the tracer model
 */
#ifndef DUMUX_TRACER_VTK_OUTPUT_FIELDS_HH
#define DUMUX_TRACER_VTK_OUTPUT_FIELDS_HH

#include <string>

namespace Dumux {

/*!
 * \ingroup TracerModel
 * \brief Adds vtk output fields specific to the tracer model
 */
class TracerVtkOutputFields
{
public:
    template <class VtkOutputModule>
    static void init(VtkOutputModule& vtk)
    {
        using VolumeVariables = typename VtkOutputModule::VolumeVariables;
        using FluidSystem = typename VolumeVariables::FluidSystem;

        // register standardized vtk output fields
        for (int phaseIdx = 0; phaseIdx < VolumeVariables::numPhases(); ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < VolumeVariables::numComponents(); ++compIdx)
            {
                vtk.addVolumeVariable( [phaseIdx, compIdx](const auto& v){  return v.moleFraction(phaseIdx, compIdx); },
                        "x_" + std::string(FluidSystem::phaseName(phaseIdx) + "^" + FluidSystem::componentName(compIdx)));
                vtk.addVolumeVariable( [phaseIdx, compIdx](const auto& v){  return v.massFraction(phaseIdx, compIdx); },
                        "X_" + std::string(FluidSystem::phaseName(phaseIdx) + "^" + FluidSystem::componentName(compIdx)));
            }
            vtk.addVolumeVariable( [](const auto& v){ return v.density(); }, "rho");
        }
    }
};

} // end namespace Dumux

#endif
