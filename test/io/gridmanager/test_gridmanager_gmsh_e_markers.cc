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
 *
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include "gridmanagertests.hh"

int main(int argc, char** argv)
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    Dumux::Parameters::init(argc, argv, "test_gridmanager_gmsh_e_markers.input");

    auto name = Dumux::getParam<std::string>("Problem.Name");
    Dumux::GridManagerTests<GRIDTYPE>::testElementMarkers("gmsh", name);

    return 0;
}
