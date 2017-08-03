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
 * \brief This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_DARCYS_LAW_HH

#include <dumux/discretization/methods.hh>

namespace Dumux
{
// forward declaration
template <class TypeTag, DiscretizationMethods Method>
class DarcysLawImplementation
{};

/*!
 * \ingroup DarcysLaw
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face. Specializations are provided
 * for the different discretization methods. These specializations
 * are found in the headers included below.
 */
template <class TypeTag>
using DarcysLaw = DarcysLawImplementation<TypeTag, GET_PROP_VALUE(TypeTag, DiscretizationMethod)>;

} // end namespace

#include <dumux/discretization/box/darcyslaw.hh>
#include <dumux/discretization/cellcentered/tpfa/darcyslaw.hh>
#include <dumux/discretization/cellcentered/mpfa/darcyslaw.hh>
#include <dumux/discretization/staggered/mimetic/darcyslaw.hh>

#endif
