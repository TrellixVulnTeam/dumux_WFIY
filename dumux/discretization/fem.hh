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
 * \ingroup Discretization
 * \brief Declares properties required for models using finite element schemes.
 */

#ifndef DUMUX_FEM_PROPERTIES_HH
#define DUMUX_FEM_PROPERTIES_HH

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/fem/fegridvariables.hh>
#include <dumux/discretization/fem/secondaryvariablesbase.hh>
#include <dumux/discretization/fem/elementboundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief Traits class for the base class of secondary variables.
 *
 * \tparam PV The type used for primary variables
 */
template<class PV>
struct SecondaryVariablesBaseTraits
{
    using PrimaryVariables = PV;
};

namespace Properties {

//! Type tag for finite-volume schemes.
// Create new type tags
namespace TTag {
//! \todo TODO Call FEModel to be more consistent with other type tags
struct FiniteElementModel { using InheritsFrom = std::tuple<GridProperties>; };
} // end namespace TTag

//! The grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FiniteElementModel>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using SV = GetPropType<TypeTag, Properties::SecondaryVariables>;
public:
    using type = FEGridVariables<GG, SV>;
};

//! Per default we use the base class of the secondary variables that
//! solely store the primary variables at an integration point
template<class TypeTag>
struct SecondaryVariables<TypeTag, TTag::FiniteElementModel>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using T = SecondaryVariablesBaseTraits<PV>;
public:
    using type = SecondaryVariablesBase<T>;
};

//! \todo TODO eliminate those? Can this have any functionality in FEM?
// //! We do not store the FVGeometry by default
// template<class TypeTag>
// struct EnableFVGridGeometryCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };
//
// //! We do not store the volume variables by default
// template<class TypeTag>
// struct EnableGridVolumeVariablesCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };
//
// //! disable flux variables data caching by default
// template<class TypeTag>
// struct EnableGridFluxVariablesCache<TypeTag, TTag::FiniteVolumeModel> { static constexpr bool value = false; };

//! Boundary types at a single degree of freedom
template<class TypeTag>
struct BoundaryTypes<TypeTag, TTag::FiniteElementModel> { using type = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; };

// TODO: bundle SolutionVector, JacobianMatrix
//       in LinearAlgebra traits

//! The type of a solution for the whole grid at a fixed time TODO: move to LinearAlgebra traits
template<class TypeTag>
struct SolutionVector<TypeTag, TTag::FiniteElementModel> { using type = Dune::BlockVector<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

//! Set the type of a global jacobian matrix from the solution types TODO: move to LinearAlgebra traits
template<class TypeTag>
struct JacobianMatrix<TypeTag, TTag::FiniteElementModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq() };
    using MatrixBlock = typename Dune::FieldMatrix<Scalar, numEq, numEq>;
public:
    using type = typename Dune::BCRSMatrix<MatrixBlock>;
};

//! Set the default for the ElementBoundaryTypes
template<class TypeTag>
struct ElementBoundaryTypes<TypeTag, TTag::FiniteElementModel> { using type = FEElementBoundaryTypes<GetPropType<TypeTag, Properties::BoundaryTypes>>; };

} // namespace Properties
} // namespace Dumux

 #endif
