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
 * \brief Helper class to get the required information on an interaction volume.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH

#include "methods.hh"
#include "facetypes.hh"

namespace Dumux
{

// Mpfa method-specific implementation of the helper class
template<class TypeTag, MpfaMethods Method, int dim>
class MpfaHelperBase {};

// dimension-specific implementation of the helper class (common for all methods)
template<class TypeTag, int dim>
class MpfaHelperImplementation {};

// Specialization for dim == 2
template<class TypeTag>
class MpfaHelperImplementation<TypeTag, 2> : public MpfaHelperBase<TypeTag,
                                                                   GET_PROP_VALUE(TypeTag, MpfaMethod),
                                                                   2>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<DimVector, dim>;

public:
    // Gets the common scv face in the outer element and the other scv face sharing the same vertex
    // orders them to form a right hand system, local indices can be deduced from on the rotational direction
    static ScvfVector getCommonAndNextScvFace(const SubControlVolumeFace& outsideScvf,
                                              const FVElementGeometry& fvGeometry,
                                              const bool clockWise)
    {
        LocalIndexType commonFaceIdx = clockWise ? 0 : 1;
        LocalIndexType nextFaceIdx = clockWise ? 1 : 0;

        ScvfVector scvfVector;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (outsideScvf.vertexIndex() == scvf.vertexIndex())
            {
                if (scvf.outsideScvIdx() == outsideScvf.insideScvIdx())
                    scvfVector[commonFaceIdx] = &scvf;
                else
                    scvfVector[nextFaceIdx] = &scvf;
            }
        }

        return scvfVector;
    }

    // calculates the inner normal vectors
    static LocalBasis calculateInnerNormals(const LocalBasis& localBasis)
    {
        static const Dune::FieldMatrix<Scalar, dim, dim> R = {{0.0, 1.0}, {-1.0, 0.0}};
        // make sure the basis forms a right hand system
        assert(isRightHandSystem(localBasis) > 0 && "Local basis does not form a right hand system");

        LocalBasis innerNormals;
        R.mv(localBasis[1], innerNormals[0]);
        R.mv(localBasis[0], innerNormals[1]);
        innerNormals[1] *= -1;

        return innerNormals;
    }

    // calculates the determinant of the local basis
    static Scalar calculateDetX(const LocalBasis& localBasis)
    {
        static const Dune::FieldMatrix<Scalar, dim, dim> R = {{0.0, 1.0}, {-1.0, 0.0}};
        // make sure the basis forms a right hand system
        assert(isRightHandSystem(localBasis) > 0 && "Local basis does not form a right hand system");

        DimVector tmp(0.0);
        R.mv(localBasis[1], tmp);

        return tmp*localBasis[0];
    }

    static bool isRightHandSystem(const LocalBasis& localBasis)
    {
        if (std::signbit(crossProduct<Scalar>(localBasis[0], localBasis[1])))
            return false;
        return true;
    }
};

// Specialization for dim == 3
template<class TypeTag>
class MpfaHelperImplementation<TypeTag, 3> : public MpfaHelperBase<TypeTag,
                                                                   GET_PROP_VALUE(TypeTag, MpfaMethod),
                                                                   3>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<DimVector, dim>;

public:

    // calculates the inner normal vectors
    static LocalBasis calculateInnerNormals(const LocalBasis& localBasis)
    {
        LocalBasis innerNormals;

        innerNormals[0] = crossProduct<Scalar>(localBasis[1], localBasis[2]);
        innerNormals[1] = crossProduct<Scalar>(localBasis[2], localBasis[0]);
        innerNormals[2] = crossProduct<Scalar>(localBasis[0], localBasis[1]);

        return innerNormals;
    }

    // calculates the determinant of the local basis
    static Scalar calculateDetX(const LocalBasis& localBasis)
    {
        assert(isRightHandSystem(localBasis) && "Local basis does not form a right hand system");
        return tripleProduct<Scalar>(localBasis[0], localBasis[1], localBasis[2]);
    }

    static bool isRightHandSystem(const LocalBasis& localBasis)
    {
        if (std::signbit(tripleProduct<Scalar>(localBasis[0], localBasis[1], localBasis[2])))
            return false;
        return true;
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specializations depending on the method and dimension are provided.
 */
template<class TypeTag>
class MpfaHelper : public MpfaHelperImplementation<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension>
{
    using ParentType = MpfaHelperImplementation<TypeTag, GET_PROP_TYPE(TypeTag, GridView)::dimension>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<DimVector, dim>;

public:
    // returns shared pointers to the two scv faces that share a vertex in the order of a right hand system
    static ScvfVector getScvFacesAtVertex(const GlobalIndexType vIdxGlobal,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry)
    {
        ScvfVector scvfVector;
        LocalBasis basisVectors;

        // The element center
        auto elementCenter = element.geometry().center();

        LocalIndexType count = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.vertexIndex() == vIdxGlobal)
            {
                scvfVector[count] = &scvf;
                basisVectors[count] = scvf.ipGlobal();
                basisVectors[count] -= elementCenter;
                count++;
            }
        }

        // We should always find dim faces sharing the vertex
        assert(count == dim);

        // sort the scv faces to form a right hand system
        if (!ParentType::isRightHandSystem(basisVectors))
            std::swap(scvfVector[0], scvfVector[1]);

        return scvfVector;
    }

    // Finds the local index in an ScvfVector that corresponds to the face that shares a facet with the given outsideScvf
    static LocalIndexType getCommonFaceLocalIndex(const SubControlVolumeFace& outsideScvf,
                                                  const ScvfVector& insideScvFaces)
    {
        for (int i = 0; i < insideScvFaces.size(); i++)
            if (insideScvFaces[i]->outsideScvIdx() == outsideScvf.insideScvIdx())
                return i;

        DUNE_THROW(Dune::InvalidStateException, "Could not find corresponding scvf in the provided vector of scvfs.");
    }

    // Returns the MpfaFaceType of an scv face
    static MpfaFaceTypes getMpfaFaceType(const Problem& problem,
                                         const Element& element,
                                         const SubControlVolumeFace& scvf)
    {
        if (!scvf.boundary())
            return MpfaFaceTypes::interior;

        auto bcTypes = problem.boundaryTypes(element, scvf);
        if (bcTypes.hasOnlyNeumann())
            return MpfaFaceTypes::neumann;
        if (bcTypes.hasOnlyDirichlet())
            return MpfaFaceTypes::dirichlet;

        // throw for outflow or mixed boundary conditions
        if (bcTypes.hasOutflow())
            DUNE_THROW(Dune::NotImplemented, "outflow BC for mpfa schemes");
        if (bcTypes.hasDirichlet() && bcTypes.hasNeumann())
            DUNE_THROW(Dune::InvalidStateException, "Mixed BC are not allowed for cellcentered schemes");

        DUNE_THROW(Dune::InvalidStateException, "unknown boundary condition type");
    }
};
} // end namespace

// The implemented helper classes need to be included here
#include <dumux/discretization/cellcentered/mpfa/omethod/helper.hh>
#include <dumux/discretization/cellcentered/mpfa/omethodfps/helper.hh>

#endif
