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
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
#ifndef DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TRACER_TEST_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup TracerTests
 * \brief Definition of the spatial parameters for the tracer problem
 */
template<class TypeTag>
class TracerTestSpatialParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimWorld>;

public:

    TracerTestSpatialParams(const Problem& problem)
    : ParentType(problem) {}

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief Define the dispersivity.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     * \param elemSol The solution for all dofs of the element
     */
    template<class ElementSolution>
    Scalar dispersivity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return 0; }

    //! Fluid properties that are spatial params in the tracer model
    //! They can possible vary with space but are usually constants

    //! fluid density
    Scalar fluidDensity(const Element &element,
                        const SubControlVolume& scv) const
    { return 1000; }

    //! fluid molar mass
    Scalar fluidMolarMass(const Element &element,
                          const SubControlVolume& scv) const
    { return 18.0; }

    Scalar fluidMolarMass(const GlobalPosition &globalPos) const
    { return 18.0; }

    //! velocity field
    GlobalPosition velocity(const SubControlVolumeFace& scvf) const
    {
        GlobalPosition vel(1e-5);
        const auto globalPos = scvf.ipGlobal();
        const auto& x = globalPos[0];
        const auto& y = globalPos[1];

        vel[0] *= x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);
        vel[1] *= -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);

        return vel;
    }

    //! velocity field
    Scalar volumeFlux(const Element &element,
                      const FVElementGeometry& fvGeometry,
                      const ElementVolumeVariables& elemVolVars,
                      const SubControlVolumeFace& scvf) const
    {
        return velocity(scvf) * scvf.unitOuterNormal() * scvf.area()
               * elemVolVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
    }
};

} // end namespace Dumux

#endif
