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
 * \ingroup FacetTests
 * \brief The spatial parameters for the compositional single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the compositional single-phase facet coupling test.
 */
template< class GridGeometry, class Scalar >
class OnePSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP< GridGeometry, Scalar, OnePSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = OnePSpatialParams< GridGeometry, Scalar >;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP< GridGeometry, Scalar, ThisType >;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type used for permeabilities
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        extrusion_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture", 1.0);
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    //! Returns the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    //! Returns the extrusion factor
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return extrusion_; }

    //! Returns the temperature at the given position
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

private:
    PermeabilityType permeability_;
    Scalar porosity_;
    Scalar extrusion_;
};

} // end namespace Dumux

#endif
