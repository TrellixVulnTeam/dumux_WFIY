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
 * \ingroup RichardsTests
 * \brief Spatial parameters for the Richards benchmarks
 */

#ifndef DUMUX_RICHARDS_BENCHMARKS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_BENCHMARKS_SPATIAL_PARAMETERS_HH

#include <dumux/common/parameters.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup RichardsTests
 * \brief Spatial parameters for the Richards benchmarks
 */
template<class GridGeometry, class Scalar>
class RichardsBenchmarkSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, RichardsBenchmarkSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RichardsBenchmarkSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsBenchmarkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("SpatialParams")
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , porosity_(getParam<Scalar>("SpatialParams.Porosity"))
    {
        // The potential rate decides about the type of the scenario.
        // See the problem file for more information.
        const auto potentialRate = getParam<Scalar>("Problem.SurfaceFluxMilliMeterPerDay");
        extrusionFactor_ = (potentialRate > 0) ? 0.1*0.1 : 0.05*0.05;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    /*!
     * \brief Returns the porosity [-] at a given location
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    /*!
     * \brief Returns the temperature [K] at a given location
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return 273.15 + 10.0; } // -> 10°C

    /*!
     * \brief Returns the extrusion factor [m^2] at a given location
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    { return extrusionFactor_; }

private:
    const PcKrSwCurve pcKrSwCurve_;
    const Scalar permeability_, porosity_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
