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
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for water management in PEM fuel cells.
 */

#ifndef DUMUX_FUELCELL_SPATIAL_PARAMS_HH
#define DUMUX_FUELCELL_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {
/*!
 * \ingroup TwoPNCTests
 * \brief Definition of the spatial parameters for water management in PEM fuel cells.
 */
template<class GridGeometry, class Scalar>
class FuelCellSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                       FuelCellSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ThisType = FuelCellSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    using PermeabilityType = DimWorldMatrix;

    FuelCellSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , K_(0)
    , pcKrSwCurve_("SpatialParams")
    {
        // intrinsic permeabilities
        K_[0][0] = 5e-11;
        K_[1][1] = 5e-11;
        temperature_ = getParam<Scalar>("Problem.InitialTemperature");
    }

    /*!
     * \brief Returns the hydraulic conductivity \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    DimWorldMatrix permeabilityAtPos(const GlobalPosition& globalPos) const
    { return K_; }

    /*!
     * \brief Defines the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The global position
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.2; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos A global coordinate vector
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The position of the center of the element
     * \return The wetting phase index
     * \note We set a hydrophobic material.
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::gasPhaseIdx; }

    /*!
     * \brief Returns the temperature in the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return temperature_;
    }

private:
    DimWorldMatrix K_;
    static constexpr Scalar eps_ = 1e-6;
    const PcKrSwCurve pcKrSwCurve_;
    Scalar temperature_;
};

} // end namespace Dumux

#endif
