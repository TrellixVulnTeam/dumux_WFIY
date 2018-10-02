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
/**
 * \file
 * \ingroup MultiDomain
 * \ingroup Geomechanics
 * \ingroup PoroElastic
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
#ifndef DUMUX_POROELASTIC_SUBPROBLEM_HH
#define DUMUX_POROELASTIC_SUBPROBLEM_HH

#include <dune/common/fmatrix.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/geomechanics/poroelastic/model.hh>
#include <dumux/geomechanics/fvproblem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>

#include "poroelasticspatialparams.hh"
#include "el2pco2tables.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class PoroElasticSubProblem;

namespace Properties {

NEW_TYPE_TAG(PoroElasticSubTypeTag, INHERITS_FROM(BoxModel, PoroElastic));
// Set the grid type
SET_TYPE_PROP(PoroElasticSubTypeTag, Grid, Dune::YaspGrid<3>);
// Set the problem property
SET_TYPE_PROP(PoroElasticSubTypeTag, Problem, Dumux::PoroElasticSubProblem<TypeTag>);

// Set the fluid system for TwoPSubProblem
SET_PROP(PoroElasticSubTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::BrineCO2<Scalar, El2P::CO2Tables>;
};

// The spatial parameters property
SET_TYPE_PROP(PoroElasticSubTypeTag, SpatialParams, PoroElasticSpatialParams< typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                              typename GET_PROP_TYPE(TypeTag, FVGridGeometry) >);

} // end namespace Properties

/*!
 * \ingroup MultiDomain
 * \ingroup Geomechanics
 * \ingroup PoroElastic
 *
 * \brief The poro-elastic sub-problem in the el2p coupled problem.
 */
template<class TypeTag>
class PoroElasticSubProblem : public GeomechanicsFVProblem<TypeTag>
{
    using ParentType = GeomechanicsFVProblem<TypeTag>;

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    //! The constructor
    PoroElasticSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup) {}

    //! The temperature in the domain
    static constexpr Scalar temperature()
    { return 273.15; }

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluate the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! Evaluate the boundary conditions for a Neumannboundary segment.
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Returns the effective fluid density
     */
    Scalar effectiveFluidDensity(const Element& element, const SubControlVolume& scv) const
    {
        // get context from coupling manager
        const auto& context = couplingManager().poroMechanicsCouplingContext();

        // here, we know that the flow problem uses cell-centered finite volumes, thus,
        // we simply take the volume variables of the scv (i.e. element) to obtain fluid properties
        const auto& facetVolVars = (*context.pmFlowElemVolVars)[scv.elementIndex()];
        Scalar wPhaseDensity = facetVolVars.density(FluidSystem::phase0Idx);
        Scalar nPhaseDensity = facetVolVars.density(FluidSystem::phase1Idx);
        Scalar Sw = facetVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = facetVolVars.saturation(FluidSystem::phase1Idx);
        return (wPhaseDensity * Sw + nPhaseDensity * Sn);
    }

    /*!
     * \brief Returns the effective pore pressure
     */
    template< class FluxVarsCache >
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    {
        // get context from coupling manager
        const auto& context = couplingManager().poroMechanicsCouplingContext();

        // here, we know that the flow problem uses cell-centered finite volumes, thus,
        // we simply take the volume variables of the element to obtain fluid properties
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto& facetVolVars = (*context.pmFlowElemVolVars)[eIdx];
        Scalar pw = facetVolVars.pressure(FluidSystem::phase0Idx);
        Scalar pn = facetVolVars.pressure(FluidSystem::phase1Idx);
        Scalar Sw = facetVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = facetVolVars.saturation(FluidSystem::phase1Idx);
        return (pw * Sw + pn * Sn);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<const CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    static constexpr Scalar eps_ = 3e-6;
};

} //end namespace

#endif
