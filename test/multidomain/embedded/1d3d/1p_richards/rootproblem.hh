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
 *
 * \brief A test problem for the one-phase root model:
 * Sap is flowing through a 1d network root xylem.
 */
#ifndef DUMUX_ROOT_PROBLEM_HH
#define DUMUX_ROOT_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "rootspatialparams.hh"

namespace Dumux {
// forward declaration
template <class TypeTag> class RootProblem;

namespace Properties {

NEW_TYPE_TAG(RootTypeTag, INHERITS_FROM(CCTpfaModel, OneP));

// Set the grid type
SET_TYPE_PROP(RootTypeTag, Grid, Dune::FoamGrid<1, 3>);

SET_BOOL_PROP(RootTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(RootTypeTag, SolutionDependentAdvection, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentMolecularDiffusion, false);
SET_BOOL_PROP(RootTypeTag, SolutionDependentHeatConduction, false);

// Set the problem property
SET_TYPE_PROP(RootTypeTag, Problem, RootProblem<TypeTag>);

// the fluid system
SET_PROP(RootTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Set the problem property
SET_TYPE_PROP(RootTypeTag, LocalResidual, OnePIncompressibleLocalResidual<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(RootTypeTag, SpatialParams, RootSpatialParams<typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                            typename GET_PROP_TYPE(TypeTag, Scalar)>);

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief Exact solution 1D-3D
 */
template <class TypeTag>
class RootProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PointSource = typename GET_PROP_TYPE(TypeTag, PointSource);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename FVGridGeometry::GlobalCoordinate;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:

    template<class SpatialParams>
    RootProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<SpatialParams> spatialParams,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, spatialParams)
    , couplingManager_(couplingManager)
    {
        //read parameters from input file
        name_ = getParam<std::string>("Problem.Name") + "_1d";
        transpirationRate_ = getParam<Scalar>("BoundaryConditions.TranspirationRate");
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto radius = this->spatialParams().radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return 273.15 + 10.0; }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The global position
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    template<class ElementVolumeVariables>
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolvars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        if (scvf.center()[2] + eps_ > this->fvGridGeometry().bBoxMax()[2])
        {
            const auto r = this->spatialParams().radius(scvf.insideScvIdx());
            values[Indices::conti0EqIdx] = transpirationRate_/(M_PI*r*r)/scvf.area();
        }
        return values;

    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

     /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param pointSource A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub-control volume within the element
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute rate mass generated or annihilate in kg/s. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];

        const auto lowDimElementIdx = this->couplingManager().pointSourceData(source.id()).lowDimElementIdx();
        const Scalar Kr = this->spatialParams().Kr(lowDimElementIdx);
        const Scalar rootRadius = this->spatialParams().radius(lowDimElementIdx);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto density = 1000;
        const Scalar sourceValue = 2* M_PI *rootRadius * Kr *(pressure3D - pressure1D)*density;
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    // \}

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (root): " << source << " (kg/s) / "
                  <<                           source*3600*24*1000 << " (g/day)" << '\n';
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(this->spatialParams().getRadii(), "radius");
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar transpirationRate_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} //end namespace Dumux

#endif
