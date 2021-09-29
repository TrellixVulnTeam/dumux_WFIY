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
 * \ingroup Common
 * \brief Base class for all finite volume problems
 */
#ifndef DUMUX_COMMON_FV_PROBLEM_HH
#define DUMUX_COMMON_FV_PROBLEM_HH

#include <memory>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/grid/common/gridenums.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/assembly/initialsolution.hh>

namespace Dumux {

/*!
 * \ingroup Common
 * \brief Base class for all finite-volume problems
 *
 * \note All quantities (regarding the units) are specified assuming a
 *       three-dimensional world. Problems discretized using 2D grids
 *       are assumed to be extruded by \f$1 m\f$ and 1D grids are assumed
 *       to have a cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class FVProblem
{
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { dim = GridView::dimension };

    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using PointSourceHelper = GetPropType<TypeTag, Properties::PointSourceHelper>;
    using PointSourceMap = std::map< std::pair<std::size_t, std::size_t>,
                                     std::vector<PointSource> >;

    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr bool isStaggered = GridGeometry::discMethod == DiscretizationMethod::staggered;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

public:
    //! Export spatial parameter type
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    //! export traits of this problem
    struct Traits
    {
        using Scalar = FVProblem::Scalar;
        using PrimaryVariables = FVProblem::PrimaryVariables;
        using NumEqVector = FVProblem::NumEqVector;
    };

    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams Spatially varying parameters
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    FVProblem(std::shared_ptr<const GridGeometry> gridGeometry,
              std::shared_ptr<SpatialParams> spatialParams,
              const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , spatialParams_(spatialParams)
    , paramGroup_(paramGroup)
    {
        // set a default name for the problem
        problemName_ = getParamFromGroup<std::string>(paramGroup, "Problem.Name");
    }

    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    FVProblem(std::shared_ptr<const GridGeometry> gridGeometry,
              const std::string& paramGroup = "")
    : FVProblem(gridGeometry,
                std::make_shared<SpatialParams>(gridGeometry),
                paramGroup)
    {}

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It could be either overwritten by the problem files, or simply
     * declared over the setName() function in the application file.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    /*!
     * \brief Set the problem name.
     *
     * This static method sets the simulation name, which should be
     * called before the application problem is declared! If not, the
     * default name "sim" will be used.
     *
     * \param newName The problem's name
     */
    void setName(const std::string& newName)
    {
        problemName_ = newName;
    }

    /*!
     * \name Boundary conditions and sources defining the problem
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub control volume
     */
    auto boundaryTypes(const Element &element,
                       const SubControlVolume &scv) const
    {
        if (!isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "boundaryTypes(..., scv) called for cell-centered method.");

        // forward it to the method which only takes the global coordinate
        return asImp_().boundaryTypesAtPos(scv.dofPosition());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    auto boundaryTypes(const Element &element,
                       const SubControlVolumeFace &scvf) const
    {
        if (isBox)
            DUNE_THROW(Dune::InvalidStateException,
                       "boundaryTypes(..., scvf) called for box method.");

        // forward it to the method which only takes the global coordinate
        return asImp_().boundaryTypesAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position of the finite volume in global coordinates
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any boundaryTypes method
        //! set Dirichlet boundary conditions everywhere for all primary variables
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face.
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     * \note used for cell-centered discretization schemes
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        // forward it to the method which only takes the global coordinate
        if (isBox)
        {
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scvf) called for box method.");
        }
        else
            return asImp_().dirichletAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param element The finite element
     * \param scv the sub control volume
     * \note used for cell-centered discretization schemes
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolume &scv) const
    {
        // forward it to the method which only takes the global coordinate
        if (!isBox && !isStaggered)
        {
            DUNE_THROW(Dune::InvalidStateException, "dirichlet(scv) called for other than box or staggered method.");
        }
        else
            return asImp_().dirichletAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the dirichlet condition ought to be
     *            set in global coordinates
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for Dirichlet conditions)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem specifies that some boundary "
                   "segments are dirichlet, but does not provide "
                   "a dirichlet() or a dirichletAtPos() method.");
    }

    /*!
     * \brief If internal Dirichlet contraints are enabled
     * Enables / disables internal (non-boundary) Dirichlet constraints. If this is overloaded
     * to return true, the assembler calls problem.hasInternalDirichletConstraint(element, scv).
     * This means you have to implement the following member function
     *
     *    bool hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const;
     *
     * which returns an indexable container of booleans defining for each equation if the corresponding dof associated
     * with the element/scv pair is constraint. If true is returned for a dof, the assembler calls
     * problem.internalDirichlet(element, scv). This means you have to additionally implement the following member function
     *
     *    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const;
     *
     * which returns the enforced Dirichlet values the dof associated with the element/scv pair.
     */
    static constexpr bool enableInternalDirichletConstraints()
    { return false; }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // forward it to the interface with only the global position
        return asImp_().neumannAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the boundary face's integration point in global coordinates
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any neumann method
        //! return no-flow Neumann boundary conditions at all Neumann boundaries
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(scv.center());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any source method
        //! return 0.0 (no source terms)
        return NumEqVector(0.0);
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the values method of the point source
     * has to return the absolute rate values in units
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const {}

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is
     * solution dependent
     *
     * \param source A single point source
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the values() method of the point sources returns
     * the absolute conserved quantity rate generated or annihilate in
     * units \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // forward to space dependent interface method
        asImp_().pointSourceAtPos(source, source.position());
    }

    /*!
     * \brief Evaluate the point sources (added by addPointSources)
     *        for all phases within a given sub-control-volume.
     *
     * This is the method for the case where the point source is space dependent
     *
     * \param pointSource A single point source
     * \param globalPos The point source position in global coordinates
     *
     * For this method, the \a values() method of the point sources returns
     * the absolute conserved quantity rate generated or annihilate in
     * units \f$ [ \textnormal{unit of conserved quantity} / s ] \f$. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    void pointSourceAtPos(PointSource& pointSource,
                          const GlobalPosition &globalPos) const {}

    /*!
     * \brief Add source term derivative to the Jacobian
     * \note Only needed in case of analytic differentiation and solution dependent sources
     */
    template<class MatrixBlock, class VolumeVariables>
    void addSourceDerivatives(MatrixBlock& block,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& volVars,
                              const SubControlVolume& scv) const {}

    /*!
     * \brief Adds contribution of point sources for a specific sub control volume
     *        to the values.
     *        Caution: Only overload this method in the implementation if you know
     *                 what you are doing.
     */
    template<class ElementVolumeVariables>
    NumEqVector scvPointSources(const Element &element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolume &scv) const
    {
        NumEqVector source(0);
        auto scvIdx = scv.indexInElement();
        auto key = std::make_pair(gridGeometry_->elementMapper().index(element), scvIdx);
        if (pointSourceMap_.count(key))
        {
            // Add the contributions to the dof source values
            // We divide by the volume. In the local residual this will be multiplied with the same
            // factor again. That's because the user specifies absolute values in kg/s.
            const auto volume = Extrusion::volume(scv)*elemVolVars[scv].extrusionFactor();

            for (const auto& ps : pointSourceMap_.at(key))
            {
                // we make a copy of the local point source here
                auto pointSource = ps;

                // Note: two concepts are implemented here. The PointSource property can be set to a
                // customized point source function achieving variable point sources,
                // see TimeDependentPointSource for an example. The second imitated the standard
                // dumux source interface with solDependentPointSource / pointSourceAtPos, methods
                // that can be overloaded in the actual problem class also achieving variable point sources.
                // The first one is more convenient for simple function like a time dependent source.
                // The second one might be more convenient for e.g. a solution dependent point source.

                // we do an update e.g. used for TimeDependentPointSource
                pointSource.update(asImp_(), element, fvGeometry, elemVolVars, scv);
                // call convienience problem interface function
                asImp_().pointSource(pointSource, element, fvGeometry, elemVolVars, scv);
                // at last take care about multiplying with the correct volume
                pointSource /= volume*pointSource.embeddings();
                // add the point source values to the local residual
                source += pointSource.values();
            }
        }

        return source;
    }

    /*!
     * \brief Compute the point source map, i.e. which scvs have point source contributions
     * \note Call this on the problem before assembly if you want to enable point sources set
     *       via the addPointSources member function.
     */
    void computePointSourceMap()
    {
        // clear the given point source maps in case it's not empty
        pointSourceMap_.clear();

        // get and apply point sources if any given in the problem
        std::vector<PointSource> sources;
        asImp_().addPointSources(sources);

        // if there are point sources calculate point source locations and save them in a map
        if (!sources.empty())
            PointSourceHelper::computePointSourceMap(*gridGeometry_, sources, pointSourceMap_, paramGroup());
    }

    /*!
     * \brief Get the point source map. It stores the point sources per scv
     */
    const PointSourceMap& pointSourceMap() const
    { return pointSourceMap_; }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     * \param sol the initial solution vector
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        assembleInitialSolution(sol, asImp_());
    }

    /*!
     * \brief Evaluate the initial value for
     * an element (for cell-centered models)
     * or vertex (for box / vertex-centered models)
     *
     * \param entity The dof entity (element or vertex)
     */
    template<class Entity>
    PrimaryVariables initial(const Entity& entity) const
    {
        static_assert(int(Entity::codimension) == 0 || int(Entity::codimension) == dim, "Entity must be element or vertex");
        return asImp_().initialAtPos(entity.geometry().center());
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        // Throw an exception (there is no reasonable default value
        // for initial values)
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "an initial() or an initialAtPos() method.");
    }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        // forward to generic interface
        return asImp_().extrusionFactorAtPos(scv.center());
    }

    /*!
     * \brief Return how much the domain is extruded at a given position.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    {
        // As a default, i.e. if the user's problem does not overload
        // any extrusion factor method, return 1.0
        return 1.0;
    }

    // \}

    //! The finite volume grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    //! The parameter group in which to retrieve runtime parameters
    const std::string& paramGroup() const
    { return paramGroup_; }

    //! Return the spatial parameters
    SpatialParams& spatialParams()
    { return *spatialParams_; }

    //! Return the spatial parameters
    const SpatialParams& spatialParams() const
    { return *spatialParams_; }

protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    //! The finite volume grid geometry
    std::shared_ptr<const GridGeometry> gridGeometry_;

    //! Spatially varying parameters
    std::shared_ptr<SpatialParams> spatialParams_;

    //! The parameter group in which to retrieve runtime parameters
    std::string paramGroup_;

    //! The name of the problem
    std::string problemName_;

    //! A map from an scv to a vector of point sources
    PointSourceMap pointSourceMap_;
};

} // end namespace Dumux

#endif
