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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \copydoc Dumux::FacetCouplingManager
 */
#ifndef DUMUX_FACETCOUPLING_MANAGER_HH
#define DUMUX_FACETCOUPLING_MANAGER_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux {

/*!
 * \brief Free function that allows the creation of a volume variables object
 *        interpolated to a given position within an element. This is the standard
 *        implementation which simply interpolates the solution to the given position
 *        and then performs a volume variables update with the interpolated solution.
 *
 * \note This assumes element-wise constant parameters for the computation of secondary
 *       variables. For heteregeneous parameter distributions a default implementation
 *       cannot be defined and an adequate overload of this function has to be provided!
 * \note For cell-centered schemes this is an unnecessary overhead because all variables
 *       are constant within the cells and a volume variables update can usually be realized
 *       more efficiently. This function is mainly to be used for the box scheme!
 */
template<class VolumeVariables, class Problem, class SolutionVector, class FVGeometry>
void makeInterpolatedVolVars(VolumeVariables& volVars,
                             const Problem& problem,
                             const SolutionVector& sol,
                             const FVGeometry& fvGeometry,
                             const typename FVGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                             const typename FVGeometry::FVGridGeometry::GridView::template Codim<0>::Entity::Geometry& elemGeom,
                             const typename FVGeometry::FVGridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate& pos)
{
    // interpolate solution and set it for each entry in element solution
    auto elemSol = elementSolution(element, sol, fvGeometry.fvGridGeometry());
    const auto centerSol = evalSolution(element, elemGeom, fvGeometry.fvGridGeometry(), elemSol, pos);
    for (int i = 0; i < fvGeometry.numScv(); ++i)
        elemSol[i] = centerSol;

    // Update volume variables with the interpolated solution. Note that this standard
    // implementation only works for element-wise constant parameters as we simply use
    // the first element scv for the vol var update. For heterogeneities within the element
    // or more complex models (e.g. 2p with interface solver) a corresponding overload
    // of this function has to be provided!
    volVars.update(elemSol, problem, element, *scvs(fvGeometry).begin());
}

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Implementation for the coupling manager between two domains of dimension d
 *        and (d-1) for models considering coupling across the bulk domain element facets.
 *        The implementations are specificto the discretization method used in the bulk
 *        domain, which is extracted automatically from the grid geometry corresponding
 *        to the provided id offset. Implementations for the different methods  have to be
 *        provided and included at the end of this file.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 * \tparam idOffset Offset added to the mapper-local domain ids for
 *                  the access to the grid quantities in grid creator
 * \tparam bulkDM Discretization method used in the bulk domain
 */
template< class MDTraits,
          class CouplingMapper,
          std::size_t idOffset = 0,
          DiscretizationMethod bulkDM = GET_PROP_TYPE(typename MDTraits::template SubDomainTypeTag<idOffset>, FVGridGeometry)::discMethod >
class FacetCouplingManagerImplementation;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that handles the coupling between an arbitrary number of sub domains
 *        in models where the coupling between d-dimensional  and (d-1)-dimensional
 *        domains occurs across the facets of the d-dimensional domain. The number
 *        of sub-domains involved is deduced from the provided multidomain traits
 *        class and specializations for the cases of two and three domains are provided below.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 */
template< class MDTraits,
          class CouplingMapper,
          std::size_t numSubDomains = MDTraits::numSubDomains >
class FacetCouplingManager;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that handles the coupling between two sub-domains in models where
 *        the coupling between the two occurs across the facets of the higher-
 *        dimensional domain. Here, we simply inherit from the discretization
 *        method-specific implementation.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 */
template< class MDTraits, class CouplingMapper >
class FacetCouplingManager<MDTraits, CouplingMapper, /*numSubDomains*/2>
: public FacetCouplingManagerImplementation<MDTraits, CouplingMapper> {};

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that handles the coupling between three sub-domains in models where
 *        the coupling between the two occurs across the facets of the d- and (d-1)-
 *        dimensional domains.
 *
 * \tparam MDTraits The multidomain traits containing the types on all sub-domains
 * \tparam CouplingMapper Class containing maps on the coupling between dofs of different grids
 */
template< class MDTraits, class CouplingMapper >
class FacetCouplingManager<MDTraits, CouplingMapper, /*numSubDomains*/3>
: public FacetCouplingManagerImplementation<MDTraits, CouplingMapper>
, public FacetCouplingManagerImplementation<MDTraits, CouplingMapper, /*idOffset*/1>
{
    using BulkFacetManager = FacetCouplingManagerImplementation<MDTraits, CouplingMapper>;
    using FacetEdgeManager = FacetCouplingManagerImplementation<MDTraits, CouplingMapper, /*idOffset*/1>;

    // convenience aliases and instances of the domain ids
    using BulkIdType = typename MDTraits::template DomainIdx<0>;
    using FacetIdType = typename MDTraits::template DomainIdx<1>;
    using EdgeIdType = typename MDTraits::template DomainIdx<2>;
    static constexpr auto bulkId = BulkIdType();
    static constexpr auto facetId = FacetIdType();
    static constexpr auto edgeId = EdgeIdType();

    // the sub-domain type tags
    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    // further types specific to the sub-problems
    template<std::size_t id> using LocalResidual = typename GET_PROP_TYPE(SubDomainTypeTag<id>, LocalResidual);
    template<std::size_t id> using PrimaryVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, PrimaryVariables);
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using GridView = typename FVGridGeometry<id>::GridView;
    template<std::size_t id> using IndexType = typename GridView<id>::IndexSet::IndexType;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVariables);
    template<std::size_t id> using ElementVolumeVariables = typename GridVariables<id>::GridVolumeVariables::LocalView;
    template<std::size_t id> using ElementFluxVariablesCache = typename GridVariables<id>::GridFluxVariablesCache::LocalView;

public:

    //! types used for coupling stencils
    template<std::size_t i, std::size_t j>
    using CouplingStencilType = typename std::conditional< (j > 1),
                                                            typename FacetEdgeManager::template CouplingStencilType<i, j>,
                                                            typename BulkFacetManager::template CouplingStencilType<i, j> >::type;

    //! the type of the solution vector
    using SolutionVector = typename MDTraits::SolutionVector;

    /*!
     * \brief Initialize the coupling manager.
     *
     * \param bulkProblem The problem to be solved on the (3d) bulk domain
     * \param facetProblem The problem to be solved on the (2d) facet domain
     * \param edgeProblem The problem to be solved on the (1d) edge domain
     * \param couplingMapper The mapper object containing the connectivity between the domains
     * \tparam curSol The current solution
     */
    void init(std::shared_ptr< Problem<bulkId> > bulkProblem,
              std::shared_ptr< Problem<facetId> > facetProblem,
              std::shared_ptr< Problem<edgeId> > edgeProblem,
              std::shared_ptr< CouplingMapper > couplingMapper,
              const SolutionVector& curSol)
    {
        BulkFacetManager::init(bulkProblem, facetProblem, couplingMapper, curSol);
        FacetEdgeManager::init(facetProblem, edgeProblem, couplingMapper, curSol);
    }

    //! Pull up functionalities from the parent classes
    using BulkFacetManager::couplingStencil;
    using FacetEdgeManager::couplingStencil;

    using BulkFacetManager::isCoupled;
    using FacetEdgeManager::isCoupled;

    using BulkFacetManager::getLowDimVolVars;
    using FacetEdgeManager::getLowDimVolVars;

    using BulkFacetManager::evalSourcesFromBulk;
    using FacetEdgeManager::evalSourcesFromBulk;

    using BulkFacetManager::evalCouplingResidual;
    using FacetEdgeManager::evalCouplingResidual;

    using BulkFacetManager::bindCouplingContext;
    using FacetEdgeManager::bindCouplingContext;

    using BulkFacetManager::updateCouplingContext;
    using FacetEdgeManager::updateCouplingContext;

    using BulkFacetManager::updateCoupledVariables;
    using FacetEdgeManager::updateCoupledVariables;

    /*!
     * \brief The coupling stencil of the bulk with the edge domain (empty stencil).
     */
    const CouplingStencilType<bulkId, edgeId>& couplingStencil(BulkIdType domainI,
                                                               const Element<bulkId>& element,
                                                               EdgeIdType domainJ) const
    { return FacetEdgeManager::getEmptyStencil(edgeId); }

    /*!
     * \brief The coupling stencil of the edge with the bulk domain (empty stencil).
     */
    const CouplingStencilType<edgeId, bulkId>& couplingStencil(EdgeIdType domainI,
                                                               const Element<edgeId>& element,
                                                               BulkIdType domainJ) const
    { return BulkFacetManager::getEmptyStencil(bulkId); }

    /*!
     * \brief updates the current solution. We have to overload this here
     *        to avoid ambiguity and update the solution in both managers
     */
    void updateSolution(const SolutionVector& sol)
    {
        BulkFacetManager::updateSolution(sol);
        FacetEdgeManager::updateSolution(sol);
    }

    /*!
     * \brief Interface for evaluating the coupling residual between the bulk and the edge domain.
     *        This is always zero as coupling only occurs between grids of codimension one. These
     *        overloads are provided by the two parent classes. However, we need this overload in
     *        order for the overload resolution not to fail.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || ((i==edgeId && j==bulkId))), int> = 0>
    typename LocalResidual<i>::ElementResidualVector
    evalCouplingResidual(Dune::index_constant<i> domainI,
                         const LocalAssembler& localAssembler,
                         Dune::index_constant<j> domainJ,
                         IndexType<j> dofIdxGlobalJ)
    {
        typename LocalResidual<i>::ElementResidualVector res(1);
        res = 0.0;
        return res;
    }

    /*!
     * \brief Interface for binding the coupling context for the facet domain. In this case
     *        we have to bind both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class Assembler >
    void bindCouplingContext(FacetIdType, const Element<facetId>& element, const Assembler& assembler)
    {
        BulkFacetManager::bindCouplingContext(facetId, element, assembler);
        FacetEdgeManager::bindCouplingContext(facetId, element, assembler);
    }

    /*!
     * \brief Interface for updating the coupling context of the facet domain. In this case
     *        we have to update both the facet -> bulk and the facet -> edge coupling context.
     */
    template< class FacetLocalAssembler >
    void updateCouplingContext(FacetIdType domainI,
                               const FacetLocalAssembler& facetLocalAssembler,
                               FacetIdType domainJ,
                               IndexType<facetId> dofIdxGlobalJ,
                               const PrimaryVariables<facetId>& priVarsJ,
                               unsigned int pvIdxJ)
    {
        BulkFacetManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
        FacetEdgeManager::updateCouplingContext(domainI, facetLocalAssembler, domainJ, dofIdxGlobalJ, priVarsJ, pvIdxJ);
    }

    /*!
     * \brief Interface for updating the coupling context between the bulk and the edge domain.
     *        We do nothing here because coupling only occurs between grids of codimension one.
     *        We have to provide this overload as the overload resolution fails otherwise.
     */
    template<std::size_t i,
             std::size_t j,
             class LocalAssembler,
             std::enable_if_t<((i==bulkId && j==edgeId) || (i==edgeId && j==bulkId)), int> = 0>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssembler& localAssembler,
                               Dune::index_constant<j> domainJ,
                               IndexType<j> dofIdxGlobalJ,
                               const PrimaryVariables<j>& priVarsJ,
                               unsigned int pvIdxJ)
    { /*do nothing here*/ }

    /*!
     * \brief Interface for updating the local views of the facet domain after updateCouplingContext
     *        the coupling context. In this case we have to forward the both managers as the facet
     *        domain is a part in both.
     */
    template< class FacetLocalAssembler, class UpdatableElementVolVars, class UpdatableFluxVarCache>
    void updateCoupledVariables(FacetIdType domainI,
                                const FacetLocalAssembler& facetLocalAssembler,
                                UpdatableElementVolVars& elemVolVars,
                                UpdatableFluxVarCache& elemFluxVarsCache)
    {
        BulkFacetManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
        FacetEdgeManager::updateCoupledVariables(domainI, facetLocalAssembler, elemVolVars, elemFluxVarsCache);
    }

    //! Return a const reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    const Problem<id>& problem() const { return BulkFacetManager::template problem<id>(); }

    //! Return a reference to bulk or facet problem
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == facetId), int> = 0>
    Problem<id>& problem() { return BulkFacetManager::template problem<id>(); }

    //! Return a const reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    const Problem<id>& problem() const { return FacetEdgeManager::template problem<id>(); }

    //! Return a reference to edge problem
    template<std::size_t id, std::enable_if_t<(id == edgeId), int> = 0>
    Problem<id>& problem() { return FacetEdgeManager::template problem<id>(); }
};

} // end namespace Dumux

// Here, we have to include all available implementations
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmanager.hh>

#endif
