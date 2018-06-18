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
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
 * \brief The interface of the coupling manager for multi domain problems
 */

#ifndef DUMUX_STAGGERED_COUPLING_MANAGER_HH
#define DUMUX_STAGGERED_COUPLING_MANAGER_HH

#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
 * \brief Base coupling manager for the staggered discretization.
 */
template<class MDTraits, class Implementation>
class StaggeredCouplingManagerBase: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using Problem = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);

    using StaggeredSubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<0>;
    using GridView = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using LocalResidual = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, LocalResidual);

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(StaggeredSubDomainTypeTag, FacePrimaryVariables);

    using CouplingStencils = std::unordered_map<std::size_t, std::vector<std::size_t> >;
    using CouplingStencil = CouplingStencils::mapped_type;

public:

    using Traits = MDTraits;

    static constexpr auto cellCenterIdx = Dune::index_constant<0>();
    static constexpr auto faceIdx = Dune::index_constant<1>();

    void init(std::shared_ptr<const Problem<0>> problem)
    {
        problemTuple_ = std::make_tuple(problem, problem);
    }

    template<class... Args>
    void init(Args&&... args)
    {
        problemTuple_ = std::make_tuple(args...);
    }

    void init(typename Traits::ProblemTuple&& problemTuple)
    {
        problemTuple_ = std::move(problemTuple);
    }


    /*!
     * \copydoc ParentType::updateCouplingContext
     *
     * \note this is a specialization for updating the cellcentered dof
     *       residual context w.r.t. staggered face dof changes
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<cellCenterIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<faceIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const FacePrimaryVariables& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \copydoc ParentType::updateCouplingContext
     *
     * \note this is a specialization for updating the staggered face dof
     *       residual context w.r.t. cellcentered dof changes
     */
    template<class LocalAssemblerI>
    void updateCouplingContext(Dune::index_constant<faceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<cellCenterIdx> domainJ,
                               std::size_t dofIdxGlobalJ,
                               const CellCenterPrimaryVariables& priVarsJ,
                               int pvIdxJ)
    {
        this->curSol()[domainJ][dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param domainJ the domain index of domain j

     * \note  this is a specialization for getting the indices of the coupled staggered face dofs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<cellCenterIdx> domainI,
                                           const Element& elementI,
                                           Dune::index_constant<faceIdx> domainJ) const
    {
        const auto& connectivityMap = problem(domainI).fvGridGeometry().connectivityMap();
        const auto eIdx = problem(domainI).fvGridGeometry().elementMapper().index(elementI);
        return connectivityMap(domainI, domainJ, eIdx);
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given subcontrolvolume face of domain i
     *
     * \param domainI the domain index of domain i
     * \param scvfI the coupled subcontrolvolume face of domain í
     * \param domainJ the domain index of domain j

     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil couplingStencil(Dune::index_constant<i> domainI,
                                           const SubControlVolumeFace& scvfI,
                                           Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "Domain i cannot be coupled to itself!");
        static_assert(AlwaysFalse<Dune::index_constant<i>>::value,
                      "The coupling manager does not implement the couplingStencil() function" );

        return CouplingStencil(); // supress compiler warning of function having no return statement
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given subcontrolvolume face of domain i
     *
     * \param domainI the domain index of domain i
     * \param scvfI the coupled subcontrolvolume face of domain í
     * \param domainJ the domain index of domain j

     * \note this is a specialization for getting the indices of the coupled cellcentered dofs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<faceIdx> domainI,
                                           const SubControlVolumeFace& scvfI,
                                           Dune::index_constant<cellCenterIdx> domainJ) const
    {
        const auto& connectivityMap = problem(domainI).fvGridGeometry().connectivityMap();
        return connectivityMap(domainI, domainJ, scvfI.index());
    }

    /*!
     * \copydoc ParentType::evalCouplingResidual
     *
     * \note this is a specialization for calculating the coupled residual for cellcentered dofs
     *       w.r.t. staggered face dof changes
     */
    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<cellCenterIdx> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        static_assert(domainI != domainJ, "Domain i cannot be coupled to itself!");
        return localAssemblerI.evalLocalResidualForCellCenter();
    }

     /*!
      * \brief evaluates the face residual of a coupled face of domain i which depends on the variables
      *        at the degree of freedom with index dofIdxGlobalJ of domain j
      *
      * \param domainI the domain index of domain i
      * \param scvfI the subcontrol volume face whose residual shall be evaluated of domain i
      * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
      * \param domainJ the domain index of domain j
      * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
      *
      * \note  the default implementation evaluates the complete face residual
      *        if only certain terms of the residual of the residual are coupled
      *        to dof with index dofIdxGlobalJ the function can be overloaded in the coupling manager
      * \return the face residual
      */
    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<faceIdx> domainI,
                                        const SubControlVolumeFace& scvfI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        static_assert(domainI != domainJ, "Domain i cannot be coupled to itself!");
        return localAssemblerI.evalLocalResidualForFace(scvfI);
    }

    //! Return a reference to the problem
    template<std::size_t id>
    const Problem<id>& problem(Dune::index_constant<id> domainIdx) const
    {
        assert(std::get<id>(problemTuple_) && "No problem set. Call init() first!");
        return *std::get<id>(problemTuple_);
    }

private:
    typename Traits::ProblemTuple problemTuple_;

};

template<class MDTraits>
class StaggeredCouplingManager : public StaggeredCouplingManagerBase<MDTraits, StaggeredCouplingManager<MDTraits>>
{
    using ParentType = StaggeredCouplingManagerBase<MDTraits, StaggeredCouplingManager<MDTraits>>;

public:
    using ParentType::ParentType;

};

} //end namespace Dumux

#endif
