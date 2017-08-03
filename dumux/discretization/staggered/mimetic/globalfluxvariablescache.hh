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
 * \brief The global object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_MIMETIC_GLOBAL_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_MIMETIC_GLOBAL_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/staggered/mimetic/fluxvariablescachefiller.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class MimeticGlobalFluxVariablesCache;

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class MimeticGlobalFluxVariablesCache<TypeTag, true>
{
    // the local class needs access to the problem
    friend StaggeredElementFluxVariablesCache<TypeTag, true>;
    // the filler class needs access to the access operators
    friend MimeticFluxVariablesCacheFiller<TypeTag>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCacheFiller = MimeticFluxVariablesCacheFiller<TypeTag>;

public:
    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    void update(Problem& problem)
    {
        problemPtr_ = &problem;

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(problem);

        const auto& globalFvGeometry = problem.model().globalFvGeometry();
        fluxVarsCache_.resize(globalFvGeometry.numScvf());
        for (const auto& element : elements(problem.gridView()))
        {
            // Prepare the geometries within the elements of the stencil
            auto fvGeometry = localView(globalFvGeometry);
            fvGeometry.bind(element);

            auto elemVolVars = localView(problem.model().curGlobalVolVars());
            elemVolVars.bind(element, fvGeometry, problem.model().curSol());

            for (auto&& scvf : scvfs(fvGeometry))
            {
                filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf);
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const MimeticGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:
    // access operators in the case of caching
    const FluxVariablesCache& operator [](IndexType scvfIdx) const
    { return fluxVarsCache_[scvfIdx]; }

    FluxVariablesCache& operator [](IndexType scvfIdx)
    { return fluxVarsCache_[scvfIdx]; }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class MimeticGlobalFluxVariablesCache<TypeTag, false>
{
    // the local class needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

public:
    // When global flux variables caching is disabled, we don't need to update the cache
    void update(Problem& problem)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementFluxVariablesCache localView(const MimeticGlobalFluxVariablesCache& global)
    { return ElementFluxVariablesCache(global); }

private:

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;
};

} // end namespace

#endif
