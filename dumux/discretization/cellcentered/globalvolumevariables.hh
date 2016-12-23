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
 * \brief The global volume variables class for cell centered models
 */
#ifndef DUMUX_DISCRETIZATION_CC_GLOBAL_VOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_CC_GLOBAL_VOLUMEVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the volume variables vector
 */
template<class TypeTag, bool enableGlobalVolVarsCache>
class CCGlobalVolumeVariables
{};

//! specialization in case of storing the volume variables
template<class TypeTag>
class CCGlobalVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/true>
{
    // The local class needs to access and change volVars
    friend typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    // The local jacobian needs to access and change volVars for derivative calculation
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    // as does the primary variable switch
    friend class PrimaryVariableSwitch<TypeTag>;
    friend typename GET_PROP_TYPE(TypeTag, PrimaryVariableSwitch);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementSolution = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    void update(Problem& problem, const SolutionVector& sol)
    {
        problemPtr_ = &problem;

        auto numScv = problem.model().globalFvGeometry().numScv();
        auto numBoundaryScvf = problem.model().globalFvGeometry().numBoundaryScvf();

        volumeVariables_.resize(numScv + numBoundaryScvf);
        for (const auto& element : elements(problem.gridView()))
        {
            auto elementSol = problem.model().elementSolution(element, sol);
            auto fvGeometry = localView(problem.model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
                volumeVariables_[scv.index()].update(elementSol,
                                                     problem,
                                                     element,
                                                     scv);

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto bcTypes = problem.boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    const auto insideScvIdx = scvf.insideScvIdx();
                    const auto& insideScv = fvGeometry.scv(insideScvIdx);
                    const auto dirichletPriVars = problem.dirichlet(element, scvf);

                    volumeVariables_[scvf.outsideScvIdx()].update(ElementSolution({dirichletPriVars}),
                                                                  problem,
                                                                  element,
                                                                  insideScv);
                }
            }
        }
    }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const CCGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

    const VolumeVariables& volVars(const IndexType scvIdx) const
    { return volumeVariables_[scvIdx]; }

    VolumeVariables& volVars(const IndexType scvIdx)
    { return volumeVariables_[scvIdx]; }
private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<VolumeVariables> volumeVariables_;
};


//! Specialization when the current volume variables are not stored globally
template<class TypeTag>
class CCGlobalVolumeVariables<TypeTag, /*enableGlobalVolVarsCache*/false>
{
    // local class needs access to the problem
    friend typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    void update(Problem& problem, const SolutionVector& sol)
    { problemPtr_ = &problem; }

    /*!
     * \brief Return a local restriction of this global object
     *        The local object is only functional after calling its bind/bindElement method
     *        This is a free function that will be found by means of ADL
     */
    friend inline ElementVolumeVariables localView(const CCGlobalVolumeVariables& global)
    { return ElementVolumeVariables(global); }

private:
    Problem& problem_() const
    { return *problemPtr_;}

    Problem* problemPtr_;
};

} // end namespace

#endif
