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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FaceCenteredStaggeredElementVolumeVariables
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_ELEMENTVOLUMEVARIABLES_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_ELEMENTVOLUMEVARIABLES_HH

#include <algorithm>
#include <cassert>
#include <vector>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the face variables vector
 */
template<class GridVolumeVariables, bool cachingEnabled>
class FaceCenteredStaggeredElementVolumeVariables
{};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of storing the face variables globally.
 */
template<class GVV>
class FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/true>
{
public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    FaceCenteredStaggeredElementVolumeVariables(const GridVolumeVariables& gridVolumeVariables)
    : gridVolumeVariablesPtr_(&gridVolumeVariables)
    , numScv_(gridVolumeVariables.problem().gridGeometry().numScv())
    {}

    //! operator for the access with an scvf
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    {
        if (scv.index() < numScv_)
            return gridVolVars().volVars(scv.index());
        else
            return boundaryVolumeVariables_[getLocalIdx_(scv.index())];
    }

    //! operator for the access with an index
    //! needed for cc methods for the access to the boundary volume variables
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    {
        if (scvIdx < numScv_)
            return gridVolVars().volVars(scvIdx);
        else
            return boundaryVolumeVariables_[getLocalIdx_(scvIdx)];
    }


    //! For compatibility reasons with the case of not storing the face vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        if (!fvGeometry.hasBoundaryScvf())
            return;

        clear_();
        boundaryVolVarIndices_.reserve(fvGeometry.gridGeometry().numBoundaryScv());
        boundaryVolumeVariables_.reserve(fvGeometry.gridGeometry().numBoundaryScv());

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() || scvf.isFrontal())
                continue;

            // check if boundary is a pure dirichlet boundary
            const auto& problem = gridVolVars().problem();
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            auto addBoundaryVolVars = [&](const auto& scvFace)
            {
                const auto& scvI = fvGeometry.scv(scvFace.insideScvIdx());
                typename VolumeVariables::PrimaryVariables pv(problem.dirichlet(element, scvFace)[scvI.directionIndex()]);
                const auto dirichletPriVars = elementSolution<FVElementGeometry>(pv);

                VolumeVariables volVars;
                volVars.update(dirichletPriVars,
                               problem,
                               element,
                               scvI);

                boundaryVolumeVariables_.emplace_back(std::move(volVars));
                boundaryVolVarIndices_.push_back(scvFace.outsideScvIdx());
            };

            if (bcTypes.hasOnlyDirichlet())
            {
                addBoundaryVolVars(scvf);
                continue;
            }

            // treat domain corners
            if (const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf); orthogonalScvf.boundary())
            {
                if (const auto orthogonalBcTypes = problem.boundaryTypes(element, orthogonalScvf); orthogonalBcTypes.hasOnlyDirichlet())
                    addBoundaryVolVars(scvf);
            }
        }
        assert(boundaryVolumeVariables_.size() == boundaryVolVarIndices_.size());
    }

    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {}


    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolumeVariablesPtr_; }


private:
    //! Clear all local storage
    void clear_()
    {
        boundaryVolVarIndices_.clear();
        boundaryVolumeVariables_.clear();
    }

    //! map a global scv index to the local storage index
    int getLocalIdx_(const std::size_t volVarIdx) const
    {
        const auto it = std::find(boundaryVolVarIndices_.begin(), boundaryVolVarIndices_.end(), volVarIdx);
        assert(it != boundaryVolVarIndices_.end() && "Could not find the current volume variables for volVarIdx!");
        return std::distance(boundaryVolVarIndices_.begin(), it);
    }

    const GridVolumeVariables* gridVolumeVariablesPtr_;
    const std::size_t numScv_;
    std::vector<std::size_t> boundaryVolVarIndices_;
    std::vector<VolumeVariables> boundaryVolumeVariables_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the face variables vector. Specialization for the case of not storing the face variables globally.
 */
template<class GVV>
class FaceCenteredStaggeredElementVolumeVariables<GVV, /*cachingEnabled*/false>
{
    using GridGeometry = std::decay_t<decltype(std::declval<GVV>().problem().gridGeometry())>;
    static constexpr auto dim = GridGeometry::GridView::dimension;
    static constexpr auto numInsideVolVars = dim * 2;
    static constexpr auto numOutsideVolVars = numInsideVolVars * 2 * (dim - 1);

public:
    //! export type of the grid volume variables
    using GridVolumeVariables = GVV;

    //! export type of the volume variables
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    FaceCenteredStaggeredElementVolumeVariables(const GridVolumeVariables& globalFacesVars) : gridVolumeVariablesPtr_(&globalFacesVars) {}

    //! const operator for the access with an scvf
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    const VolumeVariables& operator [](const SubControlVolume& scv) const
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    //! const operator for the access with an index
    const VolumeVariables& operator [](const std::size_t scvIdx) const
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! operator for the access with an scvf
    template<class SubControlVolume, typename std::enable_if_t<!std::is_integral<SubControlVolume>::value, int> = 0>
    VolumeVariables& operator [](const SubControlVolume& scv)
    { return volumeVariables_[getLocalIdx_(scv.index())]; }

    // operator for the access with an index
    VolumeVariables& operator [](const std::size_t scvIdx)
    { return volumeVariables_[getLocalIdx_(scvIdx)]; }

    //! For compatibility reasons with the case of not storing the vol vars.
    //! function to be called before assembling an element, preparing the vol vars within the stencil
    template<class FVElementGeometry, class SolutionVector>
    void bind(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
              const FVElementGeometry& fvGeometry,
              const SolutionVector& sol)
    {
        clear_();

        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();

        volVarIndices_.reserve(numInsideVolVars + numInsideVolVars);
        volumeVariables_.reserve(numInsideVolVars + numInsideVolVars);

        for (const auto& scv : scvs(fvGeometry))
        {
            for (const auto otherScvIdx : gridGeometry.connectivityMap()[scv.index()])
            {
                if (!volVarsInserted_(otherScvIdx))
                {
                    const auto& otherScv = fvGeometry.scv(otherScvIdx);
                    volVarIndices_.push_back(otherScvIdx);
                    volumeVariables_.emplace_back();
                    const auto& otherElement = gridGeometry.element(otherScv.elementIndex());
                    volumeVariables_.back().update(elementSolution(otherElement, sol, gridGeometry),
                                                   problem,
                                                   otherElement,
                                                   otherScv);
                }
            }
        }

        if (fvGeometry.hasBoundaryScvf())
        {
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (!scvf.boundary() || scvf.isFrontal())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto& problem = gridVolVars().problem();
                const auto bcTypes = problem.boundaryTypes(element, scvf);

                auto addBoundaryVolVars = [&](const auto& scvFace)
                {
                    const auto& scvI = fvGeometry.scv(scvFace.insideScvIdx());
                    typename VolumeVariables::PrimaryVariables pv(problem.dirichlet(element, scvFace)[scvI.directionIndex()]);
                    const auto dirichletPriVars = elementSolution<FVElementGeometry>(pv);

                    VolumeVariables volVars;
                    volVars.update(dirichletPriVars,
                                   problem,
                                   element,
                                   scvI);

                    volumeVariables_.emplace_back(std::move(volVars));
                    volVarIndices_.push_back(scvFace.outsideScvIdx());
                };

                if (bcTypes.hasOnlyDirichlet())
                {
                    addBoundaryVolVars(scvf);
                    continue;
                }

                // treat domain corners
                if (const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf); orthogonalScvf.boundary())
                {
                    if (const auto orthogonalBcTypes = problem.boundaryTypes(element, orthogonalScvf); orthogonalBcTypes.hasOnlyDirichlet())
                        addBoundaryVolVars(scvf);
                }
            }
        }
    }

    //! Binding of an element, prepares only the face variables of the element
    //! specialization for Staggered models
    template<class FVElementGeometry, class SolutionVector>
    void bindElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                     const FVElementGeometry& fvGeometry,
                     const SolutionVector& sol)
    {
        clear_();
        const auto& problem = gridVolVars().problem();
        const auto& gridGeometry = fvGeometry.gridGeometry();
        volVarIndices_.reserve(numInsideVolVars);

        for (const auto& scv : scvs(fvGeometry))
        {
            volVarIndices_.push_back(scv.index());
            volumeVariables_.emplace_back();
            volumeVariables_.back().update(elementSolution(element, sol, gridGeometry),
                                           problem,
                                           element,
                                           scv);
        }
    }

    //! The global volume variables object we are a restriction of
    const GridVolumeVariables& gridVolVars() const
    { return *gridVolumeVariablesPtr_; }

private:

    //! Clear all local storage
    void clear_()
    {
        volVarIndices_.clear();
        volumeVariables_.clear();
    }

    bool volVarsInserted_(const std::size_t scvIdx) const
    {
        return std::find(volVarIndices_.begin(), volVarIndices_.end(), scvIdx) != volVarIndices_.end();
    }

    int getLocalIdx_(const int scvfIdx) const
    {
        const auto it = std::find(volVarIndices_.begin(), volVarIndices_.end(), scvfIdx);
        assert(it != volVarIndices_.end() && "Could not find the current face variables for scvfIdx!");
        return std::distance(volVarIndices_.begin(), it);
    }

    const GridVolumeVariables* gridVolumeVariablesPtr_;
    std::vector<std::size_t> volVarIndices_;
    std::vector<VolumeVariables> volumeVariables_;
};

} // end namespace

#endif
