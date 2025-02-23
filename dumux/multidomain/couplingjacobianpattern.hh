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
 * \ingroup MultiDomain
 * \brief Helper function to generate Jacobian pattern for multi domain models
 */
#ifndef DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH
#define DUMUX_MUTLIDOMAIN_COUPLING_JACOBIAN_PATTERN_HH

#include <type_traits>
#include <dune/common/indices.hh>
#include <dune/istl/matrixindexset.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for cell-centered schemes
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<( (GridGeometryI::discMethod == DiscretizationMethods::cctpfa)
                                     || (GridGeometryI::discMethod == DiscretizationMethods::ccmpfa) ), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    const auto numDofsI = gridGeometryI.numDofs();
    const auto numDofsJ = gridGeometryJ.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofsI, numDofsJ);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            const auto globalI = gridGeometryI.elementMapper().index(elementI);
            for (const auto globalJ : stencil)
                pattern.add(globalI, globalJ);
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the box scheme
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::box), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    const auto numDofsI = gridGeometryI.numDofs();
    const auto numDofsJ = gridGeometryJ.numDofs();
    Dune::MatrixIndexSet pattern;
    pattern.resize(numDofsI, numDofsJ);

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        static constexpr int dim = std::decay_t<decltype(gridGeometryI.gridView())>::dimension;
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            for (std::size_t vIdxLocal = 0; vIdxLocal < elementI.subEntities(dim); ++vIdxLocal)
            {
                const auto globalI = gridGeometryI.vertexMapper().subIndex(elementI, vIdxLocal, dim);
                for (const auto globalJ : stencil)
                    pattern.add(globalI, globalJ);
            }
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on cell centers)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::staggered &&
                                    GridGeometryI::isCellCenter()), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        const auto ccGlobalI = gridGeometryI.elementMapper().index(elementI);
        for (auto&& faceGlobalJ : couplingManager.couplingStencil(domainI, elementI, domainJ))
                pattern.add(ccGlobalI, faceGlobalJ);
    }

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on faces)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::staggered &&
                                    GridGeometryI::isFace()), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    auto fvGeometry = localView(gridGeometryI);
    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        fvGeometry.bindElement(elementI);

        // loop over sub control faces
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto faceGlobalI = scvf.dofIndex();
            for (auto&& globalJ : couplingManager.couplingStencil(domainI, scvf, domainJ))
                pattern.add(faceGlobalI, globalJ);
        }
    }

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the staggered scheme (degrees of freedom on cell centers)
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::fcstaggered), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern(gridGeometryI.numDofs(), gridGeometryJ.numDofs());

    auto fvGeometry = localView(gridGeometryI);
    for (const auto& elementI : elements(gridGeometryI.gridView()))
    {
        fvGeometry.bindElement(elementI);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto globalI = scv.dofIndex();
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, scv, domainJ);
            for (const auto globalJ : stencil)
            {
                assert(globalJ < gridGeometryJ.numDofs());
                pattern.add(globalI, globalJ);

                if (gridGeometryI.isPeriodic())
                {
                    if (gridGeometryI.dofOnPeriodicBoundary(globalI))
                    {
                        const auto globalIP = gridGeometryI.periodicallyMappedDof(globalI);

                        if (globalI > globalIP)
                            pattern.add(globalIP, globalJ);
                    }
                }
            }
        }
    }

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the diamond scheme
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::fcdiamond), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern;

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        pattern.resize(gridGeometryI.numDofs(),  gridGeometryJ.numDofs());
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            for (std::size_t localFacetIndex = 0; localFacetIndex < elementI.subEntities(1); ++localFacetIndex)
            {
                const auto globalI = gridGeometryI.dofMapper().subIndex(elementI, localFacetIndex, 1);
                for (const auto globalJ : stencil)
                    pattern.add(globalI, globalJ);
            }
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

/*!
 * \ingroup MultiDomain
 * \brief Helper function to generate coupling Jacobian pattern (off-diagonal blocks)
 *        for the pq1bubble scheme
 */
template<bool isImplicit, class CouplingManager, class GridGeometryI, class GridGeometryJ, std::size_t i, std::size_t j,
         typename std::enable_if_t<(GridGeometryI::discMethod == DiscretizationMethods::pq1bubble), int> = 0>
Dune::MatrixIndexSet getCouplingJacobianPattern(const CouplingManager& couplingManager,
                                                Dune::index_constant<i> domainI,
                                                const GridGeometryI& gridGeometryI,
                                                Dune::index_constant<j> domainJ,
                                                const GridGeometryJ& gridGeometryJ)
{
    Dune::MatrixIndexSet pattern;

    // matrix pattern for implicit Jacobians
    if (isImplicit)
    {
        pattern.resize(gridGeometryI.numDofs(),  gridGeometryJ.numDofs());
        auto fvGeometry = localView(gridGeometryI);
        for (const auto& elementI : elements(gridGeometryI.gridView()))
        {
            fvGeometry.bindElement(elementI);
            const auto& stencil = couplingManager.couplingStencil(domainI, elementI, domainJ);
            for (const auto& scv : scvs(fvGeometry))
            {
                for (const auto globalJ : stencil)
                    pattern.add(scv.dofIndex(), globalJ);

            }
        }
    }

    // matrix pattern for explicit Jacobians
    // -> diagonal matrix, so coupling block is empty
    // just return the empty pattern

    return pattern;
}

} // end namespace Dumux

#endif
