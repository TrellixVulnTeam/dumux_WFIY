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
 * \ingroup Discretization
 * \brief Base class for grid geometries
 */
#ifndef DUMUX_DISCRETIZATION_BASE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_BASE_GRID_GEOMETRY_HH

#include <utility>
#include <type_traits>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/deprecated.hh>
#include <dumux/common/entitymap.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

// make the local view function available whenever we use the grid geometry
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Base class for all finite volume grid geometries
 * \tparam GV the grid view type
 * \tparam Traits traits class
 */
template<class GV, class Traits>
class BaseGridGeometry
{
    using ElementMap = EntityMap<GV, 0>;
    using ElementSet = GridViewGeometricEntitySet<GV, 0, typename Traits::ElementMapper>;
    using BoundingBoxTree = Dumux::BoundingBoxTree<ElementSet>;

    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;

    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using Element = typename GV::template Codim<0>::Entity;

public:
    //! export the grid type
    using Grid = typename GV::Grid;
    //! export the grid view type
    using GridView = GV;
    //! export the global coordinate type
    using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
    //! export the element mapper type
    using ElementMapper = typename Traits::ElementMapper;
    //! export the vertex mapper type
    using VertexMapper = typename Traits::VertexMapper;

    /*!
     * \ingroup Discretization
     * \brief Constructor computes the bounding box of the entire domain, for e.g. setting boundary conditions
     * \param gridView the grid view on which to construct the grid geometry
     */
    BaseGridGeometry(const GridView& gridView)
    : gridView_(gridView)
    , elementMapper_(makeElementMapper_(gridView))
    , vertexMapper_(makeVertexMapper_(gridView))
    , bBoxMin_(std::numeric_limits<double>::max())
    , bBoxMax_(-std::numeric_limits<double>::max())
    {
        computeGlobalBoundingBox_();
        update_();
    }

    /*!
     * \brief Update all fvElementGeometries (do this again after grid adaption)
     */
    [[deprecated("Use update(gridView) instead! Will be removed after release 3.5.")]]
    void update()
    {
        update_();
    }

    /*!
     * \brief Update all fvElementGeometries (call this after grid adaption)
     */
    void update(const GridView& gridView)
    {
        gridView_ = gridView;
        update_();
    }

    /*!
     * \brief Update all fvElementGeometries (call this after grid adaption)
     */
    void update(GridView&& gridView)
    {
        gridView_ = std::move(gridView);
        update_();
    }

    /*!
     * \brief Return the gridView this grid geometry object lives on
     */
    const GridView& gridView() const
    { return gridView_; }

    /*!
     * \brief Returns the mapper for vertices to indices for constant grids.
     */
    const VertexMapper &vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for constant grids.
     */
    const ElementMapper &elementMapper() const
    { return elementMapper_; }

    /*!
     * \brief Returns the mapper for vertices to indices for possibly adaptive grids.
     */
    VertexMapper &vertexMapper()
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices for possibly adaptive grids.
     */
    ElementMapper &elementMapper()
    { return elementMapper_; }

    /*!
     * \brief Returns the bounding box tree of the grid
     */
    const BoundingBoxTree& boundingBoxTree() const
    { return *boundingBoxTree_; }

    /*!
     * \brief Returns the element index to element map
     */
    const ElementMap& elementMap() const
    { return *elementMap_; }

    /*!
     * \brief Get an element from a global element index
     */
    Element element(GridIndexType eIdx) const
    { return elementMap()[eIdx]; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalCoordinate &bBoxMin() const
    { return bBoxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalCoordinate &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \brief Returns if the grid geometry is periodic (at all)
     */
    bool isPeriodic() const
    { return periodic_; }

    /*!
     * \brief Set the periodicity of the grid geometry
     */
    void setPeriodic(bool value = true)
    { periodic_ = value; }

private:

    //! Return an instance of the element mapper
    ElementMapper makeElementMapper_(const GridView& gridView) const
    {
        if constexpr (std::is_constructible<ElementMapper, GridView, Dune::MCMGLayout>())
            return ElementMapper(gridView, Dune::mcmgElementLayout());
        else
            return ElementMapper(gridView);
    }

    //! Return an instance of the vertex mapper
    VertexMapper makeVertexMapper_(const GridView& gridView) const
    {
        if constexpr (std::is_constructible<VertexMapper, GridView, Dune::MCMGLayout>())
            return VertexMapper(gridView, Dune::mcmgVertexLayout());
        else
            return VertexMapper(gridView);
    }

    //! Compute the bounding box of the entire domain, for e.g. setting boundary conditions
    void computeGlobalBoundingBox_()
    {
        // calculate the bounding box of the local partition of the grid view
        for (const auto& vertex : vertices(gridView_))
        {
            for (int i=0; i<dimWorld; i++)
            {
                using std::min;
                using std::max;
                bBoxMin_[i] = min(bBoxMin_[i], vertex.geometry().corner(0)[i]);
                bBoxMax_[i] = max(bBoxMax_[i], vertex.geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        if (gridView_.comm().size() > 1)
        {
            for (int i = 0; i < dimWorld; ++i)
            {
                bBoxMin_[i] = gridView_.comm().min(bBoxMin_[i]);
                bBoxMax_[i] = gridView_.comm().max(bBoxMax_[i]);
            }
        }
    }

    void update_()
    {
        // Update the mappers
        if constexpr (Deprecated::hasUpdateGridView<ElementMapper, GridView>())
            elementMapper_.update(gridView_);
        else
            Deprecated::update(elementMapper_);

        if constexpr (Deprecated::hasUpdateGridView<VertexMapper, GridView>())
            vertexMapper_.update(gridView_);
        else
            Deprecated::update(vertexMapper_);

        // Compute the bounding box of the entire domain, for e.g. setting boundary conditions
        computeGlobalBoundingBox_();

        // update element map and bounding box tree
        // always building these comes at a memory overhead but improved
        // performance and thread-safe element level access (e.g. during assembly)
        // for all simulation that use these features
        elementMap_ = std::make_shared<ElementMap>(gridView_.grid(), elementMapper_);
        boundingBoxTree_ = std::make_unique<BoundingBoxTree>(
            std::make_shared<ElementSet>(gridView_, elementMapper(), elementMap_)
        );
    }

    //! the process grid view
    GridView gridView_;

    //! entity mappers
    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    //! the bounding box tree of the grid view for efficient element intersections
    std::unique_ptr<const BoundingBoxTree> boundingBoxTree_;

    //! a map from element index to elements
    std::shared_ptr<const ElementMap> elementMap_;

    //! the bounding box of the whole domain
    GlobalCoordinate bBoxMin_;
    GlobalCoordinate bBoxMax_;

    //! if the grid geometry has periodic boundaries
    bool periodic_ = false;
};

} // end namespace Dumux

#endif
