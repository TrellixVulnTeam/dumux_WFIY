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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::StaggeredVelocityGradients
 */
#ifndef DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
#define DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH

#include <optional>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

// forward declare
namespace Dune {

template<class ct, int dim, template< int > class Ref, class Comm>
class SPGrid;
}

namespace Dumux {

namespace Detail
{

template<class Grid>
struct SupportsPeriodicity : public std::false_type {};

template<class ct, int dim, template< int > class Ref, class Comm>
struct SupportsPeriodicity<Dune::SPGrid<ct, dim, Ref, Comm>> : public std::true_type {};
}

/*!
 * \ingroup NavierStokesModel
 * \brief Helper class for calculating the velocity gradients for the Navier-Stokes model using the staggered grid discretization.
 */
class StaggeredVelocityGradients
{

public:

    /*!
     * \brief Returns the in-axis velocity gradient.
     *
     * \verbatim
     *              ---------=======                 == and # staggered half-control-volume
     *              |       #      | current scv
     *              |       #      |                 # staggered face over which fluxes are calculated
     *   vel.Opp <~~|       O~~>   x~~~~> vel.Self
     *              |       #      |                 x dof position
     *              |       #      |
     *              --------========                 -- element
     *
     *                                               O position at which gradient is evaluated (integration point)
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars, class Problem>
    static auto velocityGradII(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars,
                               const Problem& problem)
    {
        assert(scvf.isFrontal());
        // The velocities of the dof at interest and the one of the opposite scvf.
        //const auto& gridGeometry = fvGeometry.gridGeometry();
        //const auto elemSelf = fvGeometry.gridGeometry().element(scvf.insideScvIdx());
        const auto scvSelf = fvGeometry.scv(scvf.insideScvIdx());
        const auto phasefieldSelf = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(scvf.insideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(scvf.insideScvIdx())) + 1e-8
            ;
        const auto phasefieldOpposite = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(scvf.outsideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(scvf.outsideScvIdx())) + 1e-8
            ;
        const auto velocitySelf = elemVolVars[scvf.insideScvIdx()].velocity() *phasefieldSelf;
        const auto velocityOpposite = elemVolVars[scvf.outsideScvIdx()].velocity()
            *phasefieldOpposite;
        const auto distance = (fvGeometry.scv(scvf.outsideScvIdx()).dofPosition() - fvGeometry.scv(scvf.insideScvIdx()).dofPosition()).two_norm();

        return (velocityOpposite - velocitySelf) / distance * scvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient perpendicular to the orientation of our current scvf.
     *
     * \verbatim
     *              ----------------
     *              |              |outer                    || and # staggered half-control-volume (own element)
     *              |              |vel.        gradient
     *              |              |~~~~>       ------->     :: staggered half-control-volume (neighbor element)
     *              |              |             ------>
     *              | lateral scvf |              ----->     x dof position
     *              ---------######O:::::::::      ---->
     *              |      ||      |       ::       --->     -- elements
     *              |      ||      |       ::        -->
     *              |      || scv  x~~~~>  ::         ->     O position at which gradient is evaluated (integration point)
     *              |      ||      | inner ::
     *              |      ||      | vel.  ::
     *              ---------#######:::::::::
     *
     *
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars, class Problem>
    static auto velocityGradIJ(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars,
                               const Problem& problem)
    {
        assert(scvf.isLateral());
        //const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto phasefieldInner = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(scvf.insideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(scvf.insideScvIdx())) + 1e-8
            ;
        const auto phasefieldOuter = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(scvf.outsideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(scvf.outsideScvIdx())) + 1e-8
            ;
        const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity() *phasefieldInner;
        const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity() *phasefieldOuter;
        const auto distance = getDistanceIJ_(fvGeometry, scvf);
        return (outerVelocity - innerVelocity) / distance * scvf.directionSign();
    }

    /*!
     * \brief Returns the velocity gradient in line with our current scvf.
     *
     * \verbatim
     *                      ^       gradient
     *                      |  ^
     *                      |  |  ^
     *                      |  |  |  ^
     *                      |  |  |  |  ^
     *                      |  |  |  |  |  ^       || and # staggered half-control-volume (own element)
     *                      |  |  |  |  |  |
     *                                             :: staggered half-control-volume (neighbor element)
     *              ----------------
     *              |     inner    |      outer    x dof position (of own scv)
     *              |      vel.    |       vel.
     *              |       ^      |        ^      -- elements
     *              |       | lat. |        |
     *              |       | scvf |        |      O position at which gradient is evaluated (integration point)
     *              ---------######O:::::::::
     *              |      ||      ~       ::      ~ ~ orthogonal scvf
     *              |      ||      ~       ::
     *              |      || scv  x       ::
     *              |      ||      |       ::
     *              |      ||      |       ::
     *              ---------#######:::::::::
     *
     *
     * \endverbatim
     */
    template<class FVElementGeometry, class ElemVolVars, class Problem>
    static auto velocityGradJI(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const ElemVolVars& elemVolVars,
                               const Problem& problem)
    {
        assert(scvf.isLateral());
        const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
        //const auto& gridGeometry = fvGeometry.gridGeometry();
        const auto phasefieldInner = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(orthogonalScvf.insideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(orthogonalScvf.insideScvIdx())) + 1e-8
            ;
        const auto phasefieldOuter = 1.0
            //* problem.phasefield(fvGeometry.gridGeometry().element(orthogonalScvf.outsideScvIdx()),
            ///*fvGeometry,*/ fvGeometry.scv(orthogonalScvf.outsideScvIdx())) + 1e-8
            ;
        const auto innerVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity()
            *phasefieldInner;
        const auto outerVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity()
            *phasefieldOuter;
        const auto distance = getDistanceJI_(fvGeometry, scvf, orthogonalScvf);
        return (outerVelocity - innerVelocity) / distance * orthogonalScvf.directionSign();
    }

private:

    template<class FVElementGeometry>
    static auto getDistanceIJ_(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        // The scvf is on a boundary, hence there is no outer DOF.
        // We take the distance to the boundary instead.
        if (scvf.boundary())
            return getDistanceToBoundary_(insideScv, scvf);

        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        // The standard case: Our grid does not support periodic boundaries.
        if constexpr (!Detail::SupportsPeriodicity<typename FVElementGeometry::GridGeometry::Grid>())
            return getStandardDistance_(insideScv, outsideScv);
        else
        {
            const auto& gridGeometry = fvGeometry.gridGeometry();
            const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
            const auto& orthogonalInsideScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());

            // The standard case: Our grid is not periodic or the lateral scvf does not lie on a periodic boundary.
            if (!gridGeometry.isPeriodic() || !gridGeometry.dofOnPeriodicBoundary(orthogonalInsideScv.dofIndex()))
                return getStandardDistance_(insideScv, outsideScv);

            // Treating periodic boundaries is more involved:
            // 1. Consider the outside scv within the element adjacent to the other periodic boundary.
            // 2. Iterate over this scv's faces until you find the face parallel to our own scvf.
            //    This face would lie directly next to our own scvf if the grid was not periodic.
            // 3. Calculate the total distance by summing up the distances between the DOFs and the respective integration points
            //    corresponding to the inside and outside scvfs.
            auto periodicFvGeometry = localView(gridGeometry);
            const auto& periodicElement = gridGeometry.element(outsideScv.elementIndex());
            periodicFvGeometry.bindElement(periodicElement);

            for (const auto& outsideScvf : scvfs(periodicFvGeometry, outsideScv))
            {
                if (outsideScvf.directionIndex() == scvf.directionIndex() && outsideScvf.directionSign() != scvf.directionSign())
                {
                    const auto insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                    const auto outsideDistance = (outsideScv.dofPosition() - outsideScvf.ipGlobal()).two_norm();
                    return insideDistance + outsideDistance;
                }
            }
            DUNE_THROW(Dune::InvalidStateException, "Could not find scvf in periodic element");
        }
    }

    //! Get the distance between the DOFs at which the inner and outer velocities are defined.
    template<class FVElementGeometry>
    static auto getDistanceJI_(const FVElementGeometry fvGeometry,
                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                               const typename FVElementGeometry::SubControlVolumeFace& orthogonalScvf)
    {
        const auto& orthogonalInsideScv = fvGeometry.scv(orthogonalScvf.insideScvIdx());

        // The orthogonal scvf is on a boundary, hence there is no outer DOF.
        // We take the distance to the boundary instead.
        if (orthogonalScvf.boundary())
            return getDistanceToBoundary_(orthogonalInsideScv, orthogonalScvf);

        const auto& orthogonalOutsideScv = fvGeometry.scv(orthogonalScvf.outsideScvIdx());

        // The standard case: Our grid does not support periodic boundaries.
        if constexpr (!Detail::SupportsPeriodicity<typename FVElementGeometry::GridGeometry::Grid>())
            return getStandardDistance_(orthogonalInsideScv, orthogonalOutsideScv);
        else
        {
            const auto& gridGeometry = fvGeometry.gridGeometry();
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            // The standard case: Our grid is not periodic or our own DOF does not lie on a periodic boundary.
            if (!gridGeometry.isPeriodic() || !gridGeometry.dofOnPeriodicBoundary(insideScv.dofIndex()))
                return getStandardDistance_(orthogonalInsideScv, orthogonalOutsideScv);

            // Treating periodic boundaries is more involved:
            // 1. Consider the orthogonal outside scv within the element adjacent to the other periodic boundary.
            // 2. Iterate over this scv's faces until you find the orthogonal face parallel to our own orthogonal scvf.
            //    This face would lie directly next to our own orthogonal scvf if the grid was not periodic.
            // 3. Calculate the total distance by summing up the distances between the DOFs and the respective integration points
            //    corresponding to the orthogonal scvfs.
            auto periodicFvGeometry = localView(gridGeometry);
            const auto& periodicElement = gridGeometry.element(orthogonalOutsideScv.elementIndex());
            periodicFvGeometry.bindElement(periodicElement);

            for (const auto& outsideOrthogonalScvf : scvfs(periodicFvGeometry, orthogonalOutsideScv))
            {
                if (outsideOrthogonalScvf.directionIndex() == orthogonalScvf.directionIndex() && outsideOrthogonalScvf.directionSign() != orthogonalScvf.directionSign())
                {
                    const auto insideDistance = (orthogonalInsideScv.dofPosition() - orthogonalScvf.ipGlobal()).two_norm();
                    const auto outsideDistance = (orthogonalOutsideScv.dofPosition() - outsideOrthogonalScvf.ipGlobal()).two_norm();
                    return insideDistance + outsideDistance;
                }
            }
            DUNE_THROW(Dune::InvalidStateException, "Could not find scvf in periodic element");
        }
    }

    template<class SubControlVolume>
    static auto getStandardDistance_(const SubControlVolume& insideScv,
                                     const SubControlVolume& outsideScv)
    {
        return (insideScv.dofPosition() - outsideScv.dofPosition()).two_norm();
    }

    template<class SubControlVolume, class SubControlVolumeFace>
    static auto getDistanceToBoundary_(const SubControlVolume& scv,
                                       const SubControlVolumeFace& scvf)
    {
        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
    }
};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_VELOCITYGRADIENTS_HH
