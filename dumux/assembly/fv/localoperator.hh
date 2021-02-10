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
 * \ingroup Assembly
 * \brief An element-wise local operator for finite-volume schemes.
 */
#ifndef DUMUX_FV_LOCAL_OPERATOR_HH
#define DUMUX_FV_LOCAL_OPERATOR_HH

#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/istl/bvector.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise local operator for finite volume schemes.
 *        This allows for element-wise evaluation of individual terms
 *        of the equations to be solved.
 * \tparam ElementVariables local view on the grid variables
 * \tparam Op The model-specific operators
 */
template<class ElementVariables, class Op>
class FVLocalOperator
{
    // The variables required for the evaluation of the equation
    using GridVars = typename ElementVariables::GridVariables;
    using PrimaryVariables = typename GridVars::PrimaryVariables;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVars::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using NumEqVector = typename Op::NumEqVector;
    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethod::box;

public:
    //! export operators
    using Operators = Op;

    //! export the grid variables type this operator requires a local view of
    using GridVariables = GridVars;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    //! TODO: Necessary?
    using Variables = GridVars;

    //! the container storing the operator values on all dofs of an element
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    //! TODO: Necessary?
    using Residual = ElementResidualVector;

    /*!
     * \brief The constructor
     * \note The grid geometry/grid variables local views are expected to
     *       be bound to the same (the given) element
     */
    explicit FVLocalOperator(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elemVariables_(elemVars)
    {}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the terms of the local residual that do not appear in
     *        time derivatives. These are the sources and the fluxes.
     */
    ElementResidualVector evalFluxesAndSources() const
    {
        const auto& problem = elemVariables_.gridVariables().gridVolVars().problem();
        const auto& evv = elemVariables_.elemVolVars();

        // source term
        auto result = getEmptyResidual();
        for (const auto& scv : scvs(fvGeometry_))
            result[scv.localDofIndex()] -= Operators::source(problem, element_, fvGeometry_, evv, scv);

        // flux term
        for (const auto& scvf : scvfs(fvGeometry_))
            addFlux_(result, scvf);

        return result;
    }

    /*!
     * \brief Compute the storage term, i.e. the term appearing in the time derivative.
     */
    ElementResidualVector evalStorage() const
    {
        const auto& problem = elemVariables_.gridVariables().gridVolVars().problem();
        const auto& evv = elemVariables_.elemVolVars();

        // TODO: Until now, FVLocalResidual defined storage as the entire
        //       time derivative. Now it means the term above the time derivative.
        //       We should think about the correct naming here...
        // TODO: Should storage() NOT multiply with volume?? That was different until
        //       now but flux() also returns the flux multiplied with area so this should
        //       be more consistent
        auto result = getEmptyResidual();
        for (const auto& scv : scvs(fvGeometry_))
            result[scv.localDofIndex()] += operators_.storage(problem, scv, evv[scv]);

        return result;
    }

    ElementResidualVector getEmptyResidual() const
    {
        ElementResidualVector res(fvGeometry_.numScv());
        res = 0.0;
        return res;
    }

    // \}

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces. Or, should this be here at all!?

    //\}

    // \}

protected:

    //! compute and add the flux across the given face to the container (cc schemes)
    template<bool b = isBox, std::enable_if_t<!b, int> = 0>
    void addFlux_(ElementResidualVector& r, const SubControlVolumeFace& scvf) const
    {
        const auto& insideScv = fvGeometry_.scv(scvf.insideScvIdx());
        const auto localDofIdx = insideScv.localDofIndex();

        const auto& problem = elemVariables_.gridVariables().gridVolVars().problem();
        const auto& evv = elemVariables_.elemVolVars();
        const auto& efvc = elemVariables_.elemFluxVarsCache();

        if (!scvf.boundary())
            r[localDofIdx] += Operators::flux(problem, element_, fvGeometry_, evv, efvc, scvf);
        else
        {
            const auto& bcTypes = problem.boundaryTypes(element_, scvf);

            // Dirichlet boundaries
            if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                r[localDofIdx] += Operators::flux(problem, element_, fvGeometry_, evv, efvc, scvf);

            // Neumann and Robin ("solution dependent Neumann") boundary conditions
            else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
            {
                auto neumannFluxes = problem.neumann(element_, fvGeometry_, evv, efvc, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                neumannFluxes *= Extrusion::area(scvf)*evv[insideScv].extrusionFactor();
                r[localDofIdx] += neumannFluxes;
            }

            else
                DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions for cell-centered schemes. " <<
                                                 "Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
        }
    }

    //! compute and add the flux across the given face to the container (box scheme)
    template<bool b = isBox, std::enable_if_t<b, int> = 0>
    void addFlux_(ElementResidualVector& r, const SubControlVolumeFace& scvf) const
    {
        const auto& problem = elemVariables_.gridVariables().gridVolVars().problem();
        const auto& evv = elemVariables_.elemVolVars();
        const auto& efvc = elemVariables_.elemFluxVarsCache();

        // inner faces
        if (!scvf.boundary())
        {
            const auto flux = Operators::flux(problem, element_, fvGeometry_, evv, efvc, scvf);
            r[fvGeometry_.scv(scvf.insideScvIdx()).localDofIndex()] += flux;
            r[fvGeometry_.scv(scvf.outsideScvIdx()).localDofIndex()] -= flux;
        }

        // boundary faces
        else
        {
            const auto& scv = fvGeometry_.scv(scvf.insideScvIdx());
            const auto& bcTypes = problem.boundaryTypes(element_, scv);

            // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
            if (bcTypes.hasNeumann())
            {
                const auto neumannFluxes = problem.neumann(element_, fvGeometry_, evv, efvc, scvf);
                const auto area = Extrusion::area(scvf)*evv[scv].extrusionFactor();

                // only add fluxes to equations for which Neumann is set
                for (int eqIdx = 0; eqIdx < NumEqVector::size(); ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        r[scv.localDofIndex()][eqIdx] += neumannFluxes[eqIdx]*area;
            }
        }
    }

private:

    const Element& element_;                 //!< pointer to the element for which the residual is computed
    const FVElementGeometry& fvGeometry_;    //!< the local view on the finite element grid geometry
    const ElementVariables& elemVariables_;  //!< the local view on the grid variables
    Operators operators_; //!< evaluates storage/flux operators of the actual equation
};

} // end namespace Dumux

#endif
