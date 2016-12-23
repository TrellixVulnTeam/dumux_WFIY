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
 * \brief This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the box method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElemFluxVarCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;
    using Stencil = std::vector<IndexType>;

    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const IndexType phaseIdx,
                       const ElemFluxVarCache& elemFluxVarCache)
    {
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        auto insideK = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();

        // scale with correct extrusion factor
        insideK *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();

        const auto K = problem.spatialParams().meanDiffusionTensor(insideK, outsideK, scvf.unitOuterNormal());

        const auto& jacInvT = fluxVarCache.jacInvT();
        const auto& shapeJacobian = fluxVarCache.shapeJacobian();
        const auto& shapeValues = fluxVarCache.shapeValues();

        static const bool enableGravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        // evaluate gradP - rho*g at integration point
        DimVector gradP(0.0);
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            if (enableGravity)
                rho += volVars.density(phaseIdx)*shapeValues[scv.index()][0];

            // the global shape function gradient
            DimVector gradI;
            jacInvT.mv(shapeJacobian[scv.index()][0], gradI);

            gradI *= volVars.pressure(phaseIdx);
            gradP += gradI;
        }
        if (enableGravity)
        {
            // gravitational acceleration
            DimVector g(problem.gravityAtPos(scvf.center()));

            // turn gravity into a force
            g *= rho;

            // subtract from pressure gradient
            gradP -= g;
        }

        // apply the permeability and return the flux
        auto KGradP = applyPermeability_(K, gradP);
        return -1.0*(KGradP*scvf.unitOuterNormal())*scvf.area();
    }

    // This is for compatibility with the cc methods. The flux stencil info is obsolete for the box method.
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        return Stencil(0);
    }
private:
    static DimVector applyPermeability_(const DimWorldMatrix& K, const DimVector& gradI)
    {
        DimVector result(0.0);
        K.mv(gradI, result);

        return result;
    }

    static DimVector applyPermeability_(const Scalar k, const DimVector& gradI)
    {
        DimVector result(gradI);
        result *= k;
        return result;
    }
};

} // end namespace Dumux

#endif
