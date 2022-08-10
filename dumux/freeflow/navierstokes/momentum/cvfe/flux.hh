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
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CVFE_FLUXVARIABLES_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#include <dumux/freeflow/navierstokes/momentum/diamond/flux.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

template<class Problem,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class ElementFluxVariablesCache>
using NavierStokesMomentumFluxContextCvfe = NavierStokesMomentumFluxContext<Problem, FVElementGeometry, ElementVolumeVariables, ElementFluxVariablesCache>;

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the cvfe grid discretization
 */
template<class GridGeometry, class NumEqVector>
class NavierStokesMomentumFluxCvfe
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename NumEqVector::value_type;

    using Extrusion = Extrusion_t<GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    static_assert(NumEqVector::dimension == dimWorld, "Wrong dimension of velocity vector");

public:
    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    template<class Context>
    NumEqVector advectiveMomentumFlux(const Context& context) const
    {
        if (!context.problem().enableInertiaTerms())
            return NumEqVector(0.0);

        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();

        // interpolate velocity at scvf
        NumEqVector v(0.0);
        for (const auto& scv : scvs(fvGeometry))
            v.axpy(shapeValues[scv.indexInElement()][0], elemVolVars[scv].velocity());

        // get density from the problem
        const Scalar density = context.problem().density(context.element(), context.fvGeometry(), scvf);

        const auto vn = v*scvf.unitOuterNormal();
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto upwindVelocity = vn > 0 ? insideVolVars.velocity() : outsideVolVars.velocity();
        const auto downwindVelocity = vn > 0 ? outsideVolVars.velocity() : insideVolVars.velocity();
        static const auto upwindWeight = getParamFromGroup<Scalar>(context.problem().paramGroup(), "Flux.UpwindWeight");
        const auto advectiveTermIntegrand = density*vn * (upwindWeight * upwindVelocity + (1.0-upwindWeight)*downwindVelocity);

        return advectiveTermIntegrand * Extrusion::area(scvf) * insideVolVars.extrusionFactor();
    }

    /*!
     * \brief Returns the diffusive momentum flux due to viscous forces
     */
    template<class Context>
    NumEqVector diffusiveMomentumFlux(const Context& context) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];

        // interpolate velocity gradient at scvf
        Tensor gradV(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            for (int dir = 0; dir < dim; ++dir)
                gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(scv.indexInElement()));
        }

        // get viscosity from the problem
        const auto mu = context.problem().effectiveViscosity(element, fvGeometry, scvf);

        // TODO enable symmetric gradient
        static const bool enableUnsymmetrizedVelocityGradient
            = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradient", false);

        // compute -mu*gradV*n*dA
        NumEqVector diffusiveFlux = enableUnsymmetrizedVelocityGradient ?
                mv(gradV, scvf.unitOuterNormal())
                : mv(gradV + getTransposed(gradV),scvf.unitOuterNormal());

        diffusiveFlux *= -mu;

        static const bool enableDilatationTerm = getParamFromGroup<bool>(context.problem().paramGroup(), "FreeFlow.EnableDilatationTerm", false);
        if (enableDilatationTerm)
            diffusiveFlux += 2.0/3.0 * mu * trace(gradV) * scvf.unitOuterNormal();

        diffusiveFlux *= Extrusion::area(scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        return diffusiveFlux;
    }

    template<class Context>
    NumEqVector pressureContribution(const Context& context) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();

        // The pressure force needs to take the extruded scvf area into account
        const auto pressure = context.problem().pressure(element, fvGeometry, scvf);

        // The pressure contribution calculated above might have a much larger numerical value compared to the viscous or inertial forces.
        // This may lead to numerical inaccuracies due to loss of significance (cancellation) for the final residual value.
        // In the end, we are only interested in a pressure gradient between the two relevant faces so we can
        // subtract a constant reference value from the actual pressure contribution.
        const auto referencePressure = context.problem().referencePressure();

        NumEqVector pn(scvf.unitOuterNormal());
        pn *= (pressure-referencePressure)*Extrusion::area(scvf)*elemVolVars[scvf.insideScvIdx()].extrusionFactor();

        return pn;
    }
};

} // end namespace Dumux

#endif
