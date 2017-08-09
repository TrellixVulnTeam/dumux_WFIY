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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_DARCYS_LAW_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
NEW_PROP_TAG(MpfaHelper);
NEW_PROP_TAG(BoundaryInteractionVolume);
NEW_PROP_TAG(UseTpfaBoundary);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCMpfa method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCMpfa>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, AdvectionType);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    // Always use the dynamic type for vectors (compatibility with the boundary)
    using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using CoefficientVector = typename BoundaryInteractionVolume::Vector;

    static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    //! The cache used in conjunction with the mpfa Darcy's Law
    class MpfaDarcysLawCache
    {
        // We always use the dynamic types here to be compatible on the boundary
        using Stencil = typename BoundaryInteractionVolume::GlobalIndexSet;
        using PositionVector = typename BoundaryInteractionVolume::PositionVector;

    public:
        //! update cached objects
        template<class InteractionVolume>
        void updateAdvection(const InteractionVolume& iv, const SubControlVolumeFace &scvf)
        {
            const auto& localFaceData = iv.getLocalFaceData(scvf);
            // update the quantities that are equal for all phases
            advectionVolVarsStencil_ = iv.volVarsStencil();
            advectionVolVarsPositions_ = iv.volVarsPositions();
            advectionTij_ = iv.getTransmissibilities(localFaceData);

            // The neumann fluxes always have to be set per phase
            for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                phaseNeumannFluxes_[phaseIdx] = iv.getNeumannFlux(localFaceData, phaseIdx);
        }

        //! Returns the volume variables indices necessary for flux computation
        //! This includes all participating boundary volume variables. Since we
        //! do not allow mixed BC for the mpfa this is the same for all phases.
        const Stencil& advectionVolVarsStencil() const
        { return advectionVolVarsStencil_; }

        //! Returns the position on which the volume variables live. This is
        //! necessary as we need to evaluate gravity also for the boundary volvars
        const PositionVector& advectionVolVarsPositions() const
        { return advectionVolVarsPositions_; }

        //! Returns the transmissibilities associated with the volume variables
        //! All phases flow through the same rock, thus, tij are equal for all phases
        const CoefficientVector& advectionTij() const
        { return advectionTij_; }

        //! If the useTpfaBoundary property is set to false, the boundary conditions
        //! are put into the local systems leading to possible contributions on all faces
        Scalar advectionNeumannFlux(unsigned int phaseIdx) const
        { return phaseNeumannFluxes_[phaseIdx]; }

    private:
        // Quantities associated with advection
        Stencil advectionVolVarsStencil_;
        PositionVector advectionVolVarsPositions_;
        CoefficientVector advectionTij_;
        std::array<Scalar, numPhases> phaseNeumannFluxes_;
    };

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MpfaDarcysLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any advection-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            // get interaction volume from the flux vars cache filler & upate the cache
            if (problem.model().globalFvGeometry().isInBoundaryInteractionVolume(scvf))
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.boundaryInteractionVolume(), scvf);
            else
                scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.interactionVolume(), scvf);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCMpfa;

    // state the type for the corresponding cache and its filler
    using Cache = MpfaDarcysLawCache;
    using CacheFiller = MpfaDarcysLawCacheFiller;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const unsigned int phaseIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        static const bool gravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& volVarsStencil = fluxVarsCache.advectionVolVarsStencil();
        const auto& volVarsPositions = fluxVarsCache.advectionVolVarsPositions();
        const auto& tij = fluxVarsCache.advectionTij();

        const bool isInteriorBoundary = enableInteriorBoundaries && fluxVarsCache.isInteriorBoundary();
        // For interior Neumann boundaries when using Tpfa on boundaries, return the user-specified flux
        // We assume phaseIdx = eqIdx here.
        if (isInteriorBoundary
            && useTpfaBoundary
            && fluxVarsCache.interiorBoundaryDataSelf().faceType() == MpfaFaceTypes::interiorNeumann)
            return scvf.area()*
                   elemVolVars[scvf.insideScvIdx()].extrusionFactor()*
                   problem.neumann(element, fvGeometry, elemVolVars, scvf)[phaseIdx];

        // Calculate the interface density for gravity evaluation
        const auto rho = Implementation::interpolateDensity(fvGeometry, elemVolVars, scvf, fluxVarsCache, phaseIdx, isInteriorBoundary);

        // calculate Tij*pj
        Scalar flux(0.0);
        unsigned int localIdx = 0;
        for (const auto volVarIdx : volVarsStencil)
        {
            const auto& volVars = elemVolVars[volVarIdx];
            Scalar h = volVars.pressure(phaseIdx);

            // if gravity is enabled, add gravitational acceleration
            if (gravity)
            {
                // gravitational acceleration in the center of the actual element
                const auto x = volVarsPositions[localIdx];
                const auto g = problem.gravityAtPos(x);

                h -= rho*(g*x);
            }

            flux += tij[localIdx++]*h;
        }

        // if no interior boundaries are present, return the flux
        if (!enableInteriorBoundaries)
            return useTpfaBoundary ? flux : flux + fluxVarsCache.advectionNeumannFlux(phaseIdx);

        // Handle interior boundaries
        flux += Implementation::computeInteriorBoundaryContribution(problem, fvGeometry, elemVolVars, fluxVarsCache, phaseIdx, rho);

        // return overall resulting flux
        return useTpfaBoundary ? flux : flux + fluxVarsCache.advectionNeumannFlux(phaseIdx);
    }

    static Scalar interpolateDensity(const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const SubControlVolumeFace& scvf,
                                     const FluxVariablesCache& fluxVarsCache,
                                     const unsigned int phaseIdx,
                                     const bool isInteriorBoundary)
    {
        static const bool gravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        if (!gravity)
            return Scalar(0.0);
        else
        {
            // Treat interior Dirichlet boundaries differently
            if (isInteriorBoundary)
            {
                const auto& data = fluxVarsCache.interiorBoundaryDataSelf();
                if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
                    return data.facetVolVars(fvGeometry).density(phaseIdx);
            }

            // use arithmetic mean of the densities around the scvf
            if (!scvf.boundary())
            {
                Scalar rho = elemVolVars[scvf.insideScvIdx()].density(phaseIdx);
                for (auto outsideIdx : scvf.outsideScvIndices())
                    rho += elemVolVars[outsideIdx].density(phaseIdx);
                return rho/(scvf.outsideScvIndices().size()+1);
            }
            else
                return elemVolVars[scvf.outsideScvIdx()].density(phaseIdx);
        }
    }

    static Scalar computeInteriorBoundaryContribution(const Problem& problem,
                                                      const FVElementGeometry& fvGeometry,
                                                      const ElementVolumeVariables& elemVolVars,
                                                      const FluxVariablesCache& fluxVarsCache,
                                                      unsigned int phaseIdx, Scalar rho)
    {
        static const bool gravity = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);

        // obtain the transmissibilites associated with all pressures
        const auto& tij = fluxVarsCache.advectionTij();

        // the interior dirichlet boundaries local indices start after
        // the cell and the domain Dirichlet boundary pressures
        const auto startIdx = fluxVarsCache.advectionVolVarsStencil().size();

        // add interior Dirichlet boundary contributions
        Scalar flux = 0.0;
        for (auto&& data : fluxVarsCache.interiorBoundaryData())
        {
            if (data.faceType() == MpfaFaceTypes::interiorDirichlet)
            {
                Scalar h = data.facetVolVars(fvGeometry).pressure(phaseIdx);

                if (gravity)
                {
                    const auto x = fvGeometry.scvf(data.scvfIndex()).ipGlobal();
                    const auto g = problem.gravityAtPos(x);

                    h -= rho*(g*x);
                }

                flux += tij[startIdx + data.localIndexInInteractionVolume()]*h;
            }
        }

        return flux;
    }
};

} // end namespace

#endif
