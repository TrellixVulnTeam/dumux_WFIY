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
 * \brief The flux variables cache filler class for the cell-centered TPFA scheme
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_FLUXVARSCACHE_FILLER_HH
#define DUMUX_DISCRETIZATION_CCTPFA_FLUXVARSCACHE_FILLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Helper class to fill the flux var caches
 */
template<class TypeTag>
class MimeticFluxVariablesCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr bool doAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
    static constexpr bool doDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
    static constexpr bool doHeatConduction = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);

    static constexpr bool soldependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);
    static constexpr bool soldependentDiffusion = GET_PROP_VALUE(TypeTag, SolutionDependentMolecularDiffusion);
    static constexpr bool soldependentHeatConduction = GET_PROP_VALUE(TypeTag, SolutionDependentHeatConduction);

    enum ProcessIndices : unsigned int
    {
        advectionIdx,
        diffusionIdx,
        heatConductionIdx
    };

public:
    //! The constructor. Sets the problem pointer
    MimeticFluxVariablesCacheFiller(const Problem& problem) : problemPtr_(&problem) {}

    /*!
     * \brief function to fill the flux variables caches
     *
     * \param fluxVarsCacheContainer Either the element or global flux variables cache
     * \param scvfFluxVarsCache The flux var cache to be updated corresponding to the given scvf
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param elemVolVars The element volume variables
     * \param scvf The corresponding sub-control volume face
     * \param doSubCaches Array of bools indicating which sub caches have to be updated
     */
    template<class FluxVariablesCacheContainer>
    void fill(FluxVariablesCacheContainer& fluxVarsCacheContainer,
              FluxVariablesCache& scvfFluxVarsCache,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace& scvf,
              const std::array<bool, 3>& doSubCaches = std::array<bool, 3>({true, true, true}))
    {
        // fill the physics-related quantities of the caches
        if (doSubCaches[ProcessIndices::advectionIdx])
          fillAdvection(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        if (doSubCaches[ProcessIndices::diffusionIdx])
          fillDiffusion(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
        if (doSubCaches[ProcessIndices::heatConductionIdx])
          fillHeatConduction(scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf);
    }

    /*!
     * \brief function to update the flux variables caches during derivative calculation
     *
     * \copydoc fill
     */
    template<class FluxVariablesCacheContainer>
    void update(FluxVariablesCacheContainer& fluxVarsCacheContainer,
                FluxVariablesCache& scvfFluxVarsCache,
                const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolumeFace& scvf)
    {
        // array of bool with which we indicate the sub-caches which have to be
        // filled. During update, we only update solution-dependent quantities.
        static const std::array<bool, 3> updateSubCaches = []()
        {
            std::array<bool, 3> tmp;
            tmp[ProcessIndices::advectionIdx] = doAdvection && soldependentAdvection;
            tmp[ProcessIndices::diffusionIdx] = doDiffusion && soldependentDiffusion;
            tmp[ProcessIndices::heatConductionIdx] = doHeatConduction && soldependentHeatConduction;
            return tmp;
        } ();

        // forward to fill routine
        fill(fluxVarsCacheContainer, scvfFluxVarsCache, element, fvGeometry, elemVolVars, scvf, updateSubCaches);
    }

    static bool isSolutionIndependent()
    {
        static const bool isSolDependent = (doAdvection && soldependentAdvection) ||
                                           (doDiffusion && soldependentDiffusion) ||
                                           (doHeatConduction && soldependentHeatConduction);
        return !isSolDependent;
    }

private:

    const Problem& problem() const
    { return *problemPtr_; }

    //! method to fill the advective quantities
    template<bool advectionEnabled = doAdvection>
    typename std::enable_if<advectionEnabled>::type
    fillAdvection(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
        using AdvectionFiller = typename AdvectionType::CacheFiller;

        // forward to the filler for the advective quantities
        AdvectionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! do nothing if advection is not enabled
    template<bool advectionEnabled = doAdvection>
    typename std::enable_if<!advectionEnabled>::type
    fillAdvection(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {}

    //! method to fill the diffusive quantities
    template<bool diffusionEnabled = doDiffusion>
    typename std::enable_if<diffusionEnabled>::type
    fillDiffusion(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
        using DiffusionFiller = typename DiffusionType::CacheFiller;

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

        // forward to the filler of the diffusive quantities
        for (unsigned int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (unsigned int compIdx = 0; compIdx < numComponents; ++compIdx)
                if (phaseIdx != compIdx)
                    DiffusionFiller::fill(scvfFluxVarsCache, phaseIdx, compIdx, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! do nothing if diffusion is not enabled
    template<bool diffusionEnabled = doDiffusion>
    typename std::enable_if<!diffusionEnabled>::type
    fillDiffusion(FluxVariablesCache& scvfFluxVarsCache,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {}

    //! method to fill the quantities related to heat conduction
    template<bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCache& scvfFluxVarsCache,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {
        using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
        using HeatConductionFiller = typename HeatConductionType::CacheFiller;

        // forward to the filler of the diffusive quantities
        HeatConductionFiller::fill(scvfFluxVarsCache, problem(), element, fvGeometry, elemVolVars, scvf, *this);
    }

    //! do nothing if heat conduction is disabled
    template<bool heatConductionEnabled = doHeatConduction>
    typename std::enable_if<!heatConductionEnabled>::type
    fillHeatConduction(FluxVariablesCache& scvfFluxVarsCache,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf)
    {}

    const Problem* problemPtr_;
};

} // end namespace

#endif
