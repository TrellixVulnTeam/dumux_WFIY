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
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Implementation of the transmissibility laws for throats
 */
#ifndef DUMUX_PNM_THROAT_TRANSMISSIBILITY_2P_HH
#define DUMUX_PNM_THROAT_TRANSMISSIBILITY_2P_HH

#include <dumux/porenetwork/common/throatproperties.hh>
#include "emptycache.hh"

namespace Dumux::PoreNetwork {

namespace WettingLayerTransmissibility {

struct CreviceResistanceFactorZhou
{
    /*!
    * \brief Returns the crevice resistance factor used for calculating the w-phase conductance in an invaded pore throat
    *
    * \param alpha The corner half angle
    * \param theta The contact angle
    * \param f The boundary condition factor for the fluid interface (0 = free boundary)
    *
    * Zhou et al. (1997), eq. 19; Singh & Mohanty (2002), eq 8
    */
    template<class Scalar>
    static Scalar beta(const Scalar alpha, const Scalar theta, const Scalar f = 0)
    {
       using std::sin;
       using std::cos;
       using std::tan;
       const Scalar sinAlpha = sin(alpha);
       const Scalar sinSum = sin(alpha + theta);
       const Scalar cosSum = cos(alpha + theta);
       const Scalar phi1 = cosSum*cosSum + cosSum*sinSum*tan(alpha);
       const Scalar phi2 = 1 - theta/(M_PI/2 - alpha);
       const Scalar phi3 = cosSum / cos(alpha);
       const Scalar B = (M_PI/2 - alpha)*tan(alpha);

       Scalar result = 12*sinAlpha*sinAlpha*(1-B)*(1-B)*(phi1 - B*phi2)*(phi3 + f*B*phi2)*(phi3 + f*B*phi2);
       result /= (1-sinAlpha)*(1-sinAlpha)*B*B*(phi1 - B*phi2)*(phi1 - B*phi2)*(phi1 - B*phi2);
       return result;
    }
};


template<class Scalar, class CreviceResistanceFactor = CreviceResistanceFactorZhou>
struct RansohoffRadke
{
    class WettingLayerCache
    {
        using NumCornerVector = Dune::ReservedVector<Scalar, 4>;
    public:
        template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
        void fill(const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                  const ElementVolumeVariables& elemVolVars,
                  const FluxVariablesCache& fluxVarsCache,
                  const int phaseIdx)
        {
            const auto& spatialParams = problem.spatialParams();
            const auto& cornerHalfAngles = spatialParams.cornerHalfAngles(element);
            const Scalar contactAngle = spatialParams.contactAngle(element, elemVolVars);

            beta_.clear(); beta_.resize(cornerHalfAngles.size());
            for (int i = 0; i< cornerHalfAngles.size(); ++i)
                beta_[i] = CreviceResistanceFactor::beta(cornerHalfAngles[i], contactAngle);
        }

        Scalar creviceResistanceFactor(const int cornerIdx) const
        { return beta_[cornerIdx]; }

    private:
        NumCornerVector beta_;
    };

    /*!
    * \brief Returns the integral conductivity of all wetting layers occupying the corners of a throat
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar wettingLayerTransmissibility(const Element& element,
                                               const FVElementGeometry& fvGeometry,
                                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar rC = fluxVarsCache.curvatureRadius();

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto shape = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        const auto numCorners = Throat::numCorners(shape);

        // treat the wetting film layer in each corner of the throat individually (might have different corner half-angle and beta)
        Scalar result = 0.0;
        for (int i = 0; i < numCorners; ++i)
            result += fluxVarsCache.wettingLayerCrossSectionalArea(i) * rC*rC / (throatLength*fluxVarsCache.wettingLayerFlowVariables().creviceResistanceFactor(i));

        return result;
    }

    /*!
    * \brief Returns the integral conductivity of all wetting layers at critical status entry pressure
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar entryWettingLayerTransmissibility(const Element& element,
                                                    const FVElementGeometry& fvGeometry,
                                                    const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                    const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar rC = fluxVarsCache.curvatureRadiusInvasion(0);
        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto shape = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        const auto numCorners = Throat::numCorners(shape);

        // treat the wetting film layer in each corner of the throat individually (might have different corner half-angle and beta)
        Scalar result = 0.0;
        for (int i = 0; i < numCorners; ++i)
            result += fluxVarsCache.entryWettingLayerArea(i) * rC * rC / (throatLength * fluxVarsCache.wettingLayerFlowVariables().creviceResistanceFactor(i));

        return result;
    }

    /*!
    * \brief Returns the integral conductivity of all wetting layers at critical status snapoff
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar snapoffWettingLayerTransmissibility(const Element& element,
                                                      const FVElementGeometry& fvGeometry,
                                                      const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                      const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar rC = fluxVarsCache.curvatureRadiusSnapoff(0);

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto shape = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        const auto numCorners = Throat::numCorners(shape);

        // treat the wetting film layer in each corner of the throat individually (might have different corner half-angle and beta)
        Scalar result = 0.0;
        for (int i = 0; i < numCorners; ++i)
            result += fluxVarsCache.snapoffWettingLayerArea(i) * rC * rC / (throatLength * fluxVarsCache.wettingLayerFlowVariables().creviceResistanceFactor(i));

        return result;
    }

    /*!
     * \brief Returns derivative of wetting layer transmissibility at entry pressure
     *        For convinience now we calculate the derivative numerically
     *        TODO: Implement the analytical formular of this derivative
     */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar dKwdPcEntry(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                              const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar Kw = entryWettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar rCdelta = fluxVarsCache.curvatureRadiusInvasion(1);

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto shape = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        const auto numCorners = Throat::numCorners(shape);

        // treat the wetting film layer in each corner of the throat individually (might have different corner half-angle and beta)
        Scalar deltaKw = 0.0;
        for (int i = 0; i < numCorners; ++i)
            deltaKw += fluxVarsCache.epsilonEntryWettingLayerArea(i) * rCdelta * rCdelta / (throatLength * fluxVarsCache.wettingLayerFlowVariables().creviceResistanceFactor(i));

        const auto epsilonPc = fluxVarsCache.epsilonPc();
        auto result = (deltaKw - Kw)/epsilonPc;
        return result;
    }

    /*!
     * \brief Returns derivative of wetting layer transmissibility at entry pressure
     *        For convinience now we calculate the derivative numerically
     *        TODO: Implement the analytical formular of this derivative
     */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar dKwdPcSnapoff(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar Kw = snapoffWettingLayerTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar rCdelta = fluxVarsCache.curvatureRadiusSnapoff(1);

        const auto eIdx = fvGeometry.gridGeometry().elementMapper().index(element);
        const auto shape = fvGeometry.gridGeometry().throatCrossSectionShape(eIdx);
        const auto numCorners = Throat::numCorners(shape);

        // treat the wetting film layer in each corner of the throat individually (might have different corner half-angle and beta)
        Scalar deltaKw = 0.0;
        for (int i = 0; i < numCorners; ++i)
            deltaKw += fluxVarsCache.epsilonSnapoffWettingLayerArea(i) * rCdelta * rCdelta / (throatLength * fluxVarsCache.wettingLayerFlowVariables().creviceResistanceFactor(i));

        const auto epsilonPc = fluxVarsCache.epsilonPc();
        auto result = (deltaKw - Kw)/epsilonPc;
        return result;
    }
};
} // end namespace WettingLayerTransmissibility

namespace NonWettingPhaseTransmissibility {

//! TODO: evalute and maybe remove
template<class Scalar>
struct Mogensen
{
    using NonWettingPhaseCache = EmptyCache;

    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar nonWettingPhaseTransmissibility(const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                  const FluxVariablesCache& fluxVarsCache)
    {
        // Mogensen et al. (1999), does not really revover the single phase value
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const auto nPhaseIdx = fluxVarsCache.nPhaseIdx();
        const Scalar aNw = fluxVarsCache.throatCrossSectionalArea(nPhaseIdx);
        const Scalar rEff = 0.5*(sqrt(aNw / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar result = M_PI/(8*throatLength) * rEff*rEff*rEff*rEff;
        return result;
    }
};

//! TODO: evalute and maybe remove
template<class Scalar, class SinglePhaseTransmissibilityLaw>
struct Valvatne
{
    using NonWettingPhaseCache = EmptyCache;

    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar nonWettingPhaseTransmissibility(const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                  const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), also does not fully recover single-phase value, but is closer
        using std::sqrt;
        const auto nPhaseIdx = fluxVarsCache.nPhaseIdx();
        const Scalar aNw = fluxVarsCache.throatCrossSectionalArea(nPhaseIdx);
        const Scalar aTot = fluxVarsCache.throatCrossSectionalArea();

        const Scalar result = SinglePhaseTransmissibilityLaw::singlePhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache)
                              * aNw / aTot;

        return result;
    }
};

template<class Scalar>
struct BakkeOren
{
    using NonWettingPhaseCache = EmptyCache;

    /*!
    * \brief Returns the conductivity of a throat.
    *
    * See Bakke & Oren (1997), eq. 9
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar nonWettingPhaseTransmissibility(const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                  const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), quite close for single-phase value of square
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const auto nPhaseIdx = fluxVarsCache.nPhaseIdx();
        const Scalar aNw = fluxVarsCache.throatCrossSectionalArea(nPhaseIdx);
        const Scalar rEff = 0.5*(sqrt(aNw / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar result = rEff*rEff*aNw / (8.0*throatLength);
        return result;
    }

    /*!
    * \brief Returns the conductivity of a throat at entry pressure.
    *
    * See Bakke & Oren (1997), eq. 9
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar entryNonWettingPhaseTransmissibility(const Element& element,
                                                       const FVElementGeometry& fvGeometry,
                                                       const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                       const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), quite close for single-phase value of square
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar aNwCrit = fluxVarsCache.regBoundaryNonwettingThroatAreaInvasion(0);
        const Scalar rEff = 0.5 * (sqrt(aNwCrit / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar result = rEff * rEff * aNwCrit / (8.0 * throatLength);
        return result;
    }

    /*!
    * \brief Returns the conductivity of a throat at snap-off pressure.
    *
    * See Bakke & Oren (1997), eq. 9
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar snapoffNonWettingPhaseTransmissibility(const Element& element,
                                                          const FVElementGeometry& fvGeometry,
                                                          const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                                          const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), quite close for single-phase value of square
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar aNwCrit = fluxVarsCache.regBoundaryNonWettingThroatAreaSnapoff(0);
        const Scalar rEff = 0.5*(sqrt(aNwCrit / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar result = rEff*rEff*aNwCrit / (8.0*throatLength);
        return result;
    }

    /*!
    * \brief Returns the conductivity of a throat at entry pressure.
    *
    * See Bakke & Oren (1997), eq. 9
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar dKndPcEntry(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                              const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), quite close for single-phase value of square
        const Scalar Kn = entryNonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar aNwCritDelta = fluxVarsCache.regBoundaryNonwettingThroatAreaInvasion(1);
        const Scalar rEffDelta = 0.5*(sqrt(aNwCritDelta / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar KnDelta = rEffDelta*rEffDelta*aNwCritDelta / (8.0*throatLength);
        const auto epsilonPc = fluxVarsCache.epsilonPc();
        auto result = (KnDelta - Kn)/epsilonPc;
        return result;
    }

    /*!
    * \brief Returns the conductivity of a throat at snap-off pressure.
    *
    * See Bakke & Oren (1997), eq. 9
    */
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    static Scalar dKndPcSnapoff(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                const FluxVariablesCache& fluxVarsCache)
    {
        // Tora et al. (2012), quite close for single-phase value of square
        const Scalar Kn = snapoffNonWettingPhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
        using std::sqrt;
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar aNwCritDelta = fluxVarsCache.regBoundaryNonWettingThroatAreaSnapoff(1);
        const Scalar rEffDelta = 0.5*(sqrt(aNwCritDelta / M_PI) + fluxVarsCache.throatInscribedRadius());
        const Scalar KnDelta = rEffDelta*rEffDelta*aNwCritDelta / (8.0*throatLength);
        const auto epsilonPc = fluxVarsCache.epsilonPc();
        auto result = (KnDelta - Kn)/epsilonPc;
        return result;
    }
};
} // end namespace NonWettingPhaseTransmissibility
} // end namespace Dumux::PoreNetwork

#endif
