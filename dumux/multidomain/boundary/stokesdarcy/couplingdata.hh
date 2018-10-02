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
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \copydoc Dumux::StokesDarcyCouplingData
 */

#ifndef DUMUX_STOKES_DARCY_COUPLINGDATA_HH
#define DUMUX_STOKES_DARCY_COUPLINGDATA_HH

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/multidomain/couplingmanager.hh>

namespace Dumux {


/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief This structs holds a set of options which allow to modify the Stokes-Darcy coupling mechanism during runtime.
 */
struct StokesDarcyCouplingOptions
{
    /*!
     * \brief Defines which kind of averanging of diffusion coefficiencients (moleculat diffusion or thermal conductance)
     *        at the interface between free flow and porous medium shall be used.
     */
    enum class DiffusionCoefficientAveragingType
    {
        harmonic, arithmethic, ffOnly, pmOnly
    };

    /*!
     * \brief Convenience function to convert user input given as std::string to the corresponding enum class used for chosing the type
     *        of averaging of the diffusion/conduction parameter at the interface between the two domains.
     */
    static DiffusionCoefficientAveragingType stringToEnum(DiffusionCoefficientAveragingType, const std::string& diffusionCoefficientAveragingType)
    {
        if (diffusionCoefficientAveragingType == "Harmonic")
            return DiffusionCoefficientAveragingType::harmonic;
        else if (diffusionCoefficientAveragingType == "Arithmethic")
            return DiffusionCoefficientAveragingType::arithmethic;
        else if (diffusionCoefficientAveragingType == "FreeFlowOnly")
            return DiffusionCoefficientAveragingType::ffOnly;
        else if (diffusionCoefficientAveragingType == "PorousMediumOnly")
            return DiffusionCoefficientAveragingType::pmOnly;
        else
            DUNE_THROW(Dune::IOError, "Unknown DiffusionCoefficientAveragingType");
    }

};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 *        Specialization for the case of using an adapter only for the free-flow model.
 * \tparam FFFS The free-flow fluidsystem
 * \tparam PMFS The porous-medium flow fluidsystem
 */
template<class FFFS, class PMFS>
struct IsSameFluidSystem
{
    static_assert(FFFS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<typename FFFS::MultiPhaseFluidSystem, PMFS>::value;
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief This structs helps to check if the two sub models use the same fluidsystem.
 * \tparam FS The fluidsystem
 */
template<class FS>
struct IsSameFluidSystem<FS, FS>
{
    static_assert(FS::numPhases == 1, "Only single-phase fluidsystems may be used for free flow.");
    static constexpr bool value = std::is_same<FS, FS>::value; // always true
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 * \tparam hasAdapter Specifies whether an adapter class for the fluidsystem is used.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS, bool hasAdapter>
struct IndexHelper;

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that no adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, darcyIdx, FFFS, false>
{
    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i>, int coupledPhaseIdx = 0)
    { return coupledPhaseIdx; }

    /*!
     * \brief No adapter is used, just return the input index.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i>, int coupledCompdIdx)
    { return coupledCompdIdx; }
};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Helper struct to choose the correct index for phases and components. This is need if the porous-medium-flow model
          features more fluid phases than the free-flow model. Specialization for the case that a adapter is used.
 * \tparam stokesIdx The domain index of the free-flow model.
 * \tparam darcyIdx The domain index of the porous-medium-flow model.
 * \tparam FFFS The free-flow fluidsystem.
 */
template<std::size_t stokesIdx, std::size_t darcyIdx, class FFFS>
struct IndexHelper<stokesIdx, darcyIdx, FFFS, true>
{
    /*!
     * \brief The free-flow model always uses phase index 0.
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<stokesIdx>, int coupledPhaseIdx = 0)
    { return 0; }

    /*!
     * \brief The phase index of the porous-medium-flow model is given by the adapter fluidsytem (i.e., user input).
     */
    static constexpr auto couplingPhaseIdx(Dune::index_constant<darcyIdx>, int coupledPhaseIdx = 0)
    { return FFFS::multiphaseFluidsystemPhaseIdx; }

    /*!
     * \brief The free-flow model does not need any change of the component index.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<stokesIdx>, int coupledCompdIdx)
    { return coupledCompdIdx; }

    /*!
     * \brief The component index of the porous-medium-flow model is mapped by the adapter fluidsytem.
     */
    static constexpr auto couplingCompIdx(Dune::index_constant<darcyIdx>, int coupledCompdIdx)
    { return FFFS::compIdx(coupledCompdIdx); }
};

template<class MDTraits, class CouplingManager, bool enableEnergyBalance, bool isCompositional>
class StokesDarcyCouplingDataImplementation;

/*!
* \ingroup MultiDomain
* \ingroup BoundaryCoupling
* \brief Data for the coupling of a Darcy model (cell-centered finite volume)
*        with a (Navier-)Stokes model (staggerd grid).
*/
template<class MDTraits, class CouplingManager>
using StokesDarcyCouplingData = StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager,
                                                                      GET_PROP_TYPE(typename MDTraits::template SubDomainTypeTag<0>, ModelTraits)::enableEnergyBalance(),
                                                                      (GET_PROP_TYPE(typename MDTraits::template SubDomainTypeTag<0>, ModelTraits)::numComponents() > 1)>;

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief A base class which provides some common methods used for Stokes-Darcy coupling.
 */
template<class MDTraits, class CouplingManager>
class StokesDarcyCouplingDataImplementationBase
{
    using Scalar = typename MDTraits::Scalar;

    template<std::size_t id> using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;
    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ModelTraits)::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::VolumeVariables;
    template<std::size_t id> using Problem  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, Problem);
    template<std::size_t id> using FluidSystem  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FluidSystem);
    template<std::size_t id> using ModelTraits  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ModelTraits);

    static constexpr auto stokesIdx = CouplingManager::stokesIdx;
    static constexpr auto darcyIdx = CouplingManager::darcyIdx;

    static constexpr bool adapterUsed = ModelTraits<darcyIdx>::numPhases() > 1;
    using IndexHelper = Dumux::IndexHelper<stokesIdx, darcyIdx, FluidSystem<stokesIdx>, adapterUsed>;

    static constexpr int enableEnergyBalance = GET_PROP_TYPE(SubDomainTypeTag<stokesIdx>, ModelTraits)::enableEnergyBalance();
    static_assert(GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ModelTraits)::enableEnergyBalance() == enableEnergyBalance,
                  "All submodels must both be either isothermal or non-isothermal");

    static_assert(IsSameFluidSystem<FluidSystem<stokesIdx>,
                                    FluidSystem<darcyIdx>>::value,
                  "All submodels must use the same fluid system");

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    StokesDarcyCouplingDataImplementationBase(const CouplingManager& couplingmanager): couplingManager_(couplingmanager) {}

    /*!
     * \brief Returns the corresponding phase index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingPhaseIdx(Dune::index_constant<i> id, int coupledPhaseIdx = 0)
    { return IndexHelper::couplingPhaseIdx(id, coupledPhaseIdx); }

    /*!
     * \brief Returns the corresponding component index needed for coupling.
     */
    template<std::size_t i>
    static constexpr auto couplingCompIdx(Dune::index_constant<i> id, int coupledCompdIdx)
    { return IndexHelper::couplingCompIdx(id, coupledCompdIdx); }

    /*!
     * \brief Returns a reference to the coupling manager.
     */
    const CouplingManager& couplingManager() const
    { return couplingManager_; }

    /*!
     * \brief Returns the intrinsic permeability of the coupled Darcy element.
     */
    Scalar darcyPermeability(const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        const auto& stokesContext = couplingManager().stokesCouplingContext(scvf);
        return stokesContext.volVars.permeability();
    }

     /*!
     * \brief Returns the momentum flux across the coupling boundary.
     *
     * For the normal momentum coupling, the porous medium side of the coupling condition
     * is evaluated, i.e. -[p n]^pm.
     *
     */
    template<class ElementFaceVariables>
    Scalar momentumCouplingCondition(const FVElementGeometry<stokesIdx>& fvGeometry,
                                     const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                     const ElementFaceVariables& stokesElemFaceVars,
                                     const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        static constexpr auto numPhasesDarcy = GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ModelTraits)::numPhases();

        Scalar momentumFlux(0.0);
        const auto& stokesContext = couplingManager_.stokesCouplingContext(scvf);
        const auto darcyPhaseIdx = couplingPhaseIdx(darcyIdx);

        // - p_pm * n_pm = p_pm * n_ff
        const Scalar darcyPressure = stokesContext.volVars.pressure(darcyPhaseIdx);

        if(numPhasesDarcy > 1)
        {
            momentumFlux = darcyPressure;
        }
        else // use pressure reconstruction for single phase models
        {
            // v = -K/mu * (gradP + rho*g)
            const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
            const Scalar mu = stokesContext.volVars.viscosity(darcyPhaseIdx);
            const Scalar rho = stokesContext.volVars.density(darcyPhaseIdx);
            const Scalar distance = (stokesContext.element.geometry().center() - scvf.center()).two_norm();
            const Scalar g = -scvf.directionSign() * couplingManager_.problem(darcyIdx).gravity()[scvf.directionIndex()];
            const Scalar interfacePressure = ((scvf.directionSign() * velocity * (mu/darcyPermeability(scvf))) + rho * g) * distance + darcyPressure;
            momentumFlux = interfacePressure;
        }

        // normalize pressure
        if(GET_PROP_VALUE(SubDomainTypeTag<stokesIdx>, NormalizePressure))
            momentumFlux -= couplingManager_.problem(stokesIdx).initial(scvf)[Indices<stokesIdx>::pressureIdx];

        momentumFlux *= scvf.directionSign();

        return momentumFlux;
    }

    /*!
     * \brief Evaluate an advective flux across the interface and consider upwinding.
     */
    Scalar advectiveFlux(const Scalar insideQuantity, const Scalar outsideQuantity, const Scalar volumeFlow, bool insideIsUpstream) const
    {
        const Scalar upwindWeight = 1.0; //TODO use Implicit.UpwindWeight or something like Coupling.UpwindWeight?

        if(insideIsUpstream)
            return (upwindWeight * insideQuantity + (1.0 - upwindWeight) * outsideQuantity) * volumeFlow;
        else
            return (upwindWeight * outsideQuantity + (1.0 - upwindWeight) * insideQuantity) * volumeFlow;
    }

protected:

    /*!
     * \brief Returns the transmissibility used for either molecular diffusion or thermal conductivity.
     */
    template<std::size_t i, std::size_t j>
    Scalar transmissibility_(Dune::index_constant<i> domainI,
                             Dune::index_constant<j> domainJ,
                             const Scalar insideDistance,
                             const Scalar outsideDistance,
                             const Scalar avgQuantityI,
                             const Scalar avgQuantityJ,
                             const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar totalDistance = insideDistance + outsideDistance;
        if(diffCoeffAvgType == DiffusionCoefficientAveragingType::harmonic)
        {
            return harmonicMean(avgQuantityI, avgQuantityJ, insideDistance, outsideDistance)
                   / totalDistance;
        }
        else if(diffCoeffAvgType == DiffusionCoefficientAveragingType::arithmethic)
        {
            return arithmeticMean(avgQuantityI, avgQuantityJ, insideDistance, outsideDistance)
                   / totalDistance;
        }
        else if(diffCoeffAvgType == DiffusionCoefficientAveragingType::ffOnly)
            return domainI == stokesIdx
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;
        else // diffCoeffAvgType == DiffusionCoefficientAveragingType::pmOnly)
            return domainI == darcyIdx
                            ? avgQuantityI / totalDistance
                            : avgQuantityJ / totalDistance;
    }

    /*!
     * \brief Returns the distance between an scvf and the corresponding scv center.
     */
    template<class Scv, class Scvf>
    Scalar getDistance_(const Scv& scv, const Scvf& scvf) const
    {
        return (scv.dofPosition() - scvf.ipGlobal()).two_norm();
    }

    /*!
     * \brief Returns the conductive energy flux acorss the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar conductiveEnergyFlux_(Dune::index_constant<i> domainI,
                                 Dune::index_constant<j> domainJ,
                                 const FVElementGeometry<i>& fvGeometryI,
                                 const FVElementGeometry<j>& fvGeometryJ,
                                 const SubControlVolumeFace<i>& scvfI,
                                 const SubControlVolume<i>& scvI,
                                 const SubControlVolume<j>& scvJ,
                                 const VolumeVariables<i>& volVarsI,
                                 const VolumeVariables<j>& volVarsJ,
                                 const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        const Scalar insideDistance = getDistance_(scvI, scvfI);
        const Scalar outsideDistance = getDistance_(scvJ, scvfI);

        const Scalar deltaT = volVarsJ.temperature() - volVarsI.temperature();
        const Scalar tij = transmissibility_(domainI,
                                             domainJ,
                                             insideDistance,
                                             outsideDistance,
                                             thermalConductivity_(volVarsI, fvGeometryI, scvI) * volVarsI.extrusionFactor(),
                                             thermalConductivity_(volVarsJ, fvGeometryJ, scvJ) * volVarsJ.extrusionFactor(),
                                             diffCoeffAvgType);

        return -tij * deltaT;
    }

    /*!
     * \brief Returns the effective thermal conductivity (lumped parameter) within the porous medium.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar thermalConductivity_(const VolumeVariables<darcyIdx>& volVars,
                                const FVElementGeometry<darcyIdx>& fvGeometry,
                                const SubControlVolume<darcyIdx>& scv) const
    {
        using ThermalConductivityModel = typename GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ThermalConductivityModel);
        const auto& problem = this->couplingManager().problem(darcyIdx);
        return ThermalConductivityModel::effectiveThermalConductivity(volVars, problem.spatialParams(), fvGeometry.fvGridGeometry().element(scv), fvGeometry, scv);
    }

    /*!
     * \brief Returns the thermal conductivity of the fluid phase within the free flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar thermalConductivity_(const VolumeVariables<stokesIdx>& volVars,
                                const FVElementGeometry<stokesIdx>& fvGeometry,
                                const SubControlVolume<stokesIdx>& scv) const
    {
        return  volVars.effectiveThermalConductivity();
    }

private:
    const CouplingManager& couplingManager_;

};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for non-compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, false>
: public StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto stokesIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto darcyIdx = typename MDTraits::template DomainIdx<2>();
    static constexpr auto stokesCellCenterIdx = stokesIdx;
    static constexpr auto stokesFaceIdx = typename MDTraits::template DomainIdx<1>();

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVGridGeometry<id>::LocalView::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ModelTraits)::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridFaceVariables)::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::VolumeVariables;

    static_assert(GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ModelTraits)::numComponents() == GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ModelTraits)::numPhases(),
                  "Darcy Model must not be compositional");

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    Scalar massCouplingCondition(const FVElementGeometry<darcyIdx>& fvGeometry,
                                 const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                 const SubControlVolumeFace<darcyIdx>& scvf) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(scvf);
        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const Scalar darcyDensity = darcyElemVolVars[scvf.insideScvIdx()].density(couplingPhaseIdx(darcyIdx));
        const Scalar stokesDensity = darcyContext.volVars.density();
        const bool insideIsUpstream = velocity > 0.0;

        return massFlux_(velocity, darcyDensity, stokesDensity, insideIsUpstream);
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    Scalar massCouplingCondition(const FVElementGeometry<stokesIdx>& fvGeometry,
                                 const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                 const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                 const SubControlVolumeFace<stokesIdx>& scvf) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(scvf);
        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const Scalar stokesDensity = stokesElemVolVars[scvf.insideScvIdx()].density();
        const Scalar darcyDensity = stokesContext.volVars.density(couplingPhaseIdx(darcyIdx));
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return massFlux_(velocity * scvf.directionSign(), stokesDensity, darcyDensity, insideIsUpstream);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const FVElementGeometry<darcyIdx>& fvGeometry,
                                   const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                   const SubControlVolumeFace<darcyIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return energyFlux_(darcyIdx, stokesIdx, fvGeometry, darcyContext.fvGeometry, scvf,
                           darcyVolVars, stokesVolVars, velocity, insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const FVElementGeometry<stokesIdx>& fvGeometry,
                                   const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<stokesIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return energyFlux_(stokesIdx, darcyIdx, fvGeometry, stokesContext.fvGeometry, scvf,
                           stokesVolVars, darcyVolVars, velocity * scvf.directionSign(), insideIsUpstream, diffCoeffAvgType);
    }

private:

    /*!
     * \brief Evaluate the mole/mass flux across the interface.
     */
    Scalar massFlux_(const Scalar velocity,
                     const Scalar insideDensity,
                     const Scalar outSideDensity,
                     bool insideIsUpstream) const
    {
        return this->advectiveFlux(insideDensity, outSideDensity, velocity, insideIsUpstream);
    }

    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(Dune::index_constant<i> domainI,
                       Dune::index_constant<j> domainJ,
                       const FVElementGeometry<i>& insideFvGeometry,
                       const FVElementGeometry<j>& outsideFvGeometry,
                       const SubControlVolumeFace<i>& scvf,
                       const VolumeVariables<i>& insideVolVars,
                       const VolumeVariables<j>& outsideVolVars,
                       const Scalar velocity,
                       const bool insideIsUpstream,
                       const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const auto& insideScv = (*scvs(insideFvGeometry).begin());
        const auto& outsideScv = (*scvs(outsideFvGeometry).begin());

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        flux += this->advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        flux += this->conductiveEnergyFlux_(domainI, domainJ, insideFvGeometry, outsideFvGeometry, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);

        return flux;
    }

};

/*!
 * \ingroup MultiDomain
 * \ingroup BoundaryCoupling
 * \ingroup StokesDarcyCoupling
 * \brief Coupling data specialization for compositional models.
 */
template<class MDTraits, class CouplingManager, bool enableEnergyBalance>
class StokesDarcyCouplingDataImplementation<MDTraits, CouplingManager, enableEnergyBalance, true>
: public StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>
{
    using ParentType = StokesDarcyCouplingDataImplementationBase<MDTraits, CouplingManager>;
    using Scalar = typename MDTraits::Scalar;
    static constexpr auto stokesIdx = typename MDTraits::template DomainIdx<0>();
    static constexpr auto darcyIdx = typename MDTraits::template DomainIdx<2>();
    static constexpr auto stokesCellCenterIdx = stokesIdx;
    static constexpr auto stokesFaceIdx = typename MDTraits::template DomainIdx<1>();

    // the sub domain type tags
    template<std::size_t id>
    using SubDomainTypeTag = typename MDTraits::template SubDomainTypeTag<id>;

    template<std::size_t id> using FVGridGeometry = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FVGridGeometry);
    template<std::size_t id> using FVElementGeometry = typename FVGridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using SubControlVolume = typename FVGridGeometry<id>::LocalView::SubControlVolume;
    template<std::size_t id> using Indices = typename GET_PROP_TYPE(SubDomainTypeTag<id>, ModelTraits)::Indices;
    template<std::size_t id> using ElementVolumeVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::LocalView;
    template<std::size_t id> using ElementFaceVariables = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridFaceVariables)::LocalView;
    template<std::size_t id> using VolumeVariables  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, GridVolumeVariables)::VolumeVariables;
    template<std::size_t id> using FluidSystem  = typename GET_PROP_TYPE(SubDomainTypeTag<id>, FluidSystem);

    static constexpr auto numComponents = GET_PROP_TYPE(SubDomainTypeTag<stokesIdx>, ModelTraits)::numComponents();
    static constexpr auto replaceCompEqIdx = GET_PROP_TYPE(SubDomainTypeTag<stokesIdx>, ModelTraits)::replaceCompEqIdx();
    static constexpr bool useMoles = GET_PROP_TYPE(SubDomainTypeTag<stokesIdx>, ModelTraits)::useMoles();

    static_assert(GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, ModelTraits)::numComponents() == numComponents, "Submodels must use same number of components");
    static_assert(GET_PROP_VALUE(SubDomainTypeTag<darcyIdx>, UseMoles) == useMoles, "Both models must either use moles or not");
    static_assert(GET_PROP_VALUE(SubDomainTypeTag<darcyIdx>, ReplaceCompEqIdx) == replaceCompEqIdx, "Both models must use the same replaceCompEqIdx");
    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

public:
    using ParentType::ParentType;
    using ParentType::couplingPhaseIdx;
    using ParentType::couplingCompIdx;

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the Darcy domain.
     */
    NumEqVector massCouplingCondition(const FVElementGeometry<darcyIdx>& fvGeometry,
                                      const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                      const SubControlVolumeFace<darcyIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        NumEqVector flux(0.0);
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;
        const auto& outsideScv = (*scvs(darcyContext.fvGeometry).begin());

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return massFlux_(darcyIdx, stokesIdx, fvGeometry,
                         scvf, darcyVolVars, stokesVolVars,
                         outsideScv, velocity, insideIsUpstream,
                         diffCoeffAvgType);
    }

    /*!
     * \brief Returns the mass flux across the coupling boundary as seen from the free-flow domain.
     */
    NumEqVector massCouplingCondition(const FVElementGeometry<stokesIdx>& fvGeometry,
                                      const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                      const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                      const SubControlVolumeFace<stokesIdx>& scvf,
                                      const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        NumEqVector flux(0.0);
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;
        const auto& outsideScv = (*scvs(stokesContext.fvGeometry).begin());

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return massFlux_(stokesIdx, darcyIdx, fvGeometry,
                         scvf, stokesVolVars, darcyVolVars,
                         outsideScv, velocity * scvf.directionSign(),
                         insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the Darcy domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const FVElementGeometry<darcyIdx>& fvGeometry,
                                   const ElementVolumeVariables<darcyIdx>& darcyElemVolVars,
                                   const SubControlVolumeFace<darcyIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& darcyContext = this->couplingManager().darcyCouplingContext(scvf);
        const auto& darcyVolVars = darcyElemVolVars[scvf.insideScvIdx()];
        const auto& stokesVolVars = darcyContext.volVars;

        const Scalar velocity = darcyContext.velocity * scvf.unitOuterNormal();
        const bool insideIsUpstream = velocity > 0.0;

        return energyFlux_(darcyIdx, stokesIdx, fvGeometry, darcyContext.fvGeometry, scvf,
                           darcyVolVars, stokesVolVars, velocity, insideIsUpstream, diffCoeffAvgType);
    }

    /*!
     * \brief Returns the energy flux across the coupling boundary as seen from the free-flow domain.
     */
    template<bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyCouplingCondition(const FVElementGeometry<stokesIdx>& fvGeometry,
                                   const ElementVolumeVariables<stokesIdx>& stokesElemVolVars,
                                   const ElementFaceVariables<stokesIdx>& stokesElemFaceVars,
                                   const SubControlVolumeFace<stokesIdx>& scvf,
                                   const DiffusionCoefficientAveragingType diffCoeffAvgType = DiffusionCoefficientAveragingType::ffOnly) const
    {
        const auto& stokesContext = this->couplingManager().stokesCouplingContext(scvf);
        const auto& stokesVolVars = stokesElemVolVars[scvf.insideScvIdx()];
        const auto& darcyVolVars = stokesContext.volVars;

        const Scalar velocity = stokesElemFaceVars[scvf].velocitySelf();
        const bool insideIsUpstream = sign(velocity) == scvf.directionSign();

        return energyFlux_(stokesIdx, darcyIdx, fvGeometry, stokesContext.fvGeometry, scvf,
                           stokesVolVars, darcyVolVars, velocity * scvf.directionSign(), insideIsUpstream, diffCoeffAvgType);
    }

protected:

    /*!
     * \brief Evaluate the compositional mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector massFlux_(Dune::index_constant<i> domainI,
                          Dune::index_constant<j> domainJ,
                          const FVElementGeometry<i>& insideFvGeometry,
                          const SubControlVolumeFace<i>& scvf,
                          const VolumeVariables<i>& insideVolVars,
                          const VolumeVariables<j>& outsideVolVars,
                          const SubControlVolume<j>& outsideScv,
                          const Scalar velocity,
                          const bool insideIsUpstream,
                          const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector flux(0.0);

        auto moleOrMassFraction = [&](const auto& volVars, int phaseIdx, int compIdx)
        { return useMoles ? volVars.moleFraction(phaseIdx, compIdx) : volVars.massFraction(phaseIdx, compIdx); };

        auto moleOrMassDensity = [&](const auto& volVars, int phaseIdx)
        { return useMoles ? volVars.molarDensity(phaseIdx) : volVars.density(phaseIdx); };

        // treat the advective fluxes
        auto insideTerm = [&](int compIdx)
        { return moleOrMassFraction(insideVolVars, couplingPhaseIdx(domainI), compIdx) * moleOrMassDensity(insideVolVars, couplingPhaseIdx(domainI)); };

        auto outsideTerm = [&](int compIdx)
        { return moleOrMassFraction(outsideVolVars, couplingPhaseIdx(domainJ), compIdx) * moleOrMassDensity(outsideVolVars, couplingPhaseIdx(domainJ)); };


        for(int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);
            flux[domainICompIdx] += this->advectiveFlux(insideTerm(domainICompIdx), outsideTerm(domainJCompIdx), velocity, insideIsUpstream);
        }

        // treat the diffusive fluxes
        const auto& insideScv = insideFvGeometry.scv(scvf.insideScvIdx());
        auto diffusiveFlux = diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);

        if(!useMoles)
        {
            //convert everything to a mass flux
            for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                diffusiveFlux[compIdx] *= FluidSystem<i>::molarMass(compIdx);
        }

        flux += diffusiveFlux;

        // convert to total mass/mole balance, if set be user
        if(replaceCompEqIdx < numComponents)
            flux[replaceCompEqIdx] = std::accumulate(flux.begin(), flux.end(), 0.0);

        return flux;
    }

    /*!
     * \brief Returns the molecular diffusion coefficient within the free flow domain.
     */
    Scalar diffusionCoefficient_(const VolumeVariables<stokesIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return volVars.effectiveDiffusivity(phaseIdx, compIdx);
    }

    /*!
     * \brief Returns the effective diffusion coefficient within the porous medium.
     */
    Scalar diffusionCoefficient_(const VolumeVariables<darcyIdx>& volVars, int phaseIdx, int compIdx) const
    {
        using EffDiffModel = typename GET_PROP_TYPE(SubDomainTypeTag<darcyIdx>, EffectiveDiffusivityModel);
        return EffDiffModel::effectiveDiffusivity(volVars.porosity(),
                                                  volVars.saturation(phaseIdx),
                                                  volVars.diffusionCoefficient(phaseIdx, compIdx));
    }

    Scalar getComponentEnthalpy(const VolumeVariables<stokesIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<stokesIdx>::componentEnthalpy(volVars.fluidState(), 0, compIdx);
    }

    Scalar getComponentEnthalpy(const VolumeVariables<darcyIdx>& volVars, int phaseIdx, int compIdx) const
    {
        return FluidSystem<darcyIdx>::componentEnthalpy(volVars.fluidState(), phaseIdx, compIdx);
    }

    /*!
     * \brief Evaluate the diffusive mole/mass flux across the interface.
     */
    template<std::size_t i, std::size_t j>
    NumEqVector diffusiveMolecularFlux_(Dune::index_constant<i> domainI,
                                        Dune::index_constant<j> domainJ,
                                        const SubControlVolumeFace<i>& scvfI,
                                        const SubControlVolume<i>& scvI,
                                        const SubControlVolume<j>& scvJ,
                                        const VolumeVariables<i>& volVarsI,
                                        const VolumeVariables<j>& volVarsJ,
                                        const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        NumEqVector diffusiveFlux(0.0);
        const Scalar avgMolarDensity = 0.5 * volVarsI.molarDensity(couplingPhaseIdx(domainI)) + 0.5 *  volVarsJ.molarDensity(couplingPhaseIdx(domainJ));
        const Scalar insideDistance = this->getDistance_(scvI, scvfI);
        const Scalar outsideDistance = this->getDistance_(scvJ, scvfI);

        for(int compIdx = 1; compIdx < numComponents; ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            assert(FluidSystem<i>::componentName(domainICompIdx) == FluidSystem<j>::componentName(domainJCompIdx));

            const Scalar deltaMoleFrac = volVarsJ.moleFraction(couplingPhaseIdx(domainJ), domainJCompIdx) - volVarsI.moleFraction(couplingPhaseIdx(domainI), domainICompIdx);
            const Scalar tij = this->transmissibility_(domainI,
                                                       domainJ,
                                                       insideDistance,
                                                       outsideDistance,
                                                       diffusionCoefficient_(volVarsI, couplingPhaseIdx(domainI), domainICompIdx) * volVarsI.extrusionFactor(),
                                                       diffusionCoefficient_(volVarsJ, couplingPhaseIdx(domainJ), domainJCompIdx) * volVarsJ.extrusionFactor(),
                                                       diffCoeffAvgType);
            diffusiveFlux[domainICompIdx] += -avgMolarDensity * tij * deltaMoleFrac;
        }

        const Scalar cumulativeFlux = std::accumulate(diffusiveFlux.begin(), diffusiveFlux.end(), 0.0);
        diffusiveFlux[couplingCompIdx(domainI, 0)] = -cumulativeFlux;

        return diffusiveFlux;
    }

    /*!
     * \brief Evaluate the energy flux across the interface.
     */
    template<std::size_t i, std::size_t j, bool isNI = enableEnergyBalance, typename std::enable_if_t<isNI, int> = 0>
    Scalar energyFlux_(Dune::index_constant<i> domainI,
                       Dune::index_constant<j> domainJ,
                       const FVElementGeometry<i>& insideFvGeometry,
                       const FVElementGeometry<j>& outsideFvGeometry,
                       const SubControlVolumeFace<i>& scvf,
                       const VolumeVariables<i>& insideVolVars,
                       const VolumeVariables<j>& outsideVolVars,
                       const Scalar velocity,
                       const bool insideIsUpstream,
                       const DiffusionCoefficientAveragingType diffCoeffAvgType) const
    {
        Scalar flux(0.0);

        const auto& insideScv = (*scvs(insideFvGeometry).begin());
        const auto& outsideScv = (*scvs(outsideFvGeometry).begin());

        // convective fluxes
        const Scalar insideTerm = insideVolVars.density(couplingPhaseIdx(domainI)) * insideVolVars.enthalpy(couplingPhaseIdx(domainI));
        const Scalar outsideTerm = outsideVolVars.density(couplingPhaseIdx(domainJ)) * outsideVolVars.enthalpy(couplingPhaseIdx(domainJ));

        flux += this->advectiveFlux(insideTerm, outsideTerm, velocity, insideIsUpstream);

        flux += this->conductiveEnergyFlux_(domainI, domainJ, insideFvGeometry, outsideFvGeometry, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);

        auto diffusiveFlux = this->diffusiveMolecularFlux_(domainI, domainJ, scvf, insideScv, outsideScv, insideVolVars, outsideVolVars, diffCoeffAvgType);

        for (int compIdx = 0; compIdx < diffusiveFlux.size(); ++compIdx)
        {
            const int domainICompIdx = couplingCompIdx(domainI, compIdx);
            const int domainJCompIdx = couplingCompIdx(domainJ, compIdx);

            const bool insideDiffFluxIsUpstream = std::signbit(diffusiveFlux[domainICompIdx]);
            const Scalar componentEnthalpy = insideDiffFluxIsUpstream ?
                                             getComponentEnthalpy(insideVolVars, couplingPhaseIdx(domainI), domainICompIdx)
                                           : getComponentEnthalpy(outsideVolVars, couplingPhaseIdx(domainJ), domainJCompIdx);

            // always use a mass-based calculation for the energy balance
            diffusiveFlux[domainICompIdx] *= FluidSystem<i>::molarMass(domainICompIdx);

            flux += diffusiveFlux[domainICompIdx] * componentEnthalpy;
        }

        return flux;
    }
};

} // end namespace Dumux

#endif // DUMUX_STOKES_DARCY_COUPLINGDATA_HH
