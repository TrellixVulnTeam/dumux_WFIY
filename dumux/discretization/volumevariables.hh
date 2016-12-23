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
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_DISCRETIZATION_VOLUME_VARIABLES_HH
#define DUMUX_DISCRETIZATION_VOLUME_VARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(Indices);
}

// forward declaration
template <class TypeTag, bool enableEnergyBalance>
class ImplicitVolumeVariablesImplementation;

/*!
 * \ingroup ImplicitVolumeVariables
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities. The volume variables base class
 *        is specialized for isothermal and non-isothermal models.
 */
template <class TypeTag>
using ImplicitVolumeVariables = ImplicitVolumeVariablesImplementation<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup ImplicitVolumeVariables
 * \brief The isothermal base class
 */
template<class TypeTag>
class ImplicitVolumeVariablesImplementation<TypeTag, false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);

public:

    using PermeabilityType = decltype(std::declval<SpatialParams>().permeability(Element(),
                                                                                 SubControlVolume(),
                                                                                 ElementSolutionVector()));

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub-control volume
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        priVars_ = isBox ? elemSol[scv.index()] : elemSol[0];
        extrusionFactor_ = problem.extrusionFactor(element, scv, elemSol);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &priVars() const
    { return priVars_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    //! The temperature is obtained from the problem as a constant for isothermal models
    static Scalar temperature(const ElementSolutionVector &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const SubControlVolume &scv)
    {
        return problem.temperatureAtPos(scv.dofPosition());
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int phaseIdx)
    {
        return 0;
    }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

//! The non-isothermal implicit volume variables base class
template <class TypeTag>
class ImplicitVolumeVariablesImplementation<TypeTag, true>
: public ImplicitVolumeVariablesImplementation<TypeTag, false>
{
    using ParentType = ImplicitVolumeVariablesImplementation<TypeTag, false>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int temperatureIdx = Indices::temperatureIdx;
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);

public:

    using typename ParentType::PermeabilityType;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param priVars A vector containing the primary variables for the control volume
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite volume geometry for the element
     * \param scvIdx Local index of the sub control volume which is inside the element
     */
    void update(const ElementSolutionVector &elemSol,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        solidHeatCapacity_ = problem.spatialParams().solidHeatCapacity(element, scv, elemSol);
        solidDensity_ = problem.spatialParams().solidDensity(element, scv, elemSol);
        solidThermalConductivity_ = problem.spatialParams().solidThermalConductivity(element, scv, elemSol);
    }

    /*!
     * \brief Returns the total internal energy of a phase in the
     *        sub-control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar internalEnergy(const int phaseIdx) const
    { return asImp_().fluidState().internalEnergy(phaseIdx); }

    /*!
     * \brief Returns the total enthalpy of a phase in the sub-control
     *        volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar enthalpy(const int phaseIdx) const
    { return asImp_().fluidState().enthalpy(phaseIdx); }

    /*!
     * \brief Returns the total heat capacity \f$\mathrm{[J/(kg K)]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidHeatCapacity() const
    { return solidHeatCapacity_; }

    /*!
     * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the rock matrix in
     *        the sub-control volume.
     */
    Scalar solidDensity() const
    { return solidDensity_; }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of a fluid phase in
     *        the sub-control volume.
     */
    Scalar fluidThermalConductivity(const int phaseIdx) const
    { return FluidSystem::thermalConductivity(asImp_().fluidState(), phaseIdx); }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m*K)]}\f$ of the solid phase in
     *        the sub-control volume.
     */
    Scalar solidThermalConductivity() const
    { return solidThermalConductivity_; }

    //! The temperature is a primary variable for non-isothermal models
    static Scalar temperature(const ElementSolutionVector &elemSol,
                              const Problem& problem,
                              const Element &element,
                              const SubControlVolume &scv)
    {
        return isBox ? elemSol[scv.index()][temperatureIdx] : elemSol[0][temperatureIdx];
    }

    //! The phase enthalpy is zero for isothermal models
    //! This is needed for completing the fluid state
    template<class FluidState, class ParameterCache>
    static Scalar enthalpy(const FluidState& fluidState,
                           const ParameterCache& paramCache,
                           const int phaseIdx)
    {
        return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

private:
    Scalar solidHeatCapacity_;
    Scalar solidDensity_;
    Scalar solidThermalConductivity_;
};

} // end namespace Dumux

#endif
