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
 * \ingroup NavierStokesNCTests
 * \brief Density driven flow test for the multi-component staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_DENSITY_DRIVEN_NC_TEST_PROBLEM_HH
#define DUMUX_DENSITY_DRIVEN_NC_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the one-phase (Navier-)Stokes model.
 *
 * Density driven flow test for the multi-component staggered grid (Navier-)Stokes model.
 */
template <class TypeTag>
class DensityDrivenFlowProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto compIdx = 1;

public:
    DensityDrivenFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    , eps_(1e-6)
    {
        useWholeLength_ = getParam<bool>("Problem.UseWholeLength");
    }

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related to floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.1e5; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet();
        else
        {
            values.setAllNeumann();
            if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
            {
                if (useWholeLength_ || (globalPos[0] > 0.4 && globalPos[0] < 0.6))
                    values.setAllDirichlet();
            }

        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        PrimaryVariables values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            values[Indices::pressureIdx] = 1.1e+5;
            if (useWholeLength_ || (globalPos[0] > 0.4 && globalPos[0] < 0.6))
                values[Indices::conti0EqIdx + compIdx] = 1e-3;
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();

            // The resulting flux over the boundary is zero anyway (velocity is zero), but this will add some non-zero derivatives to the
            // Jacobian and makes the BC more general.
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
            values[Indices::pressureIdx] = 1.1e5;

        return values;
    }

    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    Scalar time() const
    { return timeLoop_->time(); }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<PrimaryVariables::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> values;

        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof

        if constexpr (!ParentType::isMomentumProblem())
        {
            const bool isLowerLeftCell = (scv.dofIndex() == 0);
            if (isLowerLeftCell)
                values.set(0);
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(1.1e5); }

private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    const Scalar eps_;
    bool useWholeLength_;
    TimeLoopPtr timeLoop_;
};

} // end namespace Dumux
#endif
