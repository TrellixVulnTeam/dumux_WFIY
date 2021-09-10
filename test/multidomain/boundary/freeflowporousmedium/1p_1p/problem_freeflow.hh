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
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution (Kovasznay 1948, \cite Kovasznay1948)
 */

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_NEW_HH


#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Kovasznay 1948, \cite Kovasznay1948)
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered. The set-up represents a wake behind a two-dimensional grid
 * and is chosen in a way such that an exact solution is available.
 */
template <class TypeTag>
class FreeFlowOnePTestProblem :  public NavierStokesProblem<TypeTag>
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

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    // static constexpr auto upwindSchemeOrder = getPropValue<TypeTag, Properties::UpwindSchemeOrder>();

public:
    FreeFlowOnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, "FreeFlow")
    , couplingManager_(couplingManager)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference", 0.0);

        // determine whether to simulate a vertical or horizontal flow configuration
        verticalFlow_ = problemName_.find("vertical") != std::string::npos;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        const auto& globalPos = scvf.center(); //avoid ambiguities at corners

        if (verticalFlow_)
        {
            if constexpr (ParentType::isMomentumProblem())
            {
                if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf))
                {
                    values.setCouplingNeumann(Indices::momentumYBalanceIdx);
                    values.setCouplingNeumann(Indices::momentumXBalanceIdx);
                }
                else
                    values.setAllDirichlet();
            }
            else
            {
                if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, scvf))
                    values.setAllCouplingNeumann();
                else
                    values.setAllNeumann();
            }
        }
        else // horizontal flow
        {
            if constexpr (ParentType::isMomentumProblem())
            {
                if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
                    values.setAllNeumann();
                else if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf))
                {
                    values.setCouplingNeumann(Indices::momentumYBalanceIdx);
                    values.setCouplingNeumann(Indices::momentumXBalanceIdx);
                }
                else
                    values.setAllDirichlet();
            }
            else
            {
                if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, scvf))
                    values.setAllCouplingNeumann();
                else
                    values.setAllNeumann();
            }
        }

        return values;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (verticalFlow_)
            {
                static const Scalar inletVelocity = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity");
                values[Indices::velocityYIdx] = inletVelocity * globalPos[0] * (this->gridGeometry().bBoxMax()[0] - globalPos[0]);
            }
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
        const auto& globalPos = scvf.ipGlobal();
        using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;

        if constexpr (ParentType::isMomentumProblem())
        {

            if (couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf))
            {
                values += couplingManager_->momentumCouplingCondition(CouplingManager::freeFlowMomentumIndex,
                                                                      CouplingManager::porousMediumIndex,
                                                                      fvGeometry, scvf, elemVolVars);

                values += FluxHelper::slipVelocityMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache);
            }
            else
            {
                const Scalar pressure = onLeftBoundary_(globalPos) ? deltaP_ : 0.0;
                values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, pressure, true /*zeroNormalVelocityGradient*/);
            }
        }
        else
        {
            if (couplingManager_->isCoupled(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex, scvf))
            {
                values = couplingManager_->massCouplingCondition(CouplingManager::freeFlowMassIndex, CouplingManager::porousMediumIndex,
                                                                 fvGeometry, scvf, elemVolVars);
            }
            else
            {
                using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
            }
        }

        return values;
    }

    /*!
     * \brief Returns true if the scvf lies on a porous slip boundary
     */
    bool onSlipBoundary(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        assert(scvf.isFrontal());
        return scvf.boundary() && couplingManager_->isCoupled(CouplingManager::freeFlowMomentumIndex, CouplingManager::porousMediumIndex, scvf);
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    { return couplingManager_->darcyPermeability(fvGeometry, scvf); }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    { return couplingManager_->problem(CouplingManager::porousMediumIndex).spatialParams().beaversJosephCoeffAtPos(scvf.ipGlobal()); }

    /*!
     * \brief Returns the pressure difference between inlet and outlet for the horizontal flow scenario
     */
    Scalar pressureDifference() const
    {
        assert(!verticalFlow_);
        return deltaP_;
    }

    /*!
     * \brief Returns true if a vertical flow scenario is considered
     */
    bool verticalFlow() const
    { return verticalFlow_; }

    // \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    std::string problemName_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar deltaP_;
    bool verticalFlow_;

    std::shared_ptr<CouplingManager> couplingManager_;
};
} // end namespace Dumux

#endif
