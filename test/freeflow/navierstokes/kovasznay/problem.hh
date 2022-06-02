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

#ifndef DUMUX_KOVASZNAY_TEST_PROBLEM_HH
#define DUMUX_KOVASZNAY_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

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
class KovasznayTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    KovasznayTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        //std::cout<< "upwindSchemeOrder is: " << GridGeometry::upwindStencilOrder() << "\n";
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);
    }


    KovasznayTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //std::cout<< "upwindSchemeOrder is: " << GridGeometry::upwindStencilOrder() << "\n";
        rho_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        Scalar reynoldsNumber = 1.0 / kinematicViscosity_;
        lambda_ = 0.5 * reynoldsNumber
                        - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values everywhere
        values.setAllDirichlet();

        return values;
    }

   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition & globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos);
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
        return values;
    }

   /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this is a stationary test
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        PrimaryVariables  values;
        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
            values[Indices::velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);
        }
        else
            values[Indices::pressureIdx] = pressureAtPos(globalPos);

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return analyticalSolution(globalPos, 0);
    }

    //! TODO should these be spatial params?
    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return rho_ * 0.5 * (1.0 - std::exp(2.0 * lambda_ * globalPos[0])); }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    { return rho_; }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return kinematicViscosity_; }

private:
    static constexpr Scalar eps_=1e-6;

    Scalar rho_;
    Scalar kinematicViscosity_;
    Scalar lambda_;
};
} // end namespace Dumux

#endif
