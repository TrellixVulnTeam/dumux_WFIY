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
 * \brief Test for the staggered grid Navier-Stokes model
 *        with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROBLEM_HH
#define DUMUX_SINCOS_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dune/geometry/quadraturerules.hh>


#include "../l2error.hh"

namespace Dumux {
template <class TypeTag>
class SincosTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct SincosTest {};
struct SincosTestMomentum { using InheritsFrom = std::tuple<SincosTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct SincosTestMass { using InheritsFrom = std::tuple<SincosTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SincosTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::SincosTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SincosTest> { using type = Dumux::SincosTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid.
 *
 * The 2D, incompressible Navier-Stokes equations for zero gravity and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * For the instationary case, the velocities and pressures are periodical in time. The Dirichlet boundary conditions are
 * consistent with the analytical solution and in the instationary case time-dependent.
 */
template <class TypeTag>
class SincosTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NumEqVector = typename ParentType::NumEqVector;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;


    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    SincosTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager), time_(0.0)
    {
        isStationary_ = getParam<bool>("Problem.IsStationary");
        density_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        dynamicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0) * density_;
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);

        if constexpr(ParentType::isMomentumProblem())
        {
            if (useNeumann_ && !isStationary_)
                DUNE_THROW(Dune::NotImplemented, "Neumann boundary conditions only implemented for stationary case.");
        }
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 298.0; }

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            const Scalar t = time_;
            const Scalar mu = dynamicViscosity_;
            const Scalar rho = density_;

            using std::cos;
            using std::sin;
            using Dune::power;

            if (isStationary_)
            {
                if (this->enableInertiaTerms())
                {
                    source[Indices::momentumXBalanceIdx] = (-2.0*mu*sin(y) + (1.0 - rho)*sin(x)) * cos(x);
                    source[Indices::momentumYBalanceIdx] =  (2.0*mu*sin(x) + (1.0 - rho)*sin(y)) * cos(y);
                }
                else
                {
                    source[Indices::momentumXBalanceIdx] = (-2.0*mu*sin(y) + sin(x)) * cos(x);
                    source[Indices::momentumYBalanceIdx] = ( 2.0*mu*sin(x) + sin(y)) * cos(y);
                }
            }
            else
            {
                if (this->enableInertiaTerms())
                {
                    source[Indices::momentumXBalanceIdx] = (-mu*(cos(2*t - y) - cos(2*t + y)) + 4*rho*power(sin(t), 4)*sin(x) - 4*rho*power(sin(t), 2)*sin(x)
                                                           + 4*rho*power(sin(t), 2)*sin(y) - 2*rho*sin(y) - 4*power(sin(t), 4)*sin(x) + 4*power(sin(t), 2)*sin(x))*cos(x);
                    source[Indices::momentumYBalanceIdx] =  (4*mu*sin(t)*sin(x)*cos(t) + 4*rho*power(sin(t), 4)*sin(y) - 4*rho*power(sin(t), 2)*sin(x) - 4*rho*power(sin(t), 2)*sin(y)
                                                             + 2*rho*sin(x) - 4*power(sin(t), 4)*sin(y) + 4.0*power(sin(t), 2)*sin(y))*cos(y);
                }
                else
                {

                }
            }
        }

        return source;
    }

    // \}
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

        if constexpr (ParentType::isMomentumProblem())
        {
            if (useNeumann_)
            {
                if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6
                || globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6
                || globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                {
                    values.setNeumann(Indices::velocityXIdx);
                    values.setNeumann(Indices::velocityYIdx);
                }
                else
                {
                    values.setDirichlet(Indices::velocityXIdx);
                    values.setDirichlet(Indices::velocityYIdx);
                }
            }
            else
            {
                // set Dirichlet values for the velocity everywhere
                values.setAllDirichlet();
            }
        }
        else
            values.setAllNeumann();

        return values;
    }

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
        // set fixed pressure in one cell
        std::bitset<PrimaryVariables::dimension> values;

        if (!useNeumann_ && scv.dofIndex() == 0)
        {
            assert(false);
            values.set(Indices::pressureIdx);
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.center())[Indices::pressureIdx]); }


   /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos, time_);
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

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto flux = [&](Scalar x, Scalar y)
            {
                Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
                using std::sin;
                using std::cos;
                const Scalar mu = dynamicViscosity_;
                const Scalar rho = density_;

                if (this->enableInertiaTerms())
                {
                    momentumFlux[0][0] = -2*mu*sin(x)*sin(y) + rho*sin(y)*sin(y)*cos(x)*cos(x) - 0.25*cos(2*x) - 0.25*cos(2*y);
                    momentumFlux[0][1] = -rho*(cos(2*x - 2*y) - cos(2*x + 2*y))/8.0;
                    momentumFlux[1][0] = momentumFlux[0][1];
                    momentumFlux[1][1] = 2*mu*sin(x)*sin(y) + rho*sin(x)*sin(x)*cos(y)*cos(y) - 0.25*cos(2*x) - 0.25*cos(2*y);
                }
                else
                {
                    momentumFlux[0][0] = -2*mu*sin(x)*sin(y) - 0.25*cos(2*x) - 0.25*cos(2*y);
                    momentumFlux[0][1] = 0;
                    momentumFlux[1][0] = 0;
                    momentumFlux[1][1] = 2*mu*sin(x)*sin(y) - 0.25*cos(2*x) - 0.25*cos(2*y);
                }
                return momentumFlux;
            };

            flux(scvf.ipGlobal()[0], scvf.ipGlobal()[1]).mv(scvf.unitOuterNormal(), values);
        }
        else
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = insideDensity * (this->faceVelocity(element, fvGeometry, scvf) *  scvf.unitOuterNormal());
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar t = time;
        PrimaryVariables values;

        using std::sin;
        using std::cos;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = -1.0 * cos(x) * sin(y);
            values[Indices::velocityYIdx] = sin(x) * cos(y);
        }
        else
            values[Indices::pressureIdx] = -0.25 * (cos(2.0 * x) + cos(2.0 * y));

        if (!isStationary_)
        {
            if constexpr (ParentType::isMomentumProblem())
            {
                values[Indices::velocityXIdx] *= sin(2.0 * t);
                values[Indices::velocityYIdx] *= sin(2.0 * t);
            }
            else
                values[Indices::pressureIdx] *= sin(2.0 * t) * sin(2.0 * t);
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos, time_);
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
        if (isStationary_)
            return PrimaryVariables(0.0);
        else
            return analyticalSolution(globalPos, 0.0);
    }

    /*!
     * \brief Updates the time
     */
    void updateTime(const Scalar time)
    { time_ = time; }

private:
    Scalar dynamicViscosity_;
    Scalar density_;
    Scalar time_;
    bool isStationary_;
    bool useNeumann_ = false;
};

} // end namespace Dumux

#endif
