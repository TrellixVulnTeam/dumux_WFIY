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
 * \brief The properties of the test for the staggered grid compositional Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_ANALYTICAL_NC_TEST_PROBLEM_HH
#define DUMUX_ANALYTICAL_NC_TEST_PROBLEM_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief  Test problem for the staggered grid.
 *
 * The 2D, incompressible compositional Navier-Stokes equations for zero gravity and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * The Dirichlet boundary conditions are consistent with the analytical solution.
 */
template <class TypeTag>
class AnalyticalNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr auto compIdx = 1;

public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    AnalyticalNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity");
        mu_ = rho_*nu; // dynamic viscosity
        diffusionCoeff_ = getParam<Scalar>("Component.DiffusionCoefficient");

        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        constrainCell_ = getParam<bool>("Problem.ConstrainCell",true);
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        Sources source(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        if constexpr (ParentType::isMomentumProblem())
        {
            source[Indices::momentumXBalanceIdx] = -2.0*mu_*dxxU_(x,y) - mu_*dyyU_(x,y) - mu_*dxyV_(x,y) + dxP_(x,y);
            source[Indices::momentumYBalanceIdx] = -2.0*mu_*dyyV_(x,y) - mu_*dxyU_(x,y) - mu_*dxxV_(x,y) + dyP_(x,y);

            if (this->enableInertiaTerms())
            {
                source[Indices::momentumXBalanceIdx] += rho_*dxUU_(x,y) + rho_*dyUV_(x,y);
                source[Indices::momentumYBalanceIdx] += rho_*dxUV_(x,y) + rho_*dyVV_(x,y);
            }
        }
        else
        {
            source[Indices::conti0EqIdx] = 0.0;
            source[Indices::conti0EqIdx + compIdx] = rho_*(dxUC_(x,y) + dyVC_(x,y)) - diffusionCoeff_*rho_*(dxxC_(x,y)+dyyC_(x,y));
        }

        return source;
    }

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
            if (useNeumann_)
            {
                if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6
                || globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6
                || globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                {
                    values.setNeumann(Indices::pressureIdx);
                    values.setNeumann(Indices::conti0EqIdx + compIdx);
                }
                else
                {
                    values.setDirichlet(Indices::pressureIdx);
                    values.setDirichlet(Indices::conti0EqIdx + compIdx);
                }
            }
            else
            {
                values.setDirichlet(Indices::pressureIdx);
                values.setDirichlet(Indices::conti0EqIdx + compIdx);
            }
        }

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
    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        // set fixed pressure in one cell
        std::bitset<DirichletValues::dimension> values;

        if (constrainCell_)
        {
            if (scv.dofIndex() == 15)
            {
                values.set(Indices::pressureIdx);
                values.set(Indices::conti0EqIdx + compIdx);
            }
        }
        else
        {
            GlobalPosition globalPos = scv.center();
            static const std::array<int,2> cells = getParam<std::array<int,2>>("Grid.Cells");
            GlobalPosition cellsize = {((this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]) / cells[0]),
                                    ((this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]) / cells[1])};
            if ((globalPos[1] < this->gridGeometry().bBoxMin()[1] + cellsize[1]) ||
                (globalPos[1] > this->gridGeometry().bBoxMax()[1] - cellsize[1]) ||
                (globalPos[0] > this->gridGeometry().bBoxMax()[0] - cellsize[0]))
            {
                values.set(Indices::pressureIdx);
                values.set(Indices::conti0EqIdx + compIdx);
            }
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return analyticalSolution(scv.center()); }

    /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

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
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto flux = [&](Scalar x, Scalar y)
            {
                Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);

                momentumFlux[0][0] = -2*mu_*dxU_(x,y) + p_(x,y);
                momentumFlux[0][1] = -mu_*(dyU_(x,y) + dxV_(x,y));
                momentumFlux[1][0] = -mu_*(dyU_(x,y) + dxV_(x,y));
                momentumFlux[1][1] = -2*mu_*dyV_(x,y) + p_(x,y);

                if (this->enableInertiaTerms())
                {
                    momentumFlux[0][0] += rho_*u_(x,y)*u_(x,y);
                    momentumFlux[0][1] += rho_*u_(x,y)*v_(x,y);
                    momentumFlux[1][0] += rho_*v_(x,y)*u_(x,y);
                    momentumFlux[1][1] += rho_*v_(x,y)*v_(x,y);
                }

                return momentumFlux;
            };

            flux(scvf.ipGlobal()[0], scvf.ipGlobal()[1]).mv(scvf.unitOuterNormal(), values);
        }
        else
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = insideDensity * (this->faceVelocity(element, fvGeometry, scvf) * scvf.unitOuterNormal());
            const auto flux = [&](Scalar x, Scalar y)
            {
                Dune::FieldVector<Scalar, dimWorld> componentFlux(0.0);
                componentFlux[0] = rho_*u_(x,y)*c_(x,y) - diffusionCoeff_* rho_ * dxC_(x,y);
                componentFlux[1] = rho_*v_(x,y)*c_(x,y) - diffusionCoeff_* rho_ * dyC_(x,y);
                return componentFlux;
            };
            values[Indices::conti0EqIdx + 1] = flux(scvf.ipGlobal()[0], scvf.ipGlobal()[1]) * scvf.unitOuterNormal();
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos,
                                       const Scalar& ) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        DirichletValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = u_(x,y);
            values[Indices::velocityYIdx] = v_(x,y);
        }
        else
        {
            values[Indices::pressureIdx] = p_(x,y);
            values[Indices::conti0EqIdx + compIdx] = c_(x,y);
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos, 0.0); }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    { return InitialValues(0.0); }

private:
    Scalar f_() const
    { return 1.0; }

    Scalar df_() const
    { return 0.0; }

    Scalar f1_(Scalar x) const
    { using std::cos; return -0.25 * cos(2.0 * x); }

    Scalar df1_(Scalar x) const
    { using std::sin; return 0.5 * sin(2.0 * x); }

    Scalar f2_(Scalar x) const
    { using std::cos; return -cos(x); }

    Scalar df2_(Scalar x) const
    { using std::sin; return sin(x); }

    Scalar ddf2_(Scalar x) const
    { using std::cos; return cos(x); }

    Scalar dddf2_(Scalar x) const
    { using std::sin; return -sin(x); }

    Scalar f3_(Scalar x) const
    { using std::cos; return 0.25*(cos(x)+1.0); }

    Scalar df3_(Scalar x) const
    { using std::sin; return -0.25*sin(x); }

    Scalar ddf3_(Scalar x) const
    { using std::cos; return -0.25*cos(x); }

    Scalar p_(Scalar x, Scalar y) const
    { return (f1_(x) + f1_(y)); }

    Scalar dxP_ (Scalar x, Scalar y) const
    { return df1_(x); }

    Scalar dyP_ (Scalar x, Scalar y) const
    { return df1_(y); }

    Scalar u_(Scalar x, Scalar y) const
    { return f2_(x)*df2_(y); }

    Scalar dxU_ (Scalar x, Scalar y) const
    { return df2_(x)*df2_(y); }

    Scalar dyU_ (Scalar x, Scalar y) const
    { return f2_(x)*ddf2_(y); }

    Scalar dxxU_ (Scalar x, Scalar y) const
    { return ddf2_(x)*df2_(y); }

    Scalar dxyU_ (Scalar x, Scalar y) const
    { return df2_(x)*ddf2_(y); }

    Scalar dyyU_ (Scalar x, Scalar y) const
    { return f2_(x)*dddf2_(y); }

    Scalar v_(Scalar x, Scalar y) const
    { return -f2_(y)*df2_(x); }

    Scalar dxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*ddf2_(x); }

    Scalar dyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*df2_(x); }

    Scalar dyyV_ (Scalar x, Scalar y) const
    { return -ddf2_(y)*df2_(x); }

    Scalar dxyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*ddf2_(x); }

    Scalar dxxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*dddf2_(x); }

    Scalar dxUU_ (Scalar x, Scalar y) const
    { return 2.*u_(x,y)*dxU_(x,y); }

    Scalar dyVV_ (Scalar x, Scalar y) const
    { return 2.*v_(x,y)*dyV_(x,y); }

    Scalar dxUV_ (Scalar x, Scalar y) const
    { return (v_(x,y) * dxU_(x,y)) + (u_(x,y) * dxV_(x,y)); }

    Scalar dyUV_ (Scalar x, Scalar y) const
    { return (v_(x,y) * dyU_(x,y)) + (u_(x,y) * dyV_(x,y)); }

    Scalar c_(Scalar x, Scalar y) const
    { return (f3_(x) + f3_(y)) + 1; }

    Scalar dxC_(Scalar x, Scalar) const
    { return df3_(x); }

    Scalar dyC_(Scalar, Scalar y) const
    { return df3_(y); }

    Scalar dxxC_(Scalar x, Scalar) const
    { return ddf3_(x); }

    Scalar dyyC_(Scalar, Scalar y) const
    { return ddf3_(y); }

    Scalar dxUC_(Scalar x, Scalar y) const
    { return (dxU_(x,y) * c_(x,y)) + (dxC_(x,y) * u_(x,y)); }

    Scalar dyVC_(Scalar x, Scalar y) const
    { return (dyV_(x,y) * c_(x,y)) + (dyC_(x,y) * v_(x,y)); }


    Scalar rho_;
    Scalar mu_;
    Scalar diffusionCoeff_;

    bool useNeumann_;
    bool constrainCell_;
};

} // end namespace Dumux

#endif
