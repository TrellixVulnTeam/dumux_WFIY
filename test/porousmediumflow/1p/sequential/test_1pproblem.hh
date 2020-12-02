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
 *
 * \brief test problem for the sequential one-phase model.
 */
#ifndef DUMUX_TEST_1P_PROBLEM_HH
#define DUMUX_TEST_1P_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/porousmediumflow/1p/sequential/diffusion/problem.hh>
namespace Dumux
{

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential one-phase model.
 */
template<class TypeTag>
class TestProblemOneP: public DiffusionProblem1P<TypeTag >
{
    using ParentType = DiffusionProblem1P<TypeTag>;
    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using LocalPosition = Dune::FieldVector<Scalar, dim>;


public:
    TestProblemOneP(TimeManager& timeManager, Grid& grid) :
        ParentType(grid), velocity_(*this)
    {
        delta_ = getParam<Scalar>("Problem.Delta", 1e-3);

        this->spatialParams().initialize(delta_);
    }

    /*!
    * \name Problem parameters
    */
    // \{

    /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
    std::string name() const
    {
        return "test_1p";
    }

    bool shouldWriteRestartFile() const
    { return false; }

    void addOutputVtkFields()
    {
        velocity_.calculateVelocity();
        velocity_.addOutputVtkFields(this->resultWriter());
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    //!source term [kg/(m^3 s)]
    void source(PrimaryVariables &values, const Element& element) const
    {
        values = 0;
        values = integratedSource_(element, 4);
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
    * BC can be dirichlet (pressure) or neumann (flux).
    */
    void boundaryTypes(BoundaryTypes &bcType,
            const Intersection& intersection) const
    {
        bcType.setAllDirichlet();
    }

    //! return dirichlet condition  (pressure, [Pa])
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = exact(globalPos);
    }


    //! return neumann condition  (flux, [kg/(m^2 s)])
    void neumann(PrimaryVariables &values, const Intersection& intersection) const
    {
        values = 0;
    }

private:
    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);
        using std::sin;
        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar,dim> grad(0);
        using std::sin;
        using std::cos;
        using std::atan;
        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
    }

    Scalar integratedSource_(const Element& element, int integrationPoints) const
    {
        Scalar source = 0.;
        LocalPosition localPos(0.0);
        GlobalPosition globalPos(0.0);
        Scalar halfInterval = 1.0/double(integrationPoints)/2.;
        for (int i = 1; i <= integrationPoints; i++)
        {
            for (int j = 1; j <= integrationPoints; j++)
            {
                localPos[0] = double(i)/double(integrationPoints) - halfInterval;
                localPos[1] = double(j)/double(integrationPoints) - halfInterval;
                globalPos = element.geometry().global(localPos);
                source += 1./(integrationPoints*integrationPoints) * evaluateSource_(globalPos);
            }
        }

        return source;
    }

    Scalar evaluateSource_(const GlobalPosition& globalPos) const
    {
        Scalar temp = temperatureAtPos(globalPos);
        Scalar referencePress = referencePressureAtPos(globalPos);

        using std::sin;
        using std::cos;
        using std::atan;

        Scalar pi = 4.0 * atan(1.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        Scalar dpdx = pi * cos(pi * x) * sin(pi * y);
        Scalar dpdy = pi * sin(pi * x) * cos(pi * y);
        Scalar dppdxx = -pi * pi * sin(pi * x) * sin(pi * y);
        Scalar dppdxy = pi * pi * cos(pi * x) * cos(pi * y);
        Scalar dppdyx = dppdxy;
        Scalar dppdyy = dppdxx;
        Scalar kxx = (delta_* x*x + y*y)/(x*x + y*y);
        Scalar kxy = -(1.0 - delta_) * x * y / (x*x + y*y);
        Scalar kyy = (x*x + delta_*y*y)/(x*x + y*y);
        Scalar dkxxdx = 2 * x * y*y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkyydy = 2 * x*x * y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkxydx = (1.0 - delta_) * y * (x*x - y*y) /((x*x + y*y) * (x*x + y*y));
        Scalar dkxydy = (1.0 - delta_) * x * (y*y - x*x) /((x*x + y*y) * (x*x + y*y));

        Scalar fx = dkxxdx * dpdx + kxx * dppdxx + dkxydx * dpdy + kxy * dppdyx;
        Scalar fy = dkxydy * dpdx + kxy * dppdxy + dkyydy * dpdy + kyy * dppdyy;

        return -(fx + fy) / FluidSystem::viscosity(temp, referencePress) * FluidSystem::density(temp, referencePress);
    }

    double delta_;
    FVVelocity<TypeTag, GetPropType<TypeTag, Properties::Velocity> > velocity_;
};
} //end namespace

#endif
