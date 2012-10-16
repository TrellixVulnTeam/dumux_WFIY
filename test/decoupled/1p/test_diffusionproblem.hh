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
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
#ifndef DUMUX_TEST_2P_PROBLEM_HH
#define DUMUX_TEST_2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/components/unit.hh>

#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/omethod/fvmpfaopressurepropertiessimple2p.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticpressureproperties2p.hh>
#include <dumux/decoupled/2p/diffusion/diffusionproblem2p.hh>
#include <dumux/decoupled/common/fv/fvvelocity.hh>

#include "test_diffusionspatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup IMPETtests
 */
template<class TypeTag>
class TestDiffusionProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
//// set the types for the 2PFA FV method
NEW_TYPE_TAG(FVVelocity2PTestProblem, INHERITS_FROM(FVPressureTwoP, TestDiffusionSpatialParams));
SET_TYPE_PROP(FVVelocity2PTestProblem, Problem, Dumux::TestDiffusionProblem<TypeTag>);
// Set the grid type
SET_PROP(FVVelocity2PTestProblem, Grid)
{//    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

SET_INT_PROP(FVVelocity2PTestProblem, Formulation,
        DecoupledTwoPCommonIndices::pGlobalSw);

// Set the wetting phase
SET_PROP(FVVelocity2PTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(FVVelocity2PTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(FVVelocity2PTestProblem, ProblemEnableGravity, false);


// set the types for the MPFA-O FV method
NEW_TYPE_TAG(FVMPFAOVelocity2PTestProblem, INHERITS_FROM(FVMPFAOPressurePropertiesTwoP, TestDiffusionSpatialParams));
//SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, LinearSolver, Dumux::ILUnBiCGSTABBackend<TypeTag>);
SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, LinearSolver, Dumux::SSORBiCGSTABBackend<TypeTag>);
SET_INT_PROP(FVMPFAOVelocity2PTestProblem, PreconditionerIterations, 2);
SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, Problem, Dumux::TestDiffusionProblem<TypeTag>);
// Set the grid type
SET_PROP(FVMPFAOVelocity2PTestProblem, Grid)
{//    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

SET_INT_PROP(FVMPFAOVelocity2PTestProblem, Formulation,
        DecoupledTwoPCommonIndices::pGlobalSw);

// Set the wetting phase
SET_PROP(FVMPFAOVelocity2PTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(FVMPFAOVelocity2PTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(FVMPFAOVelocity2PTestProblem, ProblemEnableGravity, false);

// set the types for the mimetic FD method
NEW_TYPE_TAG(MimeticPressure2PTestProblem, INHERITS_FROM(MimeticPressureTwoP, TestDiffusionSpatialParams));
SET_TYPE_PROP(MimeticPressure2PTestProblem, Problem, Dumux::TestDiffusionProblem<TypeTag>);
// Set the grid type
SET_PROP(MimeticPressure2PTestProblem, Grid)
{
//        typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

SET_INT_PROP(MimeticPressure2PTestProblem, Formulation,
        DecoupledTwoPCommonIndices::pGlobalSw);

// Set the wetting phase
SET_PROP(MimeticPressure2PTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(MimeticPressure2PTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(MimeticPressure2PTestProblem, ProblemEnableGravity, false);

}

/*!
 * \ingroup DecoupledProblems
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
template<class TypeTag>
class TestDiffusionProblem: public DiffusionProblem2P<TypeTag>
{
    typedef DiffusionProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) NonwettingPhase;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pGlobalIdx = Indices::pGlobalIdx,
        SwIdx = Indices::SwIdx,
        pressEqIdx = Indices::pressEqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename SolutionTypes::ScalarSolution ScalarSolution;

    TestDiffusionProblem(const GridView &gridView, const Scalar delta = 1.0) :
        ParentType(gridView), delta_(delta), velocity_(*this)
    {}

    //!for this specific problem: initialize the saturation and afterwards the model
    void init()
    {
        this->variables().initialize();
        this->spatialParams().initialize(delta_);
        for (int i = 0; i < this->gridView().size(0); i++)
        {
            this->variables().cellData(i).setSaturation(wPhaseIdx, 1.0);
        }
        this->model().initialize();
        velocity_.initialize();
    }

    /*!
    * \name Problem parameters
    */
    // \{

    bool shouldWriteRestartFile() const
    { return false; }


    void calculateFVVelocity()
    {
        velocity_.calculateVelocity();
//        velocity_.addOutputVtkFields(this->resultWriter());
    }

    //! \copydoc ParentType::addOutputVtkFields()
    void addOutputVtkFields()
    {
        ScalarSolution *exactPressure = this->resultWriter().allocateManagedBuffer(this->gridView().size(0));

        ElementIterator eIt = this->gridView().template begin<0>();
        ElementIterator eItEnd = this->gridView().template end<0>();
        for(;eIt != eItEnd; ++eIt)
        {
            (*exactPressure)[this->elementMapper().map(*eIt)][0] = exact(eIt->geometry().center());
        }

        this->resultWriter().attachCellData(*exactPressure, "exact pressure");

        this->spatialParams().addOutputVtkFields(this->resultWriter());

        return;
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

    void source(PrimaryVariables &values,const Element& element) const
    {
        values = 0;

        values[wPhaseIdx] = integratedSource_(element, 4);
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
    *
    * BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
    */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {

          bcTypes.setAllDirichlet();
    }

    //! set dirichlet condition  (saturation [-])
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
            values[pGlobalIdx] = exact(globalPos);
            values[SwIdx] = 1.0;
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    Scalar exact (const GlobalPosition& globalPos) const
    {
        Scalar pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
        {
        Dune::FieldVector<Scalar,dim> grad(0);
        Scalar pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
        }

private:

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

        return -(fx + fy) / WettingPhase::viscosity(temp, referencePress) * WettingPhase::density(temp, referencePress);
    }

    Scalar delta_;
    Dumux::FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity) > velocity_;
};
} //end namespace

#endif
