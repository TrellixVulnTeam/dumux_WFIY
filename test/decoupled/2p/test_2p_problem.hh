// $Id: test_2p_problem.hh 3783 2010-06-24 11:33:53Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_TEST_2P_PROBLEM_HH
#define DUMUX_TEST_2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/oil.hh>

#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/fvmpfaovelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include <dumux/decoupled/2p/transport/fv/capillarydiffusion.hh>
#include <dumux/decoupled/2p/transport/fv/gravitypart.hh>

#include "test_2p_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class Test2PProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TwoPTestProblem, INHERITS_FROM(DecoupledTwoP, MPFAProperties, Transport));

// Set the grid type
SET_PROP(TwoPTestProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

// Set the problem property
SET_PROP(TwoPTestProblem, Problem)
{
public:
    typedef Dumux::Test2PProblem<TTAG(TwoPTestProblem)> type;
};

// Set the model properties
SET_PROP(TwoPTestProblem, SaturationModel)
{
    typedef Dumux::FVSaturation2P<TTAG(TwoPTestProblem)> type;
};
SET_TYPE_PROP(TwoPTestProblem, DiffusivePart, Dumux::CapillaryDiffusion<TypeTag>);
SET_TYPE_PROP(TwoPTestProblem, ConvectivePart, Dumux::GravityPart<TypeTag>);

SET_PROP(TwoPTestProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(TwoPTestProblem)> type;
//    typedef Dumux::FVMPFAOVelocity2P<TTAG(TwoPTestProblem)> type;
};

//SET_INT_PROP(TwoPTestProblem, VelocityFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

//SET_INT_PROP(TwoPTestProblem, PressureFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureGlobal);

// Set the wetting phase
SET_PROP(TwoPTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(TwoPTestProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::Test2PSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(TwoPTestProblem, EnableGravity, false);

SET_SCALAR_PROP(TwoPTestProblem, CFLFactor, 0.95);
}

/*!
 * \ingroup DecoupledProblems
 */
template<class TypeTag = TTAG(TwoPTestProblem)>
class Test2PProblem: public IMPESProblem2P<TypeTag, Test2PProblem<TypeTag> >
{
typedef Test2PProblem<TypeTag> ThisType;
typedef IMPESProblem2P<TypeTag, ThisType> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
};

typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
Test2PProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
{
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
const char *name() const
{
    return "test2p";
}

bool shouldWriteRestartFile() const
{
    return false;
}

/*!
 * \brief Returns the temperature within the domain.
 *
 * This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
{
    return 273.15 + 10; // -> 10°C
}

// \}


Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
{
    return 1e5; // -> 10°C
}

std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
{
    return std::vector<Scalar>(2, 0.0);
}

typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if ((globalPos[0] < eps_))
    return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}

BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] > (upperRight_[0] - eps_) || globalPos[0] < eps_)
//    if (globalPos[0] < eps_)
    return Dumux::BoundaryConditions::dirichlet;
    else
    return Dumux::BoundaryConditions::neumann;
}

Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
        {
            const Element& element = *(intersection.inside());

            Scalar pRef = referencePressure(globalPos, element);
            Scalar temp = temperature(globalPos, element);
            Scalar sat = this->variables().satElement(element);

            FluidState fluidState;
            fluidState.update(sat, pRef, pRef, temp);
            return (2e5 + (upperRight_[dim-1] - globalPos[dim-1]) * FluidSystem::phaseDensity(wPhaseIdx, temp, pRef, fluidState) * this->gravity().two_norm());
        }
        else
        return 2e5;
    }
    // all other boundaries
    return 2e5;
}

Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    return 0.8;
    // all other boundaries
    return 0.2;
}

std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const
{
    std::vector<Scalar> neumannFlux(2, 0.0);
    if (globalPos[0] > upperRight_[0] - eps_)
    {
        neumannFlux[nPhaseIdx] = 3e-4;
    }
    return neumannFlux;
}

Scalar neumannSat(const GlobalPosition& globalPos, const Intersection& intersection, Scalar factor) const
{
    if (globalPos[0] > upperRight_[0] - eps_)
    return factor;
    return 0;
}

Scalar initSat(const GlobalPosition& globalPos, const Element& element) const
{
    return 0.2;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
