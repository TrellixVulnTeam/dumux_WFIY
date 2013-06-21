//$Id: el2pproperties.hh 5404 2011-03-11 10:30:54Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines the properties required for the two phase linear-elastic model.
 *
 * This class inherits from the properties of the two-phase model and
 * from the properties of the simple linear-elastic model
 */

#ifndef DUMUX_ELASTIC2P_PROPERTY_DEFAULTS_HH
#define DUMUX_ELASTIC2P_PROPERTY_DEFAULTS_HH

#include "el2pproperties.hh"

#include "el2pmodel.hh"
#include "el2pindices.hh"
#include "el2plocalresidual.hh"
#include "el2plocaljacobian.hh"
#include "el2pfluxvariables.hh"
#include "el2pelementvolumevariables.hh"
#include "el2pvolumevariables.hh"
#include "el2plocaloperator.hh"
#include "el2passembler.hh"
#include "el2pnewtoncontroller.hh"
#include "el2pindices.hh"
#include <dumux/implicit/box/boxpropertydefaults.hh>
#include <dumux/implicit/2p/2ppropertydefaults.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/common/function.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>

namespace Dumux
{

//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////

namespace Properties
{
SET_INT_PROP(BoxElasticTwoP, NumEq, 5); //!< set the number of equations to 5
SET_INT_PROP(BoxElasticTwoP, NumPhases, 2); //!< The number of fluid phases in the elastic 2p model is 2

//! Use the elastic local jacobian operator for the two-phase linear-elastic model
SET_TYPE_PROP(BoxElasticTwoP,
              LocalResidual,
              ElTwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxElasticTwoP, Model, ElTwoPModel<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(BoxElasticTwoP, ElementVolumeVariables, Dumux::ElTwoPElementVolumeVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxElasticTwoP, VolumeVariables, ElTwoPVolumeVariables<TypeTag>);

//! Set the default formulation to pWsN
SET_INT_PROP(BoxElasticTwoP,
             Formulation,
             0);

//! The indices required by the two-phase linear-elastic model

SET_PROP(BoxElasticTwoP, Indices)
{
    typedef ElTwoPIndices<TypeTag> type;
};

//! The FluxVariables required by the two-phase linear-elastic model
SET_TYPE_PROP(BoxElasticTwoP, FluxVariables, ElTwoPFluxVariables<TypeTag>);

//! the default upwind factor. Default 1.0, i.e. fully upwind...
SET_SCALAR_PROP(BoxElasticTwoP, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(BoxElasticTwoP, ImplicitMobilityUpwindWeight, 1.0);

//! enable gravity by default
SET_BOOL_PROP(BoxElasticTwoP, ProblemEnableGravity, true);


//! Enable evaluation of shape function gradients at the sub-control volume center by default
// Used for the computation of the pressure gradients
SET_BOOL_PROP(BoxElasticTwoP, EvalGradientsAtSCVCenter, true);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxElasticTwoP, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};


// use the SuperLU linear solver by default
SET_TYPE_PROP(BoxElasticTwoP, LinearSolver, Dumux::SuperLUBackend<TypeTag> );

// set the grid operator
SET_PROP(BoxElasticTwoP, GridOperator)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsTrafo) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, LocalOperator) LocalOperator;


    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:

    typedef Dune::PDELab::GridOperator<GridFunctionSpace,
            GridFunctionSpace,
            LocalOperator,
            Dune::PDELab::ISTLBCRSMatrixBackend<1, 1>,
            Scalar, Scalar, Scalar,
            ConstraintsTrafo,
            ConstraintsTrafo,
            true
            > type;
};

SET_PROP(BoxElasticTwoP, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
public:
    typedef typename GridOperator::Traits::Jacobian type;
};

SET_PROP(BoxElasticTwoP, SolutionVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridOperator) GridOperator;
public:
    typedef typename GridOperator::Traits::Domain type;
};

SET_PROP(BoxElasticTwoP, PressureGridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureFEM)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
         dim = GridView::dimension};
public:
    typedef Dune::PDELab::NoConstraints Constraints;

    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints, Dune::PDELab::ISTLVectorBackend<1> >
        ScalarGridFunctionSpace;

    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq-dim, Dune::PDELab::GridFunctionSpaceBlockwiseMapper>
        type;

    typedef typename type::template ConstraintsContainer<Scalar>::Type
        ConstraintsTrafo;
};

SET_PROP(BoxElasticTwoP, DisplacementGridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DisplacementFEM)) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::NoConstraints Constraints;

    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints, Dune::PDELab::ISTLVectorBackend<1> >
        ScalarGridFunctionSpace;

    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, dim, Dune::PDELab::GridFunctionSpaceBlockwiseMapper>
        type;

    typedef typename type::template ConstraintsContainer<Scalar>::Type
        ConstraintsTrafo;
};

SET_PROP(BoxElasticTwoP, GridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureGridFunctionSpace)) PressureGFS;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DisplacementGridFunctionSpace)) DisplacementGFS;
    typedef Dune::PDELab::GridFunctionSpaceLexicographicMapper GFSMapper;
public:
    typedef Dune::PDELab::NoConstraints Constraints;

    typedef void ScalarGridFunctionSpace;

    typedef Dune::PDELab::CompositeGridFunctionSpace<GFSMapper, PressureGFS, DisplacementGFS> type;

    typedef typename type::template ConstraintsContainer<Scalar>::Type
        ConstraintsTrafo;
};

SET_PROP(BoxElasticTwoP, ConstraintsTrafo)
{private:
  typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type type;
};

// set the grid function space for the sub-models
SET_TYPE_PROP(BoxElasticTwoP, Constraints, Dune::PDELab::NoConstraints);

SET_TYPE_PROP(BoxElasticTwoP, JacobianAssembler, Dumux::PDELab::El2PAssembler<TypeTag>);

SET_PROP(BoxElasticTwoP, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxElasticTwoP, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(BoxElasticTwoP, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(BoxElasticTwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxElasticTwoP, ImplicitEnableJacobianRecycling, false);
// enable partial reassembling by default
SET_BOOL_PROP(BoxElasticTwoP, ImplicitEnablePartialReassemble, false);

SET_TYPE_PROP(BoxElasticTwoP, NewtonController, ElTwoPNewtonController<TypeTag>);

SET_PROP(BoxElasticTwoP, LocalOperator)
{
    typedef Dumux::PDELab::El2PLocalOperator<TypeTag> type;
};

//! use the local FEM space associated with cubes by default
SET_PROP(BoxElasticTwoP, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};
/*!
 * \brief A vector of primary variables.
 */
SET_PROP(BoxElasticTwoP, PrimaryVariables)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
public:
    typedef Dune::FieldVector<Scalar, numEq> type;
};

//! The local jacobian operator
SET_TYPE_PROP(BoxElasticTwoP, LocalJacobian, Dumux::ElTwoPLocalJacobian<TypeTag>);

//! set number of equations of the mathematical model as default
SET_INT_PROP(BoxElasticTwoP, LinearSolverBlockSize, 1);

// write the stress and displacement output according to rock mechanics sign convention (compressive stresses > 0)
SET_BOOL_PROP(BoxElasticTwoP, VtkRockMechanicsSignConvention, true);

// \}
}
}

#endif

