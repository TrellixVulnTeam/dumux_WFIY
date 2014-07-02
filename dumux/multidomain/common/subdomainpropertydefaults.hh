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
 * \brief Specify default properties required in the subdomains of dune-multidomain
 */
#ifndef DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH
#define DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH

#include <dune/grid/multidomaingrid.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/finiteelementmap/q1fem.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>

#include "subdomainproperties.hh"
#include "multidomainproperties.hh"
#include "multidomainlocaloperator.hh"
#include <dumux/multidomain/couplinglocalresiduals/boxcouplinglocalresidual.hh>

namespace Dumux
{

namespace Properties
{

// Specifies the grid type for the subdomains
SET_PROP(SubDomain, Grid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomain;
    typedef typename GET_PROP_TYPE(MultiDomain, Grid) HostGrid;
    typedef typename Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,4> MDGridTraits;
    typedef typename Dune::MultiDomainGrid<HostGrid, MDGridTraits> Grid;
public:
    typedef typename Grid::SubDomainGrid type;
};

// set the default BaseLocalResidual to BoxCouplingLocalResidual
SET_TYPE_PROP(SubDomain, BaseLocalResidual, BoxCouplingLocalResidual<TypeTag>);

// set the local operator used for submodels
SET_TYPE_PROP(SubDomain, LocalOperator,
              Dumux::PDELab::MultiDomainLocalOperator<TypeTag>);

// use the time manager for the coupled problem in the sub problems
SET_PROP(SubDomain, TimeManager)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomainTypeTag;
public:
    typedef typename GET_PROP_TYPE(MultiDomainTypeTag, TimeManager) type;
};

// set the grid functions space for the sub-models
SET_PROP(SubDomain, ScalarGridFunctionSpace)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalFEMSpace) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
 public:
    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints,
        Dune::PDELab::ISTLVectorBackend<1> > type;
};

// set the grid functions space for the sub-models
SET_PROP(SubDomain, GridFunctionSpace)
{private:
    typedef typename GET_PROP_TYPE(TypeTag, ScalarGridFunctionSpace) ScalarGridFunctionSpace;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
public:
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq, Dune::PDELab::GridFunctionSpaceBlockwiseMapper> type;
};

// set the grid function space for the sub-models
SET_TYPE_PROP(SubDomain, Constraints, Dune::PDELab::NoConstraints);

// set the grid functions space for the sub-models
SET_PROP(SubDomain, ConstraintsTrafo)
{private:
  typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
public:
    typedef typename GridFunctionSpace::template ConstraintsContainer<Scalar>::Type type;
};

// set the grid operator space used for submodels
// DEPRECATED: use GridOperator instead
SET_PROP(SubDomain, GridOperatorSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsTrafo) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, LocalOperator) LocalOperator;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:
    typedef Dune::PDELab::GridOperatorSpace<GridFunctionSpace,
        GridFunctionSpace,
        LocalOperator,
        ConstraintsTrafo,
        ConstraintsTrafo,
        Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
        true
        > type;
};

// use the local FEM space associated with cubes by default
SET_PROP(SubDomain, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type;
};

SET_PROP(SubDomain, ParameterTree)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomainTypeTag;
    typedef typename GET_PROP(MultiDomainTypeTag, ParameterTree) ParameterTree;
public:
    typedef typename ParameterTree::type type;

    static type &tree()
    { return ParameterTree::tree(); }

    static type &compileTimeParams()
    { return ParameterTree::compileTimeParams(); }

    static type &runTimeParams()
    { return ParameterTree::runTimeParams(); }

    static type &deprecatedRunTimeParams()
    { return ParameterTree::deprecatedRunTimeParams(); }

    static type &unusedNewRunTimeParams()
    { return ParameterTree::unusedNewRunTimeParams(); }
};

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH
