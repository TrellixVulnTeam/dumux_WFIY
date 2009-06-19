// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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

#ifndef DUNE_MIMETICOPERATOR_HH
#define DUNE_MIMETICOPERATOR_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/functions/CRfunction.hh"
#include"dumux/shapefunctions/CRshapefunctions.hh"
//#include"localstiffnessext.hh"
#include"CRoperator.hh"

namespace Dune
{
/*! @brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - Grid    A grid type
  - Scalar   The field type used in the elements of the stiffness matrix
  - numEq    number of degrees of freedom per node (system size)
*/
template<class Grid, class Scalar, class GridView, class Communication, int numEq>
class MimeticOperatorAssembler : public CROperatorAssembler<Grid, Scalar, GridView, Communication, numEq>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum {dim=GridView::dimension};
    typedef P0Function<GridView,Scalar,2*dim> VType;
    typedef BlockVector< FieldVector<Scalar,1> > PType;
    typedef typename GridView::template Codim<0>::Iterator Iterator;
    typedef MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;

public:

    MimeticOperatorAssembler (const Grid& grid, const GridView& gridView, Communication lcomm)
    : CROperatorAssembler<Grid, Scalar, GridView, Communication, numEq>(grid, gridView, lcomm), elementMapper(gridView)
    {}

    template<class LocalStiffness>
    void calculatePressure (LocalStiffness& loc, CRFunction<Grid,Scalar,GridView,Communication,numEq>& u,
                            VType& velocity, PType& pressure)
    {
        // run over all level elements
        Iterator eendit = this->gridview.template end<0>();
        for (Iterator it = this->gridview.template begin<0>(); it!=eendit; ++it)
        {
            // get access to shape functions for CR elements
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::CRShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
                sfs = Dune::CRShapeFunctions<Scalar,Scalar,dim>::general(gt,1);

            int elemId = elementMapper.map(*it);

            // get local to global id map and pressure traces
            Dune::FieldVector<Scalar,2*dim> pressTrace(0);
            for (int k = 0; k < sfs.size(); k++)
            {
                pressTrace[k] = (*u)[this->facemapper.map(*it, k, 1)];
            }

            // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
            // The matrix W developed here corresponds to one element-associated
            // block of the matrix B^{-1} there.
            Dune::FieldVector<Scalar,2*dim> faceVol(0);
            Dune::FieldMatrix<Scalar,2*dim,2*dim> W(0);
            Dune::FieldVector<Scalar,2*dim> c(0);
            Dune::FieldMatrix<Scalar,2*dim,2*dim> Pi(0);
            Dune::FieldVector<Scalar,2*dim> F(0);
            Scalar dinv = 0;
            Scalar qmean = 0;
            loc.assembleElementMatrices(*it, faceVol, W, c, Pi, dinv, F, qmean);

            pressure[elemId] = dinv*(qmean + (F*pressTrace));

            Dune::FieldVector<Scalar,2*dim> v(0);
            for (int i = 0; i < 2*dim; i++)
                for (int j = 0; j < 2*dim; j++)
                    v[i] += W[i][j]*faceVol[j]*(pressure[elemId] - pressTrace[j]);

            (*velocity)[elemId] = v;
        }
    }

private:
    ElementMapper elementMapper;
};
}
#endif
