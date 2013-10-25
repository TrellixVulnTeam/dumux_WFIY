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
 * \brief spatial parameters for the test problem for diffusion models.
 */
#ifndef TEST_1P_SPATIALPARAMS_HH
#define TEST_1P_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fvspatialparams1p.hh>

namespace Dumux
{

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for 1-p diffusion models.
 */
template<class TypeTag>
class TestOnePSpatialParams: public FVSpatialParamsOneP<TypeTag>
{
    typedef FVSpatialParamsOneP<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:

    const FieldMatrix& intrinsicPermeability (const Element& element) const
    {
        return permeability_[indexSet_.index(element)];
    }

    double porosity(const Element& element) const
    {
        return 0.2;
    }

    void initialize(const double delta)
    {
        delta_ = delta;
        permeability_.resize(gridView_.size(0));

        ElementIterator eIt = gridView_.template begin<0>();
        ElementIterator eEndIt = gridView_.template end<0>();
        for(;eIt != eEndIt; ++eIt)
        {
            setPermeability_(permeability_[indexSet_.index(*eIt)], eIt->geometry().center());
        }

    }

    TestOnePSpatialParams(const GridView& gridView)
    : ParentType(gridView), gridView_(gridView), indexSet_(gridView.indexSet())
    { }

private:
    void setPermeability_(FieldMatrix& perm, const GlobalPosition& globalPos) const
    {
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        perm[0][0] = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        perm[0][1] = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        perm[1][0] = perm[0][1];
        perm[1][1] = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
    }

    const GridView gridView_;
    const IndexSet& indexSet_;
    std::vector<FieldMatrix> permeability_;
    double delta_;
};

} // end namespace
#endif
