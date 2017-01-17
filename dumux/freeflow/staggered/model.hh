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
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
 */

#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

// #include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup NavierStokesModel
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */
template<class TypeTag >
class NavierStokesModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry) GlobalFVGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:



    /*!
     * \brief \copybrief Dumux::ImplicitModel::addOutputVtkFields
     *
     * Specialization for the NavierStokesModel, adding the pressure and
     * the process rank to the VTK writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        // TODO: implement vtk output properly, account for 3d
        using  VectorField = Dune::BlockVector<Dune::FieldVector<double, dimWorld> >;

        // create the required scalar fields
        const auto numElements = this->gridView_().size(0);
        auto *p = writer.allocateManagedBuffer(numElements);
        auto *delP = writer.allocateManagedBuffer(numElements);

        auto *v_x_pos = writer.allocateManagedBuffer(numElements);
        auto *v_x_neg = writer.allocateManagedBuffer(numElements);
        auto *v_y_pos = writer.allocateManagedBuffer(numElements);
        auto *v_y_neg = writer.allocateManagedBuffer(numElements);

        VectorField *velocity = writer.template allocateManagedBuffer<double, dimWorld>(numElements);

        // initialize velocity field
        for (unsigned int i = 0; i < numElements; ++i)
        {
            (*velocity)[i] = Scalar(0);
        }

       auto *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            auto eIdx = this->elementMapper().index(element);
            (*rank)[eIdx] = this->gridView_().comm().rank();

           // get the local fv geometry
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(this->curGlobalVolVars());
            elemVolVars.bindElement(element, fvGeometry, this->curSol());

            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                auto dofIdxGlobal = scv.dofIndex();

                (*p)[dofIdxGlobal] = volVars.pressure();
                (*delP)[dofIdxGlobal] = volVars.pressure() - 1.1e5;

                GlobalPosition velocityVector(0.0);
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    auto& origFaceVars = this->curGlobalFaceVars().faceVars(scvf.dofIndexSelf());
                    auto dirIdx = scvf.directionIndex();

                    velocityVector[dirIdx] += 0.5*origFaceVars.velocity();

                    if(scvf.unitOuterNormal()[dirIdx] > 0.0)
                    {
                        if(dirIdx == 0)
                            (*v_x_pos)[dofIdxGlobal] = origFaceVars.velocity();
                        if(dirIdx == 1)
                            (*v_y_pos)[dofIdxGlobal] = origFaceVars.velocity();
                    }
                    else
                    {
                        if(dirIdx == 0)
                            (*v_x_neg)[dofIdxGlobal] = origFaceVars.velocity();
                        if(dirIdx == 1)
                            (*v_y_neg)[dofIdxGlobal] = origFaceVars.velocity();
                    }
                }
                (*velocity)[dofIdxGlobal] = velocityVector;
            }
        }
        writer.attachDofData(*p, "p", isBox);
        writer.attachDofData(*delP, "delP", isBox);
        writer.attachDofData(*v_x_pos, "v_x_pos", isBox);
        writer.attachDofData(*v_x_neg, "v_x_neg", isBox);
        writer.attachDofData(*v_y_pos, "v_y_pos", isBox);
        writer.attachDofData(*v_y_neg, "v_y_neg", isBox);
        writer.attachCellData(*rank, "process rank");
        writer.attachDofData(*velocity,  "velocity", isBox, dim);
    }

    auto velocity(const Element& element) const
    {
        GlobalPosition velocityVector(0.0);

        // get the local fv geometry
        auto fvGeometry = localView(this->globalFvGeometry());
        fvGeometry.bindElement(element);
        for (auto&& scv : scvs(fvGeometry))
        {
            for (auto&& scvf : scvfs(fvGeometry))
            {
                auto& origFaceVars = this->curGlobalFaceVars().faceVars(scvf.dofIndexSelf());
                auto dirIdx = scvf.directionIndex();
                velocityVector[dirIdx] += 0.5*origFaceVars.velocity();
            }
        }
        return velocityVector;
    }


};
}

#include "propertydefaults.hh"

#endif
