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

#ifndef DUMUX_MIMETIC2P_HH
#define DUMUX_MIMETIC2P_HH

/*!
 * \file
 *
 * \brief Local stiffness matrix for the diffusion equation discretized by mimetic FD
 */

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include<dumux/common/boundaryconditions.hh>
#include "localstiffness.hh"

#include <dune/common/dynvector.hh>

namespace Dumux
{
/*!
 *  \ingroup Mimetic2p
 */
/**
 * @brief compute local stiffness matrix for conforming finite elements for the full 2-phase pressure equation
 *
 */

//! A class for computing local stiffness matrices
/*! A class for computing local stiffness matrix for the full 2-phase pressure equation
 */
template<class TypeTag>
class MimeticTwoPLocalStiffness: public LocalStiffness<TypeTag, 1>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // grid types
    enum
    {
        dim = GridView::dimension
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressureEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNw,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
        pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation),
        //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
        saturationType = GET_PROP_VALUE(TypeTag, SaturationFormulation)
    };

    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum
    {
        m = 1
    };
    enum
    {
        size = 2 * dim
    };

    //! Constructor
    MimeticTwoPLocalStiffness(Problem& problem,  bool levelBoundaryAsDirichlet,
            const GridView& gridView, bool procBoundaryAsDirichlet = true) :
            problem_(problem), gridView_(gridView), maxError_(0), timeStep_(1)
    {
        ErrorTermFactor_ = GET_PARAM(TypeTag, Scalar, ImpetErrorTermFactor);
        ErrorTermLowerBound_ = GET_PARAM(TypeTag, Scalar, ImpetErrorTermLowerBound);
        ErrorTermUpperBound_ = GET_PARAM(TypeTag, Scalar, ImpetErrorTermUpperBound);

        density_[wPhaseIdx] = 0.0;
        density_[nPhaseIdx] = 0.0;
        viscosity_[wPhaseIdx] =0.0;
        viscosity_[nPhaseIdx] = 0.0;
    }

    void initialize()
    {
        ElementIterator element = problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

        int size = gridView_.size(0);
        rhs_.resize(size , 0.);
        W_.resize(size);
    }

    void reset()
    {
        rhs_ = 0;
        maxError_ = 0;
        timeStep_ = 1;
    }

    void setErrorInfo(Scalar maxErr, Scalar dt)
    {
        maxError_ = maxErr;
        timeStep_ = dt;
    }

    //! assemble local stiffness matrix for given element and order
    /*! On exit the following things have been done:
     - The stiffness matrix for the given entity and polynomial degree has been assembled and is
     accessible with the mat() method.
     - The boundary conditions have been evaluated and are accessible with the bc() method
     - The right hand side has been assembled. It contains either the value of the essential boundary
     condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
     @param[in]  element a codim 0 entity reference
     @param[in]  k order of CR basis (only k = 1 is implemented)
     */
    void assemble(const Element& element, int k = 1)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // clear assemble data
        for (unsigned int i = 0; i < numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
            for (unsigned int j = 0; j < numFaces; j++)
                this->A[i][j] = 0;
        }

        assembleV(element, k);
        assembleBC(element, k);
    }

    // TODO/FIXME: this is only valid for linear problems where
    // the local stiffness matrix is independend of the current
    // solution. We need to implement this properly, but this
    // should at least make the thing compile...
    typedef Dune::FieldVector<Scalar, m> VBlockType;
    void assemble(const Element &cell, const Dune::BlockVector<VBlockType>& localSolution, int orderOfShapeFns = 1)
    {
        assemble(cell, orderOfShapeFns);
    }

    //! assemble only boundary conditions for given element
    /*! On exit the following things have been done:
     - The boundary conditions have been evaluated and are accessible with the bc() method
     - The right hand side contains either the value of the essential boundary
     condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
     @param[in]  element a codim 0 entity reference
     @param[in]  k order of CR basis
     */
    void assembleBoundaryCondition(const Element& element, int k = 1)
    {
        unsigned int numFaces = element.template count<1>();
        this->setcurrentsize(numFaces);

        // clear assemble data
        for (unsigned int i = 0; i < numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
        }

        assembleBC(element, k);
    }

    template<class Vector>
    void completeRHS(const Element& element, Dune::FieldVector<int, 2*dim>& local2Global, Vector& f)
    {
        int globalIdx = problem_.variables().index(element);
        unsigned int numFaces = element.template count<1>();

        Dune::FieldVector<Scalar, 2 * dim> F(0.);
        Scalar dInv = 0.;
        computeReconstructionMatrices(element, W_[globalIdx], F, dInv);

        for (unsigned int i = 0; i < numFaces; i++)
        {
            if (!this->bc(i).isDirichlet(pressureEqIdx))
                f[local2Global[i]][0] += (dInv * F[i] * rhs_[globalIdx]);
        }
    }

    Scalar constructPressure(const Element& element, Dune::FieldVector<Scalar,2*dim>& pressTrace)
    {
        int globalIdx = problem_.variables().index(element);
        Scalar volume = element.geometry().volume();

        PrimaryVariables sourceVec(0.0);
        problem_.source(sourceVec, element);
        Scalar qmean = volume * (sourceVec[wPhaseIdx]/density_[wPhaseIdx] + sourceVec[nPhaseIdx]/density_[nPhaseIdx]);
        qmean += rhs_[globalIdx];

        Dune::FieldVector<Scalar, 2 * dim> F(0.);
        Scalar dInv = 0.;
        computeReconstructionMatrices(element, W_[globalIdx], F, dInv);

        return (dInv*(qmean + (F*pressTrace)));
    }

    void constructVelocity(const Element& element, Dune::FieldVector<Scalar,2*dim>& vel,
                           Dune::FieldVector<Scalar,2*dim>& pressTrace, Scalar press)
    {
        int globalIdx = problem_.variables().index(element);

        Dune::FieldVector<Scalar, 2 * dim> faceVol(0);
        IntersectionIterator isEndIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
        {
            faceVol[isIt->indexInInside()] = isIt->geometry().volume();
        }

        vel = 0;
        for (int i = 0; i < 2*dim; i++)
            for (int j = 0; j < 2*dim; j++)
                vel[i] += W_[globalIdx][i][j]*faceVol[j]*(press - pressTrace[j]);
    }

    void constructVelocity(const Element& element, int faceIdx, Scalar& vel,
                           Dune::FieldVector<Scalar,2*dim>& pressTrace, Scalar press)
    {
        int globalIdx = problem_.variables().index(element);

        Dune::FieldVector<Scalar, 2 * dim> faceVol(0);
        IntersectionIterator isEndIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
        {
            faceVol[isIt->indexInInside()] = isIt->geometry().volume();
        }

        vel = 0;
            for (int j = 0; j < 2*dim; j++)
                vel += W_[globalIdx][faceIdx][j]*faceVol[j]*(press - pressTrace[j]);
    }

    void computeReconstructionMatrices(const Element& element, 
    								   const Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim>& W, 
    								   Dune::FieldVector<Scalar, 2 * dim>& F, Scalar& dInv)
    {
        Dune::FieldVector<Scalar, 2 * dim> c(0);
        Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> Pi(0);

        IntersectionIterator isEndIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
        {
            // local number of facet
            int idx = isIt->indexInInside();

            Scalar faceVol = isIt->geometry().volume();

            // Corresponding to the element under consideration,
            // calculate the part of the matrix C coupling velocities and element pressures.
            // This is just a row vector of size numFaces.
            // scale with volume
            c[idx] = faceVol;
            // Set up the element part of the matrix \Pi coupling velocities
            // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
            Pi[idx][idx] = faceVol;
        }

        // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
        Dune::FieldVector<Scalar, 2 * dim> Wc(0);
        W.umv(c, Wc);
        dInv = 1.0 / (c * Wc);

        // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
        F = 0;
        Pi.umv(Wc, F);
    }

    void assembleElementMatrices(const Element& element, Dune::FieldVector<Scalar, 2 * dim>& faceVol,
    Dune::FieldVector<Scalar, 2 * dim>& F,
    Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim>& Pi,
    Scalar& dInv,
    Scalar& qmean);


    const Problem& problem() const
    {
        return problem_;
    }

private:
    void assembleV(const Element& element, int k = 1);

    void assembleBC(const Element& element, int k = 1);

    Scalar evaluateErrorTerm(CellData& cellData)
    {
        //error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        // error reduction routine: volumetric error is damped and inserted to right hand side
        Scalar sat = 0;
        switch (saturationType)
        {
        case Sw:
            sat = cellData.saturation(wPhaseIdx);
            break;
        case Sn:
            sat = cellData.saturation(nPhaseIdx);
            break;
        }

        Scalar error = (sat > 1.0) ? sat - 1.0 : 0.0;
        if (sat < 0.0)
        {
            error = sat;
        }
        error /= timeStep_;

        Scalar errorAbs = std::abs(error);

        if ((errorAbs * timeStep_ > 1e-6) && (errorAbs > ErrorTermLowerBound_ * maxError_)
                && (!problem_.timeManager().willBeFinished()))
        {
            return ErrorTermFactor_ * error;
        }
        return 0.0;
    }

private:
    Problem& problem_;
    GridView gridView_;
    Dune::DynamicVector<Scalar> rhs_;
    std::vector<Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> > W_;
    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
};

template<class TypeTag>
void MimeticTwoPLocalStiffness<TypeTag>::assembleV(const Element& element, int)
{
    unsigned int numFaces = element.template count<1>();
    this->setcurrentsize(numFaces);

    int globalIdx = problem_.variables().index(element);

    // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
    // The matrix W developed here corresponds to one element-associated
    // block of the matrix B^{-1} there.
    Dune::FieldVector<Scalar, 2 * dim> faceVol(0);
    Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> Pi(0.);
    Dune::FieldVector<Scalar, 2 * dim> F(0.);
    Scalar dInv = 0.;
    Scalar qmean = 0.;
    this->assembleElementMatrices(element, faceVol, F, Pi, dInv, qmean);

    // Calculate the element part of the matrix Pi W Pi^T.
    Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> PiWPiT(W_[globalIdx]);
    PiWPiT.rightmultiply(Pi);
    PiWPiT.leftmultiply(Pi);

    // Calculate the element part of the matrix F D^{-1} F^T.
    Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> FDinvFT(0);
    for (unsigned int i = 0; i < numFaces; i++)
        for (unsigned int j = 0; j < numFaces; j++)
            FDinvFT[i][j] = dInv * F[i] * F[j];

    // Calculate the element part of the matrix S = Pi W Pi^T - F D^{-1} F^T.
    for (unsigned int i = 0; i < numFaces; i++)
        for (unsigned int j = 0; j < numFaces; j++)
            this->A[i][j] = PiWPiT[i][j] - FDinvFT[i][j];

    // Calculate the source term F D^{-1} f
    // NOT WORKING AT THE MOMENT
    Scalar factor = dInv * qmean;
    for (unsigned int i = 0; i < numFaces; i++)
        this->b[i] = F[i] * factor;

    //        std::cout << "faceVol = " << faceVol << std::endl << "W = " << std::endl << W << std::endl
    //              << "c = " << c << std::endl << "Pi = " << std::endl << Pi << std::endl
    //              << "dinv = " << dinv << std::endl << "F = " << F << std::endl;
    //        std::cout << "dinvF = " << dinvF << ", q = " << qmean
    //             << ", b = " << this->b[0] << ", " << this->b[1] << ", " << this->b[2] << ", " << this->b[3] << std::endl;
}

template<class TypeTag>
void MimeticTwoPLocalStiffness<TypeTag>::assembleElementMatrices(const Element& element,
        Dune::FieldVector<Scalar, 2 * dim>& faceVol,
        Dune::FieldVector<Scalar, 2 * dim>& F,
        Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim>& Pi,
        Scalar& dInv,
        Scalar& qmean)
{
    unsigned int numFaces = element.template count<1>();
    this->setcurrentsize(numFaces);

    // get global coordinate of cell center
    Dune::FieldVector<Scalar, dim> centerGlobal = element.geometry().center();

    int globalIdx = problem_.variables().index(element);

    CellData& cellData = problem_.variables().cellData(globalIdx);

    // cell volume
    Scalar volume = element.geometry().volume();

    // build the matrices R and ~N
    Dune::FieldMatrix<Scalar, 2 * dim, dim> R(0), N(0);

    //       std::cout << "element " << elemId << ": center " << centerGlobal << std::endl;

    //collect information needed for calculation of fluxes due to capillary-potential (pc + gravity!)
    Scalar gravPotFace[2*dim];

    Scalar gravPot = (problem_.bBoxMax() - centerGlobal) * problem_.gravity() * (density_[nPhaseIdx] - density_[wPhaseIdx]);

    int i = -1;
    IntersectionIterator isEndIt = gridView_.iend(element);
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
    {
        // local number of facet
        i = isIt->indexInInside();

        Dune::FieldVector<Scalar, dim> faceGlobal = isIt->geometry().center();
        faceVol[i] = isIt->geometry().volume();

        // get normal vector
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = isIt->centerUnitOuterNormal();

        N[i] = unitOuterNormal;

        for (int k = 0; k < dim; k++)
            // move origin to the center of gravity
            R[i][k] = faceVol[i] * (faceGlobal[k] - centerGlobal[k]);

        gravPotFace[i] = (problem_.bBoxMax() - faceGlobal) * problem_.gravity() * (density_[nPhaseIdx] - density_[wPhaseIdx]);
    }

    // proceed along the lines of Algorithm 1 from
    // Brezzi/Lipnikov/Simonicini M3AS 2005
    // (1) orthonormalize columns of the matrix R
    Scalar norm = R[0][0] * R[0][0];
    for (unsigned int i = 1; i < numFaces; i++)
        norm += R[i][0] * R[i][0];
    norm = sqrt(norm);
    for (unsigned int i = 0; i < numFaces; i++)
        R[i][0] /= norm;
    Scalar weight = R[0][1] * R[0][0];
    for (unsigned int i = 1; i < numFaces; i++)
        weight += R[i][1] * R[i][0];
    for (unsigned int i = 0; i < numFaces; i++)
        R[i][1] -= weight * R[i][0];
    norm = R[0][1] * R[0][1];
    for (unsigned int i = 1; i < numFaces; i++)
        norm += R[i][1] * R[i][1];
    norm = sqrt(norm);
    for (unsigned int i = 0; i < numFaces; i++)
        R[i][1] /= norm;
    if (dim == 3)
    {
        Scalar weight1 = R[0][2] * R[0][0];
        Scalar weight2 = R[0][2] * R[0][1];
        for (unsigned int i = 1; i < numFaces; i++)
        {
            weight1 += R[i][2] * R[i][0];
            weight2 += R[i][2] * R[i][1];
        }
        for (unsigned int i = 0; i < numFaces; i++)
            R[i][2] -= weight1 * R[i][0] + weight2 * R[i][1];
        norm = R[0][2] * R[0][2];
        for (unsigned int i = 1; i < numFaces; i++)
            norm += R[i][2] * R[i][2];
        norm = sqrt(norm);
        for (unsigned int i = 0; i < numFaces; i++)
            R[i][2] /= norm;
    }
    //      std::cout << "~R =\dim" << R;

    // (2) Build the matrix ~D
    Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> D(0);
    for (unsigned int s = 0; s < numFaces; s++)
    {
        Dune::FieldVector<Scalar, 2 * dim> es(0);
        es[s] = 1;
        for (unsigned int k = 0; k < numFaces; k++)
        {
            D[k][s] = es[k];
            for (unsigned int i = 0; i < dim; i++)
            {
                D[k][s] -= R[s][i] * R[k][i];
            }
        }
    }

//    std::cout<<"D = "<<D<<"\n";

    // eval diffusion tensor, ASSUMING to be constant over each cell
    Dune::FieldMatrix<Scalar, dim, dim> mobility(problem_.spatialParams().intrinsicPermeability(element));
    mobility *= (cellData.mobility(wPhaseIdx) + cellData.mobility(nPhaseIdx));

    Scalar traceMob = mobility[0][0];
    for (int i = 1; i < dim; i++)
        traceMob += mobility[i][i];
    D *= 2.0 * traceMob / volume;
    //      std::cout << "u~D =\dim" << D;

    // (3) Build the matrix W = Minv
    Dune::FieldMatrix<Scalar, 2 * dim, dim> NK(N);
    NK.rightmultiply(mobility);

    for (unsigned int i = 0; i < numFaces; i++)
    {
        for (unsigned int j = 0; j < numFaces; j++)
        {
            W_[globalIdx][i][j] = NK[i][0] * N[j][0];
            for (unsigned int k = 1; k < dim; k++)
                W_[globalIdx][i][j] += NK[i][k] * N[j][k];
        }
    }

    W_[globalIdx] /= volume;
    W_[globalIdx] += D;

    //       std::cout << "W = \dim" << W;

    // Now the notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
    // The matrix W developed so far corresponds to one element-associated
    // block of the matrix B^{-1} there.

    // Corresponding to the element under consideration,
    // calculate the part of the matrix C coupling velocities and element pressures.
    // This is just a row vector of size numFaces.
    // scale with volume
    Dune::FieldVector<Scalar, 2 * dim> c(0);
    for (unsigned int i = 0; i < numFaces; i++)
        c[i] = faceVol[i];

    // Set up the element part of the matrix \Pi coupling velocities
    // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
    for (unsigned int i = 0; i < numFaces; i++)
        Pi[i][i] = faceVol[i];

    // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
    Dune::FieldVector<Scalar, 2 * dim> Wc(0);
    W_[globalIdx].umv(c, Wc);
    dInv = 1.0 / (c * Wc);

    // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
    F = 0;
    Pi.umv(Wc, F);
    //      std::cout << "Pi = \dim" << Pi << "c = " << c << ", F = " << F << std::endl;

    //accumulate fluxes due to capillary potential (pc + gravity!)
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
    {
        int idx = isIt->indexInInside();

        Scalar fracFlow = 0;

        Scalar flux = 0;
        for (int j = 0; j < 2 * dim; j++)
            flux += W_[globalIdx][idx][j] * faceVol[j] * (gravPot - gravPotFace[j]);

        //it is enough to evaluate the capillary/gravity flux for neighbors -> not needed anyway at the boundaries!
        if (isIt->neighbor())
        {
            ElementPointer neighbor = isIt->outside();
            int globalIdxNeighbor = problem_.variables().index(*neighbor);
            if (flux >= 0.)
            {
                switch (pressureType)
                {
                case pw:
                {
                    fracFlow = cellData.fracFlowFunc(nPhaseIdx);
                    break;
                }
                case pn:
                {
                    fracFlow = -cellData.fracFlowFunc(wPhaseIdx);
                    break;
                }
                }


                rhs_[globalIdx] -= (faceVol[idx] * fracFlow * flux);
                rhs_[globalIdxNeighbor] += (faceVol[idx] * fracFlow * flux);
            }
        }
        else if (isIt->boundary())
        {
            BoundaryTypes bctype;
            problem_.boundaryTypes(bctype, *isIt);

            if (bctype.isDirichlet(pressureEqIdx))
            {
                if (flux > 0. || !bctype.isDirichlet(satEqIdx))
                {
                    switch (pressureType)
                    {
                    case pw:
                    {
                        fracFlow = cellData.fracFlowFunc(nPhaseIdx);
                        break;
                    }
                    case pn:
                    {
                        fracFlow = -cellData.fracFlowFunc(wPhaseIdx);
                        break;
                    }
                    }
                }
                else if (flux < 0. && bctype.isDirichlet(satEqIdx))
                {
                    PrimaryVariables boundValues(0.0);
                    problem_.dirichlet(boundValues, *isIt);

                    Scalar krw = MaterialLaw::krw(problem_.spatialParams().materialLawParams(element),
                            boundValues[saturationIdx]);
                    Scalar krn = MaterialLaw::krn(problem_.spatialParams().materialLawParams(element),
                            boundValues[saturationIdx]);

                    switch (pressureType)
                    {
                    case pw:
                    {
                        fracFlow = krn / viscosity_[nPhaseIdx]
                                / (krw / viscosity_[wPhaseIdx] + krn / viscosity_[nPhaseIdx]);
                        break;
                    }
                    case pn:
                    {
                        fracFlow = -krw / viscosity_[wPhaseIdx]
                                / (krw / viscosity_[wPhaseIdx] + krn / viscosity_[nPhaseIdx]);
                        break;
                    }
                    }
                }

                rhs_[globalIdx] -= (faceVol[idx] * fracFlow * flux);
            }
        }
    }

    PrimaryVariables sourceVec(0.0);
    problem_.source(sourceVec, element);
    qmean = volume * (sourceVec[wPhaseIdx]/density_[wPhaseIdx] + sourceVec[nPhaseIdx]/density_[nPhaseIdx]);

        qmean += evaluateErrorTerm(cellData) * volume;
}

template<class TypeTag>
void MimeticTwoPLocalStiffness<TypeTag>::assembleBC(const Element& element, int k)
{
    // evaluate boundary conditions via intersection iterator
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    IntersectionIterator endIsIt = gridView_.iend(element);
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != endIsIt; ++isIt)
    {
        int faceIndex = isIt->indexInInside();
        if (isIt->boundary())
        {
            problem_.boundaryTypes(this->bctype[faceIndex], *isIt);
            PrimaryVariables boundValues(0.0);

            if (this->bctype[faceIndex].isNeumann(pressureEqIdx))
            {
                problem_.neumann(boundValues, *isIt);
                Scalar J = (boundValues[wPhaseIdx]/density_[wPhaseIdx] + boundValues[nPhaseIdx]/density_[nPhaseIdx]);
                this->b[faceIndex] -= J * isIt->geometry().volume();
            }
            else if (this->bctype[faceIndex].isDirichlet(pressureEqIdx))
            {
                problem_.dirichlet(boundValues, *isIt);
                if (pressureType == pw)
                    this->b[faceIndex] = boundValues[pressureIdx] +
                    (problem_.bBoxMax() - isIt->geometry().center()) * problem_.gravity() * density_[wPhaseIdx];
                else if (pressureType == pn)
                    this->b[faceIndex] = boundValues[pressureIdx] +
                    (problem_.bBoxMax() - isIt->geometry().center()) * problem_.gravity() * density_[nPhaseIdx];
                else
                    this->b[faceIndex] = boundValues[pressureIdx];
            }
        }
    }

}

/** @} */
}
#endif
