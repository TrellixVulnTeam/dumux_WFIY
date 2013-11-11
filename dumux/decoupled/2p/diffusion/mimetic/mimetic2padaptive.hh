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

#ifndef DUMUX_MIMETIC2PADAPTIVE_HH
#define DUMUX_MIMETIC2PADAPTIVE_HH

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
#include <dumux/common/intersectionmapper.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

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
class MimeticTwoPLocalStiffnessAdaptive: public LocalStiffness<TypeTag, 1>
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
        pressEqIdx = Indices::pressureEqIdx,
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
        pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation), //!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
        saturationType = GET_PROP_VALUE(TypeTag, SaturationFormulation) //!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
    };

    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    typedef Dumux::IntersectionMapper<GridView> IntersectionMapper;

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
    MimeticTwoPLocalStiffnessAdaptive(Problem& problem,  bool levelBoundaryAsDirichlet,
            const GridView& gridView, const IntersectionMapper& intersectionMapper,bool procBoundaryAsDirichlet = true) :
            problem_(problem), gridView_(gridView), intersectionMapper_(intersectionMapper), maxError_(0), timeStep_(1)
    {
        int size = gridView_.size(0);
        rhs_.resize(size , 0.);
        W_.resize(size);
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

        adapt();
    }

    void adapt()
    {
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

    int numberOfFaces(int globalIdx)
    {
        return W_[globalIdx].N();
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
        int numFaces = intersectionMapper_.size(element);

        this->setcurrentsize(numFaces);

        // clear assemble data
        for (int i = 0; i < numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
            for (int j = 0; j < numFaces; j++)
                this->A[i][j] = 0;
        }

        assembleV(element, numFaces, k);
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
        int numFaces = intersectionMapper_.size(element);

        // clear assemble data
        for (int i = 0; i < numFaces; i++)
        {
            this->b[i] = 0;
            this->bctype[i].reset();
        }

        assembleBC(element, k);
    }

    template<class Vector>
    void completeRHS(const Element& element, std::vector<int>& local2Global, Vector& f)
    {
        int globalIdx = problem_.variables().index(element);

        int numFaces = numberOfFaces(globalIdx);

        Dune::DynamicVector<Scalar> F(numFaces);
        Scalar dInv = 0.;
        computeReconstructionMatrices(element, W_[globalIdx], F, dInv);

        for (int i = 0; i < numFaces; i++)
        {
//            if (!this->bctype[i].isSet())
                f[local2Global[i]][0] += (dInv * F[i] * rhs_[globalIdx]);
        }
    }

    Scalar constructPressure(const Element& element, Dune::DynamicVector<Scalar>& pressTrace)
    {
        int globalIdx = problem_.variables().index(element);

        int numFaces = numberOfFaces(globalIdx);

        Scalar volume = element.geometry().volume();

        PrimaryVariables sourceVec(0.0);
        problem_.source(sourceVec, element);
        Scalar qmean = volume * (sourceVec[wPhaseIdx]/density_[wPhaseIdx] + sourceVec[nPhaseIdx]/density_[nPhaseIdx]);
        qmean += rhs_[globalIdx];

        Dune::DynamicVector<Scalar> F(numFaces, 0.);
        Scalar dInv = 0.;
        computeReconstructionMatrices(element, W_[globalIdx], F, dInv);

        return (dInv*(qmean + (F*pressTrace)));
    }

    void constructVelocity(const Element& element, Dune::DynamicVector<Scalar>& vel, Dune::DynamicVector<Scalar>& pressTrace, Scalar press)
    {
        int globalIdx = problem_.variables().index(element);

        int numFaces = numberOfFaces(globalIdx);

        Dune::DynamicVector<Scalar> faceVol(numFaces);

        Dune::FieldVector<Scalar, 2*dim> faceVolumeReal(0.0);

        int idx = 0;
        IntersectionIterator isEndIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
        {
            faceVol[idx] = isIt->geometry().volume();
            int indexInInside = isIt->indexInInside();
            faceVolumeReal[indexInInside] += faceVol[idx];

            idx++;
        }

        vel = 0;
        for (int i = 0; i < numFaces; i++)
            for (int j = 0; j < numFaces; j++)
                vel[i] += W_[globalIdx][i][j]*faceVol[j]*(press - pressTrace[j]);

        for (int i = 0; i < numFaces; i++)
        {
                vel[i] *= faceVol[i]/faceVolumeReal[intersectionMapper_.maplocal(globalIdx, i)];
        }

    }

    void computeReconstructionMatrices(const Element& element, const Dune::DynamicMatrix<Scalar>& W, Dune::DynamicVector<Scalar>& F, Scalar& dInv)
    {
        int globalIdx = problem_.variables().index(element);

        int numFaces = numberOfFaces(globalIdx);

        Dune::DynamicVector<Scalar> c(numFaces);
        Dune::DynamicMatrix<Scalar> Pi(numFaces, numFaces);

        int idx = 0;
        IntersectionIterator isEndIt = gridView_.iend(element);
        for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
        {
            Scalar faceVol = isIt->geometry().volume();

            // Corresponding to the element under consideration,
            // calculate the part of the matrix C coupling velocities and element pressures.
            // This is just a row vector of size numFaces.
            // scale with volume
            c[idx] = faceVol;
            // Set up the element part of the matrix \Pi coupling velocities
            // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
            Pi[idx][idx] = faceVol;
            idx++;
        }

        // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
        Dune::DynamicVector<Scalar> Wc(numFaces);
        W.umv(c, Wc);
        dInv = 1.0 / (c * Wc);

        // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
        F = 0;
        Pi.umv(Wc, F);
    }

    void assembleElementMatrices(const Element& element,
            Dune::DynamicVector<Scalar>& faceVol,
            Dune::DynamicVector<Scalar>& F,
            Dune::DynamicMatrix<Scalar>& Pi,
    Scalar& dInv,
    Scalar& qmean);


    const Problem& problem() const
    {
        return problem_;
    }

private:
    void assembleV(const Element& element, int numFaces, int k = 1);

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
    const GridView& gridView_;
    const IntersectionMapper& intersectionMapper_;
    Dune::DynamicVector<Scalar> rhs_;
    std::vector<Dune::DynamicMatrix<Scalar> > W_;
    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
};

template<class TypeTag>
void MimeticTwoPLocalStiffnessAdaptive<TypeTag>::assembleV(const Element& element, int numFaces, int k)
{
    int globalIdx = problem_.variables().index(element);

    // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
    // The matrix W developed here corresponds to one element-associated
    // block of the matrix B^{-1} there.
    W_[globalIdx].resize(numFaces, numFaces);
    Dune::DynamicVector<Scalar> faceVol(numFaces);
    Dune::DynamicMatrix<Scalar> Pi(numFaces, numFaces);
    Dune::DynamicVector<Scalar> F(numFaces);
    Scalar dInv = 0.;
    Scalar qmean = 0.;
    this->assembleElementMatrices(element, faceVol, F, Pi, dInv, qmean);

    // Calculate the element part of the matrix Pi W Pi^T.
    Dune::DynamicMatrix<Scalar> PiWPiT(W_[globalIdx]);
    PiWPiT.rightmultiply(Pi);
    PiWPiT.leftmultiply(Pi);

    // Calculate the element part of the matrix F D^{-1} F^T.
    Dune::DynamicMatrix<Scalar> FDinvFT(numFaces, numFaces);
    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < numFaces; j++)
            FDinvFT[i][j] = dInv * F[i] * F[j];

    // Calculate the element part of the matrix S = Pi W Pi^T - F D^{-1} F^T.
    for (int i = 0; i < numFaces; i++)
        for (int j = 0; j < numFaces; j++)
            this->A[i][j] = PiWPiT[i][j] - FDinvFT[i][j];

    // Calculate the source term F D^{-1} f
    // NOT WORKING AT THE MOMENT
    Scalar factor = dInv * qmean;
    for (int i = 0; i < numFaces; i++)
        this->b[i] = F[i] * factor;

    //        std::cout << "faceVol = " << faceVol << std::endl << "W = " << std::endl << W << std::endl
    //              << "c = " << c << std::endl << "Pi = " << std::endl << Pi << std::endl
    //              << "dinv = " << dinv << std::endl << "F = " << F << std::endl;
    //        std::cout << "dinvF = " << dinvF << ", q = " << qmean
    //             << ", b = " << this->b[0] << ", " << this->b[1] << ", " << this->b[2] << ", " << this->b[3] << std::endl;
}

template<class TypeTag>
void MimeticTwoPLocalStiffnessAdaptive<TypeTag>::assembleElementMatrices(const Element& element,
        Dune::DynamicVector<Scalar>& faceVol,
        Dune::DynamicVector<Scalar>& F,
        Dune::DynamicMatrix<Scalar>& Pi,
        Scalar& dInv,
        Scalar& qmean)
{
    // get global coordinate of cell center
    const Dune::FieldVector<Scalar, dim>& centerGlobal = element.geometry().center();

    int globalIdx = problem_.variables().index(element);

    CellData& cellData = problem_.variables().cellData(globalIdx);

    int numFaces = intersectionMapper_.size(globalIdx);

    // cell volume
    Scalar volume = element.geometry().volume();

    // build the matrices R and ~N
    Dune::DynamicMatrix<Scalar> R(numFaces, dim, 0.);
    Dune::DynamicMatrix<Scalar> N(numFaces, dim, 0.);

    //       std::cout << "element " << elemId << ": center " << centerGlobal << std::endl;;

    //collect information needed for calculation of fluxes due to capillary-potential (pc + gravity!)
    std::vector<Scalar> pcPotFace(numFaces);
    Scalar pcPot = (problem_.bBoxMax() - element.geometry().center()) * problem_.gravity() * (density_[nPhaseIdx] - density_[wPhaseIdx]);

    int idx = 0;
    IntersectionIterator isEndIt = gridView_.iend(element);
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
    {
        // local number of facet

        Dune::FieldVector<Scalar, dim> faceGlobal(isIt->geometry().center());
        faceVol[idx] = isIt->geometry().volume();

        // get normal vector
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = isIt->centerUnitOuterNormal();


        for (int k = 0; k < dim; k++)
        {
            // move origin to the center of gravity
            R[idx][k] = faceVol[idx] * (faceGlobal[k] - centerGlobal[k]);
            N[idx][k] = unitOuterNormal[k];
        }

        pcPotFace[idx] = (problem_.bBoxMax() - faceGlobal) * problem_.gravity() * (density_[nPhaseIdx] - density_[wPhaseIdx]);

        idx++;
    }

    // proceed along the lines of Algorithm 1 from
    // Brezzi/Lipnikov/Simonicini M3AS 2005
    // (1) orthonormalize columns of the matrix R (Gram-Schmidt orthonormalization)
    Scalar norm = R[0][0] * R[0][0];
    for (int i = 1; i < numFaces; i++)
        norm += R[i][0] * R[i][0];
    norm = sqrt(norm);
    for (int i = 0; i < numFaces; i++)
        R[i][0] /= norm;
    Scalar weight = R[0][1] * R[0][0];
    for (int i = 1; i < numFaces; i++)
        weight += R[i][1] * R[i][0];
    for (int i = 0; i < numFaces; i++)
        R[i][1] -= weight * R[i][0];
    norm = R[0][1] * R[0][1];
    for (int i = 1; i < numFaces; i++)
        norm += R[i][1] * R[i][1];
    norm = sqrt(norm);
    for (int i = 0; i < numFaces; i++)
        R[i][1] /= norm;
    if (dim == 3)
    {
        Scalar weight1 = R[0][2] * R[0][0];
        Scalar weight2 = R[0][2] * R[0][1];
        for (int i = 1; i < numFaces; i++)
        {
            weight1 += R[i][2] * R[i][0];
            weight2 += R[i][2] * R[i][1];
        }
        for (int i = 0; i < numFaces; i++)
            R[i][2] -= weight1 * R[i][0] + weight2 * R[i][1];
        norm = R[0][2] * R[0][2];
        for (int i = 1; i < numFaces; i++)
            norm += R[i][2] * R[i][2];
        norm = sqrt(norm);
        for (int i = 0; i < numFaces; i++)
            R[i][2] /= norm;
    }
    //      std::cout << "~R =\dim" << R;

    // (2) Build the matrix ~D
    Dune::DynamicMatrix<Scalar> D(numFaces, numFaces, 0.);
    for (int s = 0; s < numFaces; s++)
    {
        Dune::DynamicVector<Scalar> es(numFaces, 0.);
        es[s] = 1;
        for (int k = 0; k < numFaces; k++)
        {
            D[k][s] = es[k];
            for (int i = 0; i < dim; i++)
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
    Dune::DynamicMatrix<Scalar> NK(N);
    NK.rightmultiply(mobility);

//    std::cout<<"NK = "<<NK<<"\n";

    for (int i = 0; i < numFaces; i++)
    {
        for (int j = 0; j < numFaces; j++)
        {
            W_[globalIdx][i][j] = NK[i][0] * N[j][0];
            for (int k = 1; k < dim; k++)
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
    Dune::DynamicVector<Scalar> c(numFaces);
    for (int i = 0; i < numFaces; i++)
        c[i] = faceVol[i];

    // Set up the element part of the matrix \Pi coupling velocities
    // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
    for (int i = 0; i < numFaces; i++)
        Pi[i][i] = faceVol[i];

    // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
    Dune::DynamicVector<Scalar> Wc(numFaces);
    W_[globalIdx].umv(c, Wc);
    dInv = 1.0 / (c * Wc);

    // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
    F = 0;
    Pi.umv(Wc, F);
    //      std::cout << "Pi = \dim" << Pi << "c = " << c << ", F = " << F << std::endl;

    //accumulate fluxes due to capillary potential (pc + gravity!)
    idx = 0;
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != isEndIt; ++isIt)
    {
            Scalar fracFlow = 0;

            Scalar flux = 0;
            for (int j = 0; j < numFaces; j++)
                flux += W_[globalIdx][idx][j] * faceVol[j] * (pcPot - pcPotFace[j]);

            if (isIt->neighbor())
            {

            ElementPointer neighbor = isIt->outside();
            int globalIdxNeighbor = problem_.variables().index(*neighbor);
            if (flux > 0.)
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

        if (bctype.isDirichlet(pressEqIdx))
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
        idx++;
    }

    PrimaryVariables sourceVec(0.0);
    problem_.source(sourceVec, element);
    qmean = volume * (sourceVec[wPhaseIdx]/density_[wPhaseIdx] + sourceVec[nPhaseIdx]/density_[nPhaseIdx]);

    qmean += evaluateErrorTerm(cellData) * volume;
}

template<class TypeTag>
void MimeticTwoPLocalStiffnessAdaptive<TypeTag>::assembleBC(const Element& element, int k)
{
    // evaluate boundary conditions via intersection iterator
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    unsigned int faceIndex = 0;
    IntersectionIterator endIsIt = gridView_.iend(element);
    for (IntersectionIterator isIt = gridView_.ibegin(element); isIt != endIsIt; ++isIt)
    {
        if (isIt->neighbor())
        {

        }
        else if (isIt->boundary())
        {
            problem_.boundaryTypes(this->bctype[faceIndex], *isIt);
            PrimaryVariables boundValues(0.0);

            if (this->bctype[faceIndex].isNeumann(pressEqIdx))
            {
                problem_.neumann(boundValues, *isIt);
                Scalar J = (boundValues[wPhaseIdx]/density_[wPhaseIdx] + boundValues[nPhaseIdx]/density_[nPhaseIdx]);
                this->b[faceIndex] -= J * isIt->geometry().volume();
            }
            else if (this->bctype[faceIndex].isDirichlet(pressEqIdx))
            {
                problem_.dirichlet(boundValues, *isIt);
                if (pressureType == pw)
                    this->b[faceIndex] = boundValues[pressureIdx] + (problem_.bBoxMax() - isIt->geometry().center()) * problem_.gravity() * density_[wPhaseIdx];
                else if (pressureType == pn)
                    this->b[faceIndex] = boundValues[pressureIdx] + (problem_.bBoxMax() - isIt->geometry().center()) * problem_.gravity() * density_[nPhaseIdx];
                else
                    this->b[faceIndex] = boundValues[pressureIdx];
            }
        }
        faceIndex++;
    }
}

/** @} */
}
#endif
