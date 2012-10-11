// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_FVPRESSURE2P2C_ADAPTIVE_HH
#define DUMUX_FVPRESSURE2P2C_ADAPTIVE_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
// include 2p mpfa pressure model
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal2pfaboundpressure2padaptive.hh>

#include <dumux/decoupled/2p2c/fvpressure2p2c.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Benjamin Faigle, Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{
namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(GridAdaptEnableMultiPointFluxApproximation);
NEW_PROP_TAG(GridAdaptEnableSecondHalfEdge);
}

//! The finite volume model for the solution of the compositional pressure equation
/*! \ingroup Adaptive2p2c
 *  Provides a Finite Volume implementation for the pressure equation of a compressible
 *  system with two components. An IMPES-like method is used for the sequential
 *  solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 *  Isothermal conditions and local thermodynamic
 *  equilibrium are assumed.  Gravity is included.
 *  \f[
         c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
          = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 *  \f]
 *  where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 *  \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * This adaptive implementation uses its own initialization and assembling methods for matrix and right-hand-side vector,
 * but solution is done by the base class FVPressure. Fluxes near hanging nodes can be
 * calculated using an \a mpfa method.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2P2CAdaptive
: public FVPressure2P2C<TypeTag>
{
    //the model implementation
    typedef typename GET_PROP_TYPE(TypeTag, PressureModel) Implementation;
    typedef FVPressure2P2C<TypeTag> ParentType;
    typedef FVPressure<TypeTag> BaseType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };
    enum
    {
        rhs = BaseType::rhs, matrix = BaseType::matrix,
    };

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dim+1> TransmissivityMatrix;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, GET_PROP_VALUE(TypeTag, NumPhases)> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    // the typenames used for the stiffness matrix and solution vector
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) RHSVector;

protected:
    Problem& problem()
    {
        return this->problem_;
    }
    const Problem& problem() const
    {
        return this->problem_;
    }

public:
    //initializes the matrix to store the system of equations
    void initializeMatrix();
    //function which assembles the system of equations to be solved
    void assemble(bool first);

    void getMpfaFlux(const IntersectionIterator&, const CellData&);

    // mpfa transmissibilities
    int computeTransmissibilities(const IntersectionIterator&,
                                    IntersectionIterator&,
                                    TransmissivityMatrix&,
                                    TransmissivityMatrix&,
                                    GlobalPosition&,
                                    int&);

    //! Adapt primary variables vector after adapting the grid
    void adaptPressure()
    {
        int gridSize = problem().gridView().size(0);
        this->pressure().resize(gridSize);

        for(int i=0; i< gridSize; i++)
        {
            this->pressure()[i]
                  = problem().variables().cellData(i).pressure(this->pressureType);
        }
    }

    //! Constructs a FVPressure2P2CAdaptive object
    /**
     * \param problem a problem class object
     */
    FVPressure2P2CAdaptive(Problem& problem) : FVPressure2P2C<TypeTag>(problem),
            problem_(problem), pressureModelAdaptive2p_(0)
    {
        enableVolumeIntegral = GET_PARAM_FROM_GROUP(TypeTag,bool, Impet, EnableVolumeIntegral);
        enableMPFA= GET_PARAM_FROM_GROUP(TypeTag,bool, GridAdapt, EnableMultiPointFluxApproximation);

        enableSecondHalfEdge = GET_PARAM_FROM_GROUP(TypeTag,bool, GridAdapt, EnableSecondHalfEdge);
        // prepare rotation matrix: R_ is initialized as zero matrix
        R_=0.;
        // evaluate matrix R
        if (dim == 2)
        {
            R_[0][1] = 1;
            R_[1][0] = -1;
        }
    }

private:
    Problem& problem_;
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

protected:
    // mpfa transmissibilities from mpfal2pfa
    int transmissibilityAdapter_(const IntersectionIterator&,
            const IntersectionIterator&,
            IntersectionIterator&,
            GlobalPosition&,
            TransmissivityMatrix&);

    //! Matrix for vector rotation used in mpfa
    DimMatrix R_;
    bool enableVolumeIntegral; //!< Enables the volume integral of the pressure equation
    bool enableMPFA; //!< Enables mpfa method to calculate the fluxes near hanging nodes
    bool enableSecondHalfEdge; //!< If possible, 2 interaction volumes are used for the mpfa method near hanging nodes
    //! The 2p Mpfa pressure module, that is only used for the calulation of transmissibility of the second interaction volumes
    FVMPFAL2PFABoundPressure2PAdaptive<TypeTag>* pressureModelAdaptive2p_;
};


//! initializes the matrix to store the system of equations
/*! In comparison with the Tpfa method, an mpfa uses a larger flux stencil, hence more
 * matrix entries are required if not only the unique interaction region on the hanging
 * nodes are considered. The method checks weather the additonally regarded cells through
 * mpfa are already "normal" neighbors for fluxes through other interfaces, or if they
 * need to be added.
 */
template<class TypeTag>
void FVPressure2P2CAdaptive<TypeTag>::initializeMatrix()
{
    int gridSize_ = problem().gridView().size(0);
    // update RHS vector, matrix
    this->A_.setSize (gridSize_,gridSize_); //
    this->f_.resize(gridSize_);

    // determine matrix row sizes
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        // initialize row size
        int rowSize = 1;

        // set perimeter to zero
        cellDataI.perimeter() = 0;

        // prepare storage for all found 3rd cells
        std::vector<int> foundAdditionals;

        int numberOfIntersections = 0;
        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            cellDataI.perimeter()
                    += isIt->geometry().volume();
            numberOfIntersections++;
            if (isIt->neighbor())
            {
                rowSize++;

                // if mpfa is used, more entries might be needed if both halfedges are regarded
                if (enableMPFA && enableSecondHalfEdge && isIt->outside()->level()!=eIt.level())
                {
                    GlobalPosition globalPos3(0.);
                    int globalIdx3=-1;
                    TransmissivityMatrix T(0.);
                    IntersectionIterator additionalIsIt = isIt;
                    TransmissivityMatrix additionalT(0.);
                    // compute Transmissibilities: also examines subcontrolvolume information
                    int halfedgesStored
                        = problem().variables().getMpfaData(*isIt, additionalIsIt, T, additionalT, globalPos3, globalIdx3);
                    if (halfedgesStored == 0)
                        halfedgesStored = problem().pressureModel().computeTransmissibilities(isIt,additionalIsIt, T,additionalT, globalPos3, globalIdx3 );

                    if(halfedgesStored == 2)
                    {
                        bool increaseRowSize = true;
                        //check if additional cell is ordinary neighbor of eIt
                        for (IntersectionIterator isIt2 = problem().gridView().template ibegin(*eIt); isIt2 != isItEnd; ++isIt2)
                        {
                            if(!isIt2->neighbor())
                                continue;
                            if(additionalIsIt->outside() == isIt2->outside() )
                                increaseRowSize = false;
                        }
                        //also check if additional cell was already used for another interaction triangle
                        for(int i =0; i<foundAdditionals.size(); i++)
                            if(foundAdditionals[i] == problem().variables().index(*additionalIsIt->outside()))
                                    increaseRowSize = false;

                        if (increaseRowSize)
                        {
                            rowSize++;
                            foundAdditionals.push_back(problem().variables().index(*additionalIsIt->outside()));
                        }
                    }
                }
            }
        }

        cellDataI.fluxData().resize(numberOfIntersections);
        this->A_.setrowsize(globalIdxI, rowSize);
    }
    this->A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem().variables().index(*eIt);

        // add diagonal index
        this->A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem().variables().index(*outside);

                // add off diagonal index
                this->A_.addindex(globalIdxI, globalIdxJ);

                // if mpfa is used, more entries might be needed if both halfedges are regarded
                if (enableMPFA && enableSecondHalfEdge && isIt->outside()->level()!=eIt.level())
                {
                    GlobalPosition globalPos3(0.);
                    int globalIdx3=-1;
                    TransmissivityMatrix T(0.);
                    IntersectionIterator additionalIsIt = isIt;
                    TransmissivityMatrix additionalT(0.);
                    // compute Transmissibilities: also examines subcontrolvolume information
                    int halfedgesStored
                        = problem().variables().getMpfaData(*isIt, additionalIsIt, T, additionalT, globalPos3, globalIdx3);
                    if (halfedgesStored == 0)
                        halfedgesStored = problem().pressureModel().computeTransmissibilities(isIt,additionalIsIt, T,additionalT, globalPos3, globalIdx3 );
                    // add off diagonal index if 2 half-edges regarded
                    if(halfedgesStored == 2)
                        this->A_.addindex(globalIdxI, problem().variables().index(*additionalIsIt->outside()));
                }
            }
    }
    this->A_.endindices();

    return;
}

//! function which assembles the system of equations to be solved
/*! This function assembles the Matrix and the RHS vectors to solve for
 * a pressure field with an Finite-Volume Discretization in an implicit
 * fashion. Compared to the method in FVPressure, this implementation
 * calculates fluxes near hanging nodes with the mpfa method using the method
 * getMpfaFlux(). Matrix and Right-hand-side entries are done therein.
 *
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2CAdaptive<TypeTag>::assemble(bool first)
{
    if(first)
    {
        BaseType::assemble(true);
        return;
    }

    initializeMatrix();
    // initialization: set matrix A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    ElementIterator eItEnd = problem().gridView().template end<0>();
    for (ElementIterator eIt = problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell information
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        Dune::FieldVector<Scalar, 2> entries(0.);

        /*****  source term ***********/
        problem().pressureModel().getSource(entries, *eIt, cellDataI, first);
        this->f_[globalIdxI] += entries[rhs];

        /*****  flux term ***********/
        // iterate over all faces of the cell
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            /************* handle interior face *****************/
            if (isIt->neighbor())
            {
                ElementPointer elementNeighbor = isIt->outside();

                int globalIdxJ = problem().variables().index(*elementNeighbor);

                //check for hanging nodes
                //take a hanging node never from the element with smaller level!
                bool haveSameLevel = (eIt->level() == elementNeighbor->level());
                //calculate only from one side, but add matrix entries for both sides
                if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce) && (globalIdxI > globalIdxJ) && haveSameLevel)
                    continue;

                entries = 0;
                //check for hanging nodes
                if(!haveSameLevel && enableMPFA)
                {
                    problem().pressureModel().getMpfaFlux(isIt, cellDataI);
                }
                else
                {
                    problem().pressureModel().getFlux(entries, *isIt, cellDataI, first);

                    //set right hand side
                    this->f_[globalIdxI] -= entries[rhs];

                    // set diagonal entry
                    this->A_[globalIdxI][globalIdxI] += entries[matrix];

                        // set off-diagonal entry
                    this->A_[globalIdxI][globalIdxJ] -= entries[matrix];

                    if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce))
                    {
                        this->f_[globalIdxJ] += entries[rhs];
                        this->A_[globalIdxJ][globalIdxJ] += entries[matrix];
                        this->A_[globalIdxJ][globalIdxI] -= entries[matrix];
                    }
                }

            } // end neighbor

            /************* boundary face ************************/
            else
            {
                entries = 0;
                problem().pressureModel().getFluxOnBoundary(entries, *isIt, cellDataI, first);

                //set right hand side
                this->f_[globalIdxI] += entries[rhs];
                // set diagonal entry
                this->A_[globalIdxI][globalIdxI] += entries[matrix];
            }
        } //end interfaces loop
//        printmatrix(std::cout, this->A_, "global stiffness matrix", "row", 11, 3);

        /*****  storage term ***********/
        entries = 0;
        problem().pressureModel().getStorage(entries, *eIt, cellDataI, first);
        this->f_[globalIdxI] += entries[rhs];
//         set diagonal entry
        this->A_[globalIdxI][globalIdxI] += entries[matrix];
    } // end grid traversal
//    printmatrix(std::cout, this->A_, "global stiffness matrix after assempling", "row", 11,3);
//    printvector(std::cout, this->f_, "right hand side", "row", 10);
    return;
}

//! Compute flux over an irregular interface using a \a mpfa method
/** A mpfa l-method is applied to calculate fluxes near hanging nodes, using:
 * \f[
      - \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
        \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha} \frac{\partial v_{t}}{\partial C^{\kappa}}
      + \frac{ V_i}{U_i} \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
       \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
          \sum_{\kappa} X^{\kappa}_{\alpha} \frac{\frac{\partial v_{t,j}}{\partial C^{\kappa}_j}-\frac{\partial v_{t,i}}{\partial C^{\kappa}_i}}{\Delta x}
    \f]
 *
 * We provide two options: Calculating the flux expressed by twice the flux
 * through the one unique interaction region on the hanging node if one
 * halfedge is stored (eg on boundaries). Or using the second interaction
 * region covering neighboring cells.
 * Due to that, the matrix and rhs entries are filled up within this function.
 * \param intersectionIterator Iterator of the intersection between cell I and J
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FVPressure2P2CAdaptive<TypeTag>::getMpfaFlux(const IntersectionIterator& intersectionIterator,
                                                    const CellData& cellDataI)
{
    // acess Cell I
    ElementPointer elementPointerI = intersectionIterator->inside();
    int globalIdxI = problem().variables().index(*elementPointerI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementPointerI->geometry().volume();
    Scalar perimeter = cellDataI.perimeter();

    const GlobalPosition& gravity_ = problem().gravity();

    // get absolute permeability
    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(*elementPointerI));

    // access neighbor
    ElementPointer neighborPointer = intersectionIterator->outside();
    int globalIdxJ = problem().variables().index(*neighborPointer);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    DimMatrix permeabilityJ
        = problem().spatialParams().intrinsicPermeability(*neighborPointer);

    // compute vectorized permeabilities
    DimMatrix meanPermeability(0);
    Dumux::harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

    Dune::FieldVector<Scalar, dim> permeability(0);
    meanPermeability.mv(unitDistVec, permeability);

    // get average density for gravity flux
    Scalar rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
    Scalar rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));

    // reset potential gradients
    Scalar potentialW = 0;
    Scalar potentialNW = 0;

    // determine volume derivatives in neighbor
    if (!cellDataJ.hasVolumeDerivatives())
        this->volumeDerivatives(globalPosNeighbor, *neighborPointer);

    Scalar dv_dC1 = (cellDataJ.dv(wPhaseIdx)
                + cellDataI.dv(wPhaseIdx)) / 2; // dV/dm1= dv/dC^1
    Scalar dv_dC2 = (cellDataJ.dv(nPhaseIdx)
                + cellDataI.dv(nPhaseIdx)) / 2;

    Scalar graddv_dC1 = (cellDataJ.dv(wPhaseIdx)
                            - cellDataI.dv(wPhaseIdx)) / dist;
    Scalar graddv_dC2 = (cellDataJ.dv(nPhaseIdx)
                            - cellDataI.dv(nPhaseIdx)) / dist;


//                    potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                    potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
//
//                    densityW = (potentialW > 0.) ? densityWI : densityWJ;
//                    densityNW = (potentialNW > 0.) ? densityNWI : densityNWJ;
//
//                    densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                    densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
    //jochen: central weighting for gravity term
    Scalar densityW = rhoMeanW;
    Scalar densityNW = rhoMeanNW;

        // Prepare MPFA
        /** get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int globalIdx3=-1;
        TransmissivityMatrix T(0.);
            // prepare second half-edge
            IntersectionIterator additionalIsIt = intersectionIterator;
            TransmissivityMatrix additionalT(0.);

        int halfedgesStored = problem().variables().getMpfaData(*intersectionIterator,
                                            additionalIsIt, T, additionalT,
                                            globalPos3, globalIdx3 );
	    if (halfedgesStored == 0)
            halfedgesStored = problem().pressureModel().computeTransmissibilities(intersectionIterator,
                                                                    additionalIsIt, T, additionalT,
                                                                    globalPos3, globalIdx3 );

        if(!halfedgesStored)
            Dune::dgrave << "something went wrong getting mpfa data on cell " << globalIdxI << std::endl;

        // shortcurts mpfa case
        CellData& cellData3 = problem().variables().cellData(globalIdx3);
        Scalar temp1 = globalPos * gravity_;
        Scalar temp2 = globalPosNeighbor * gravity_;
        Scalar temp3 = globalPos3 * gravity_;

        potentialW = (cellDataI.pressure(wPhaseIdx)-temp1*densityW) * T[2]
                     + (cellDataJ.pressure(wPhaseIdx)-temp2*densityW) * T[0]
                     + (cellData3.pressure(wPhaseIdx)-temp3*densityW) * T[1];
        potentialNW = (cellDataI.pressure(nPhaseIdx)-temp1*densityNW) * T[2]
                     + (cellDataJ.pressure(nPhaseIdx)-temp2*densityNW) * T[0]
                     + (cellData3.pressure(nPhaseIdx)-temp3*densityNW) * T[1];
        // regard second half edge, if there is one
        if(halfedgesStored == 2)
        {
            int AdditionalIdx = problem().variables().index(*(additionalIsIt->outside()));
            CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);
            potentialW += (cellDataI.pressure(wPhaseIdx)-temp1*densityW) * additionalT[2]
                            + (cellDataJ.pressure(wPhaseIdx)-temp2*densityW) * additionalT[0]
                            + (cellDataAdditional.pressure(wPhaseIdx)
                                -(additionalIsIt->outside()->geometry().center()*gravity_)
                              *densityW) * additionalT[1];
            potentialNW += (cellDataI.pressure(nPhaseIdx)-temp1*densityNW) * additionalT[2]
                            + (cellDataJ.pressure(nPhaseIdx)-temp2*densityNW) * additionalT[0]
                            + (cellDataAdditional.pressure(nPhaseIdx)
                                -(additionalIsIt->outside()->geometry().center()*gravity_)
                              *densityNW) * additionalT[1];
        }


        // initialize convenience shortcuts
        Scalar lambdaW(0.), lambdaN(0.);
        Scalar dV_w(0.), dV_n(0.);        // dV_a = \sum_k \rho_a * dv/dC^k * X^k_a
        Scalar gV_w(0.), gV_n(0.);        // multipaper eq(3.3) line 3 analogon dV_w


        //do the upwinding of the mobility depending on the phase potentials
        if (potentialW >= 0.)
        {
            dV_w = (dv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
            lambdaW = cellDataI.mobility(wPhaseIdx);
            gV_w = (graddv_dC1 * cellDataI.massFraction(wPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataI.massFraction(wPhaseIdx, nCompIdx));
            dV_w *= cellDataI.density(wPhaseIdx);
            gV_w *= cellDataI.density(wPhaseIdx);
        }
        else
        {
            dV_w = (dv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
            lambdaW = cellDataJ.mobility(wPhaseIdx);
            gV_w = (graddv_dC1 * cellDataJ.massFraction(wPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataJ.massFraction(wPhaseIdx, nCompIdx));
            dV_w *= cellDataJ.density(wPhaseIdx);
            gV_w *= cellDataJ.density(wPhaseIdx);
        }
        if (potentialNW >= 0.)
        {
            dV_n = (dv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
            lambdaN = cellDataI.mobility(nPhaseIdx);
            gV_n = (graddv_dC1 * cellDataI.massFraction(nPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataI.massFraction(nPhaseIdx, nCompIdx));
            dV_n *= cellDataI.density(nPhaseIdx);
            gV_n *= cellDataI.density(nPhaseIdx);
        }
        else
        {
            dV_n = (dv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                    + dv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
            lambdaN = cellDataJ.mobility(nPhaseIdx);
            gV_n = (graddv_dC1 * cellDataJ.massFraction(nPhaseIdx, wCompIdx)
                    + graddv_dC2 * cellDataJ.massFraction(nPhaseIdx, nCompIdx));
            dV_n *= cellDataJ.density(nPhaseIdx);
            gV_n *= cellDataJ.density(nPhaseIdx);
        }

    /** compute matrix entry: advective fluxes */
    /* extend T with other matrix entries and assemble to A_    */
    this->A_[globalIdxI][globalIdxJ] += (lambdaW * dV_w + lambdaN * dV_n) * T[0];
    this->A_[globalIdxI][globalIdx3] += (lambdaW * dV_w + lambdaN * dV_n) * T[1];
    this->A_[globalIdxI][globalIdxI] += (lambdaW * dV_w + lambdaN * dV_n) * T[2];

    // add gravity to RHS vector
    this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n) * temp2 * T[0];
    this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n) * temp3 * T[1];
    this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n) * temp1 * T[2];

    // weithing accounts for the fraction of the subcontrol volume
    Scalar weightingFactor = volume / perimeter;    // transforms flux through area A -> V * A/perimeter
    if(enableVolumeIntegral) // switch off volume integral for mpfa case
    {
        // correct for area integral
        this->A_[globalIdxI][globalIdxJ] -= weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*T[0];
        this->A_[globalIdxI][globalIdx3] -= weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*T[1];
        this->A_[globalIdxI][globalIdxI] -= weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*T[2];

        // add gravity to RHS vector
        this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n) * temp2 * T[0];
        this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n) * temp3 * T[1];
        this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n) * temp1 * T[2];
    }

    // capillary pressure flux
    Scalar pcGradient = cellDataI.capillaryPressure() * T[2]
                       + cellDataJ.capillaryPressure() * T[0]
                       + cellData3.capillaryPressure() * T[1];

    if (this->pressureType == pw)
        pcGradient *= + lambdaN * dV_n - enableVolumeIntegral * weightingFactor * lambdaN * gV_n;
    else if (this->pressureType == pn)
        pcGradient *= - lambdaW * dV_w + enableVolumeIntegral * weightingFactor * lambdaW * gV_w;

    this->f_[globalIdxI] += pcGradient;

    // include second half-edge
    if(halfedgesStored == 2)
    {
        int AdditionalIdx = problem().variables().index(*(additionalIsIt->outside()));
        CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);

        /* extend T with other matrix entries and assemble to A_    */
        this->A_[globalIdxI][globalIdxJ] += (lambdaW * dV_w + lambdaN * dV_n) * additionalT[0];
        this->A_[globalIdxI][AdditionalIdx] +=  (lambdaW * dV_w + lambdaN * dV_n) * additionalT[1];
        this->A_[globalIdxI][globalIdxI] +=  (lambdaW * dV_w + lambdaN * dV_n) * additionalT[2];

        // add gravity to RHS vector
        this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n) * temp2 * additionalT[0];
        this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n)
                * (additionalIsIt->outside()->geometry().center()*gravity_) * additionalT[1];
        this->f_[globalIdxI] += (densityW * lambdaW * dV_w + densityNW * lambdaN * dV_n) * temp1 * additionalT[2];

        if(enableVolumeIntegral) // switch off volume integral for mpfa case
        {
            // correct for area integral
            this->A_[globalIdxI][globalIdxJ] -=  weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*additionalT[0];
            this->A_[globalIdxI][AdditionalIdx] -= weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*additionalT[1];
            this->A_[globalIdxI][globalIdxI] -= weightingFactor * (lambdaW * gV_w + lambdaN * gV_n)*additionalT[2];

            // add gravity to RHS vector
            this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n) * temp2 * additionalT[0];
            this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n)
                    * (additionalIsIt->outside()->geometry().center()*gravity_) * additionalT[1];
            this->f_[globalIdxI] -= weightingFactor * (densityW * lambdaW * gV_w + densityNW * lambdaN * gV_n) * temp1 * additionalT[2];
        }

        // capillary pressure flux
        pcGradient = cellDataI.capillaryPressure() * additionalT[2]
                           + cellDataJ.capillaryPressure() * additionalT[0]
                           + cellDataAdditional.capillaryPressure() * additionalT[1];
        if (this->pressureType == pw)
            pcGradient *= + lambdaN * dV_n - enableVolumeIntegral * weightingFactor * lambdaN * gV_n;
        else if (this->pressureType == pn)
            pcGradient *= - lambdaW * dV_w + enableVolumeIntegral * weightingFactor * lambdaW * gV_w;

        this->f_[globalIdxI] += pcGradient;
    }
}

//! Computes the transmissibility coefficients for the MPFA-l method
/*!  Indices used in a interaction volume of the MPFA-l method
 * \verbatim
               ___________________________________________________
               | nux: cell geometry     | nx: face normal        |
               |                        | =integrationOuterNormal|
               |                        |                        |
               |  nextIt________________|                        |
               |                        |    nu4    3            |
               |                        |      \;.´ |            |
               |                        |     .´ <--|nu3         |
               |   elementnumber        |  .´ ^n1   |            |
               |            1...........x_____|_____|____________|
               |             `.    |nu5 !`.   _nu7  |            |
               |                `.`'´   !   `./` <--|            |
       lastIt  |                  `._   !     `. nu2|            |
            \  |                 nu6/`. !nu1^   ` . |            |
             \ |                       `!---|------`2            |
              `|                        |                        |
               |\ isIt__________________|->n12                   |
               | \                      |                        |
               |_face14_________________|________________________|


     _
     /` is a normal directing upwards to the right,
     while a normal directing downwards to the right looks like \;
\endverbatim
*
* \param isIt Iterator to the current intersection
* \param additionalIntersectionIt Iterator to the additional interface included in the second interaction region
* \param T Transmissitivity matrix of the first unique interaction region
* \param additionalT Transmissitivity matrix of the second non-unique interaction region
* \param globalPos3 Unique interaction region: Position of the 3rd cell, the other neighbor on the hanging node
* \param globalIdx3 Unique interaction region: Index of the 3rd cell, the other neighbor on the hanging node
*/
template <class TypeTag>
int FVPressure2P2CAdaptive<TypeTag>::computeTransmissibilities(const IntersectionIterator& isIt,
        IntersectionIterator& additionalIntersectionIt,
        TransmissivityMatrix& T,
        TransmissivityMatrix& additionalT,
        GlobalPosition& globalPos3,
        int& globalIdx3)
{
    // get geometry information of cellI = cell1, cellJ = cell2
    ElementPointer eIt = isIt->inside();
    ElementPointer neighborPointer = isIt->outside();
    GlobalPosition globalPos1 = eIt->geometry().center();
    GlobalPosition globalPos2 = neighborPointer->geometry().center();
    DimMatrix K1(problem().spatialParams().intrinsicPermeability(*eIt));
    DimMatrix K2(problem().spatialParams().intrinsicPermeability(*neighborPointer));
//    int debugIdx =  problem().variables().index(*isIt->outside());

    /* 1. get geometrical information of interaction triangle   */
    // geometry and Data of face IJ in nomenclature of mpfa
    GlobalPosition globalPosFace12 = isIt->geometry().center();
    GlobalPosition integrationOuterNormaln12 = isIt->centerUnitOuterNormal();
    integrationOuterNormaln12 *= isIt->geometry().volume() / 2.0; // TODO: 2.0 only in 2D


    // nextIs points to next intersection
    IntersectionIterator nextIs = isIt;
    ++nextIs;
    if (nextIs== problem().gridView().template iend(*eIt))
        nextIs = problem().gridView().template ibegin(*eIt);

    // get last intersection : --intersection does not exist
    // paceingIt loops one IS bevore prevIs
    IntersectionIterator prevIs = problem().gridView().template ibegin(*eIt);
    IntersectionIterator paceingIt = prevIs;
    for (++paceingIt; paceingIt != problem().gridView().template iend(*eIt); ++paceingIt)
    {
//        if((paceingIt->neighbor()))
//            debugIdx = problem().variables().index(*paceingIt->outside());
//        if((prevIs->neighbor()))
//            debugIdx = problem().variables().index(*prevIs->outside());

        if (!paceingIt->neighbor())  // continue if no neighbor found
            ++prevIs;   // we investigate next paceingIt -> prevIs is also increased
        else if (paceingIt->outside() == isIt->outside())  // we already found prevIs
                break;
        else if (paceingIt == problem().gridView().template iend(*eIt))
                prevIs = paceingIt; // this could only happen if isIt is begin, so prevIs has to be last.
        else
            ++prevIs;   // we investigate next paceingIt -> prevIs is also increased
    }

//        GlobalPosition prevIsCoords = prevIs->geometry().center();
//        GlobalPosition nextIsCoords = nextIs->geometry().center();

    // search for face13, face23
    IntersectionIterator face13 = isIt; // as long as isIt13=isIt, it is still not found!
    IntersectionIterator face23 = isIt; // as long as face23=isIt, it is still not found!
    // store other intersection for the other interaction region for the other half-edge

    IntersectionIterator isIt23End = problem().gridView().iend(*neighborPointer);
    for (IntersectionIterator isIt23 = problem().gridView().ibegin(*neighborPointer); isIt23 != isIt23End; ++isIt23)
    {
//        if((isIt23->neighbor()))
//            debugIdx = problem().variables().index(*isIt23->outside());
        // stop search if found
        if( (face13->outside() != isIt->outside()))
            break;

        if(!(isIt23->neighbor()))
            continue;

        // either prevIs or nextIs is face13, it is that with common interface with cell2
        // investigate if prevIs points to cell 3
        if (prevIs->neighbor())
        {
//            if((prevIs->neighbor()))
//                debugIdx = problem().variables().index(*prevIs->outside());

            if (prevIs->outside() == isIt23->outside())
            {
                face23 = isIt23;
                face13 = prevIs;
                additionalIntersectionIt = nextIs;
            }
        }
        // investigate if nextIs points to cell 3
        if (nextIs->neighbor())
        {
//            if((nextIs->neighbor()))
//                debugIdx = problem().variables().index(*nextIs->outside());

            if (nextIs->outside() == isIt23->outside())
            {
                face23 = isIt23;
                face13 = nextIs;
                additionalIntersectionIt = prevIs;
            }
        }
    }
    if(face13->outside() == isIt->outside()) //means isIt13 not found yet
        Dune::dgrave << "is 13 not found!!!" << std::endl;

    // get information of cell3
    globalPos3 = face13->outside()->geometry().center();
    globalIdx3 = problem().variables().index(*(face13->outside()));
    // get absolute permeability of neighbor cell 3
    DimMatrix K3(problem().spatialParams().intrinsicPermeability(*(face13->outside())));


    // get the intersection node /bar^{x_3} between 'isIt' and 'isIt13', denoted as 'corner123'
    // initialization of corner123
    GlobalPosition corner123(0.);

    // get the global coordinate of corner123
    for (int i = 0; i < isIt->geometry().corners(); ++i)
    {
        for (int j = 0; j < face13->geometry().corners(); ++j)
        {
            if (face13->geometry().corner(j) == isIt->geometry().corner(i))
            {
                corner123 = face13->geometry().corner(j);

                // stop outer (i) and inner (j) loop
                i = isIt->geometry().corners();
                break;
            }
        }
    }

    /** 3. Calculate omega, chi for matrices  **/
    // center of face in global coordinates, i.e., the midpoint of edge 'isIt24'
    GlobalPosition globalPosFace23 = face23->geometry().center();

    // get face volume
    Scalar face23vol = face23->geometry().volume();

    // get outer normal vector scaled with half volume of face 'isIt24'
    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln23 = face23->centerUnitOuterNormal();
    integrationOuterNormaln23 *= face23vol / 2.0; //TODO: 2.0 only for 2D

    // compute normal vectors nu1-nu7 in triangle R for first half edge
    GlobalPosition nu1(0);
    R_.umv(globalPosFace12-globalPos2 ,nu1);    //globalPos2 = globalPosNeighbor

    GlobalPosition nu2(0);
    R_.umv(globalPos2-globalPosFace23, nu2);

    GlobalPosition nu3(0);
    R_.umv(globalPosFace23-globalPos3, nu3);

    GlobalPosition nu4(0);
    R_.umv(globalPos3-corner123, nu4);

    GlobalPosition nu5(0);
    R_.umv(corner123-globalPos1, nu5);  //globalPos1 = globalPos

    GlobalPosition nu6(0);
    R_.umv(globalPos1-globalPosFace12, nu6);

    GlobalPosition nu7(0);
    R_.umv(corner123-globalPos2, nu7);

    // compute T, i.e., the area of quadrilateral made by normal vectors 'nu'
    GlobalPosition Rnu2(0);
    R_.umv(nu2, Rnu2);
    Scalar T1 = nu1 * Rnu2;

    GlobalPosition Rnu4(0);
    R_.umv(nu4, Rnu4);
    Scalar T2 = nu3 * Rnu4;

    GlobalPosition Rnu6(0);
    R_.umv(nu6, Rnu6);
    Scalar T3 = nu5 * Rnu6;

    // compute components needed for flux calculation, denoted as 'omega' and 'chi'
    GlobalPosition K2nu1(0);
    K2.umv(nu1, K2nu1); // permeabilityJ = K2 = perm of cell 2 around x_1
    GlobalPosition K2nu2(0);
    K2.umv(nu2, K2nu2);
    GlobalPosition K3nu3(0);
    K3.umv(nu3, K3nu3);
    GlobalPosition K3nu4(0);
    K3.umv(nu4, K3nu4);
    GlobalPosition K1nu5(0);
    K1.umv(nu5, K1nu5); // permeabilityI = K1 = perm of cell 1 around x_3
    GlobalPosition K1nu6(0);
    K1.umv(nu6, K1nu6);
    GlobalPosition Rnu1(0);
    R_.umv(nu1, Rnu1);

    double omega111 = /*lambda2 */ (integrationOuterNormaln23 * K2nu1)/T1;
    double omega112 = /*lambda2 */ (integrationOuterNormaln23 * K2nu2)/T1;
    double omega211 = /*lambda2 */ (integrationOuterNormaln12 * K2nu1)/T1;
    double omega212 = /*lambda2 */ (integrationOuterNormaln12 * K2nu2)/T1;
    double omega123 = /*lambda3 */ (integrationOuterNormaln23 * K3nu3)/T2;
    double omega124 = /*lambda3 */ (integrationOuterNormaln23 * K3nu4)/T2;
    double omega235 = /*lambda1 */ (integrationOuterNormaln12 * K1nu5)/T3;
    double omega236 = /*lambda1 */ (integrationOuterNormaln12 * K1nu6)/T3;
    double chi711 = (nu7 * Rnu1)/T1;
    double chi712 = (nu7 * Rnu2)/T1;

    /** 4. Calculate A, B, C, D and solve for T **/
    // compute transmissibility matrix T = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega111;
    C[0][1] = -omega112;
    C[1][0] = -omega211;
    C[1][1] = -omega212;

    D[0][0] = omega111 + omega112;
    D[1][0] = omega211 + omega212;

    A[0][0] = omega111 - omega124 - omega123*chi711;
    A[0][1] = omega112 - omega123*chi712;
    A[1][0] = omega211 - omega236*chi711;
    A[1][1] = omega212 - omega235 - omega236*chi712;

    B[0][0] = omega111 + omega112 + omega123*(1.0 - chi711 - chi712);
    B[0][1] = -omega123 - omega124;
    B[1][0] = omega211 + omega212 + omega236*(1.0 - chi711 - chi712);
    B[1][2] = -omega235 - omega236;

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));
    T = D[1];
    if(!enableSecondHalfEdge )//or fabs(isIt->centerUnitOuterNormal()[0])<0.5) // [0]<0.5 => switch off vertical 2hes
    {
        T *= 2;
        // set your map entry
        problem().variables().storeMpfaData(*isIt, T, globalPos3, globalIdx3);
        return 1; // indicates that only 1 halfedge was regarded
    }
    else
    {
        // derive additional T for second half edge

        // get the intersection node /bar^{x_3} between 'isIt' and 'aditionalIntersection', denoted as 'corner1245'
        // initialization of corner1245
        GlobalPosition corner1245(0.);

        // get the global coordinate of corner1245, which is not connected with face13
        IntersectionIterator tempIntersection = face13;
        bool corner1245found = false;
        // ensure iterator increases over local end
        if (tempIntersection== problem().gridView().template iend(*eIt))
            tempIntersection = problem().gridView().template ibegin(*eIt);
        while (!corner1245found)
        {
            ++tempIntersection;

            // ensure iterator increases over local end
            if (tempIntersection== problem().gridView().template iend(*eIt))
                tempIntersection = problem().gridView().template ibegin(*eIt);
            // enshure we do not arrive at isIt
            if (tempIntersection == isIt)
                continue;

            // loop over both corners of is
            for (int i = 0; i < isIt->geometry().corners(); ++i)
            {
                // test if a corner of additionalIntersectionIt also lies on is
                for (int j = 0; j < tempIntersection->geometry().corners(); ++j)
                {
                    if (tempIntersection->geometry().corner(j) == isIt->geometry().corner(i))
                    {
                        corner1245 = tempIntersection->geometry().corner(j);
                        additionalIntersectionIt = tempIntersection;

                        // stop outer (i) and inner (j) loop
                        i = isIt->geometry().corners();
                        // stop also Intersection loop
                        corner1245found = true;
                        break;
                    }
                }
            }
        }
        if (!additionalIntersectionIt->neighbor())// or (additionalIntersectionIt->outside().level() == isIt->inside().level()))
        {
            T *= 2;
            // set your map entry
            problem().variables().storeMpfaData(*isIt, T, globalPos3, globalIdx3);
            return 1;
        }

        // use Markus mpfa-l implementation
        this->transmissibilityAdapter_(isIt, face23, additionalIntersectionIt,
                                        corner1245, additionalT);

        // store transmissivity data for second half edge
        problem().variables().storeMpfaData(*isIt, additionalIntersectionIt, T,additionalT, globalPos3, globalIdx3);
        // return that information about two interaction volumes have been stored
        return 2;
    }
}

//! An adapter to use the traditional 2p implementation for the second interaction region
/*!
 * The second interaction region consists of 4 cells, hence the traditional non-adaptive
 * implementation to calculate the transmissibility coefficients is used. This, however,
 * uses its own naming conventions and local Indices used in the work of I. Aavatsmark.
 *
 * Indices used in a interaction volume of the MPFA-l method
 * \verbatim
                 |                        |                        |
                 |            4-----------3-----------3            |
                 |            | --> nu43  |  nu34 <-- |            |
                 |            | |nu41    1|--> n43   ||nu32        |
                 |            | v   ^     |0     ^   v|            |
                 |____________4__0__|n14__|__n23_|_1__2____________|
                 |            |    1    0 |     0     |            |
                 |            | ^         |1   nu23 ^ |            |
                 |            | |nu14    0|--> n12  | |            |
                 |            | -->nu12   |   nu21<-- |            |
                 |        J = 1-----------1-----------2 = I        |
                 |                        |                        |
    \endverbatim.
* \param isIt Iterator to the current intersection
* \param face23_2p2cnaming Iterator to the second interface included in the unique interaction region = face23
* \param additionalIntersectionIt Iterator to the intersection included in the second interaction region
* \param corner1234 Center point of non-unique interaction region.
* \param additionalT Transmissibility Matrix of non unique interaction region.
*/
template<class TypeTag>
int FVPressure2P2CAdaptive<TypeTag>::transmissibilityAdapter_(const IntersectionIterator& isIt,
                                const IntersectionIterator& face23_2p2cnaming,
                                IntersectionIterator& additionalIntersectionIt,
                                GlobalPosition& corner1234,
                                TransmissivityMatrix& additionalT)
{
    /****** find all 4 faces *********/
    // store additionalIntersection, which is face23 in 2p mpfa-l naming scheme
    IntersectionIterator isIt23 = additionalIntersectionIt;
    IntersectionIterator isIt14 = isIt; //has to be initialized with something

    // search for 4th intersection that connects to corner1245, which is needed
    // for markus implementation of the mpfa

    // reuse corner1245found bool
    bool face14found = false;
    // repeat search procedure with isIt and face23 (in 2p2c naming)
    IntersectionIterator tempIntersection = face23_2p2cnaming;

    while (!face14found)
    {
        ++tempIntersection;

        // ensure iterator increases over local end of neighbor J
        if (tempIntersection== problem().gridView().template iend(*isIt->outside()))
            tempIntersection = problem().gridView().template ibegin(*isIt->outside());

        if(!tempIntersection->neighbor())
            continue;

        // enshure we hace not arrived at isIt (but seen from other side)
        if (tempIntersection->outside() == isIt->inside())
            continue; // restart loop to recheck after next increase for local end

        // loop over both corners of is
        for (int i = 0; i < isIt->geometry().corners(); ++i)
        {
            // test if a corner of additionalIntersectionIt also lies on is
            for (int j = 0; j < tempIntersection->geometry().corners(); ++j)
            {
                if (tempIntersection->geometry().corner(j) == isIt->geometry().corner(i))
                {
//                        // test if this tempIntersection is better than additionalIntersectionIt
//                        if (tempIntersection->outside().level() > additionalIntersectionIt->outside().level())
//                        {
//                            // select this as element for second half edge
//                            additionalIntersectionIt = tempIntersection;
//                            corner1245 = additionalIntersectionIt->geometry().corner(j);
//                        }
                    isIt14 = tempIntersection;
                    // stop outer (i) and inner (j) loop
                    i = isIt->geometry().corners();
                    // stop also Intersection loop
                    face14found = true;
                    break;
                }
            }
        }
    }
    /**** end find 4 faces **/


    if(!pressureModelAdaptive2p_)
        pressureModelAdaptive2p_= new FVMPFAL2PFABoundPressure2PAdaptive<TypeTag>(problem()) ;
    // create Interaction Volume object
    Dumux::FVMPFALInteractionVolume<TypeTag> interactionVolume;

    interactionVolume.setCenterPosition(corner1234);

    //***************   store pointer 1
    interactionVolume.setSubVolumeElement(*isIt->outside(), 0);
    interactionVolume.setIndexOnElement(isIt->indexInOutside(), 0, 0);
    interactionVolume.setIndexOnElement(isIt14->indexInInside(), 0, 1);

    // center of face in global coordinates, i.e., the midpoint of edge 'isIt12'
    const GlobalPosition& globalPosFace12 = isIt->geometry().center();

    // get face volume
    Scalar faceVol12 = isIt->geometry().volume() / 2.0;

    // get outer normal vector scaled with half volume of face 'isIt12'
    Dune::FieldVector<Scalar, dimWorld> unitOuterNormal12 = isIt->centerUnitOuterNormal();
    unitOuterNormal12 *=-1;

    // center of face in global coordinates, i.e., the midpoint of edge 'isIt14'
    const GlobalPosition& globalPosFace41 = isIt14->geometry().center();

    // get face volume
    Scalar faceVol41 = isIt14->geometry().volume() / 2.0;

    // get outer normal vector scaled with half volume of face 'isIt14': for numbering of n see Aavatsmark, Eigestad
    Dune::FieldVector<Scalar, dimWorld> unitOuterNormal14 = isIt14->centerUnitOuterNormal();

    interactionVolume.setNormal(unitOuterNormal12, 0, 0);
    interactionVolume.setNormal(unitOuterNormal14, 0, 1);
    //get the normals of from cell 2 and 4
    unitOuterNormal14 *= -1;
    unitOuterNormal12 *= -1;
    interactionVolume.setFaceArea(faceVol12, 0, 0);
    interactionVolume.setFaceArea(faceVol41, 0, 1);
    interactionVolume.setFacePosition(globalPosFace12, 0, 0);
    interactionVolume.setFacePosition(globalPosFace41, 0, 1);

    // access neighbor cell 2 of 'isIt12'
    ElementPointer elementPointer2 = isIt->inside();

    //****************     store pointer 2
    interactionVolume.setSubVolumeElement(elementPointer2, 1);
//    interactionVolume.setIndexOnElement(isIt->indexInInside(), 1, 1);
    interactionVolume.setNormal(unitOuterNormal12, 1, 1);
    interactionVolume.setFaceArea(faceVol12, 1, 1);
    interactionVolume.setFacePosition(globalPosFace12, 1, 1);


    //****************     data for cell 4
    ElementPointer elementPointer4 = isIt14->outside();

    //store pointer 4
    interactionVolume.setSubVolumeElement(elementPointer4, 3);
//    interactionVolume.setIndexOnElement(globalIdx4, 3, 0);

    interactionVolume.setNormal(unitOuterNormal14, 3, 0);
    interactionVolume.setFaceArea(faceVol41, 3, 0);
    interactionVolume.setFacePosition(globalPosFace41, 3, 0);

    //****************     data for cell 3
    ElementPointer elementPointer3 = isIt23->outside();
    //store pointer 3
    interactionVolume.setSubVolumeElement(elementPointer3, 2);

    GlobalPosition globalPosFace23 = isIt23->geometry().center();
    Scalar faceVol23 = isIt23->geometry().volume() / 2.0;
    // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
    GlobalPosition unitOuterNormal23 = isIt23->centerUnitOuterNormal();

    interactionVolume.setNormal(unitOuterNormal23, 1, 0);
    unitOuterNormal23 *= -1;
    interactionVolume.setNormal(unitOuterNormal23, 2, 1);
    interactionVolume.setFaceArea(faceVol23, 1, 0);
    interactionVolume.setFaceArea(faceVol23, 2, 1);
    Scalar dummy=NAN;
    interactionVolume.setFaceArea(dummy, 2, 0);
    interactionVolume.setFaceArea(dummy, 3, 1);
    interactionVolume.setFacePosition(globalPosFace23, 1, 0);
    interactionVolume.setFacePosition(globalPosFace23, 2, 1);

//    int globalIdx1 = problem().variables().index(*isIt->outside());
//    int globalIdx2 = problem().variables().index(*elementPointer2);
//    int globalIdx3 = problem().variables().index(*elementPointer3);
//    int globalIdx4 = problem().variables().index(*isIt14->outside());

    Dune::FieldVector<Scalar, dim> unity(1.);
    std::vector<Dune::FieldVector<Scalar, dim> > lambda(4, unity);

    Dune::FieldMatrix<Scalar,dim,2*dim-dim+1> T;
    int triangleType = pressureModelAdaptive2p_->calculateTransmissibility(
            T, interactionVolume, lambda,
            0, 1, 2, 3);

    // 3.decide which triangle (which transmissibility coefficients) to use
    if (triangleType == pressureModelAdaptive2p_->rightTriangle)
    {
        additionalIntersectionIt = isIt23;
        // Translate flux from 2p mpfa-l local indexing to
        // 2p2c naming with i and j (reverse direction)
        // and adjust for the fact that it is right Triangle
        additionalT[0] = -T[1][2];
        additionalT[1] = -T[1][1];
        additionalT[2] = -T[1][0];
    }
    else if (triangleType == pressureModelAdaptive2p_->leftTriangle)
    {
        additionalIntersectionIt = isIt14;
        // Translate flux from 2p mpfa-l local indexing to
        // 2p2c naming with i and j (reverse direction)
        // and adjust for the fact that it is left Triangle
        additionalT[0] = -T[1][0];
        additionalT[1] = -T[1][1];
        additionalT[2] = -T[1][2];

    }
    else
        DUNE_THROW(Dune::MathError, "No transmissivity for second half edge found!");

    return triangleType;
}


}//end namespace Dumux
#endif
