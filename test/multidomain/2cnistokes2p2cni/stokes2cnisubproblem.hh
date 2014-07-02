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
/**
 * \file
 * \ingroup Stokes2cniProblems
 * \brief Non-isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 */
#ifndef DUMUX_STOKES2CNI_SUBPROBLEM_HH
#define DUMUX_STOKES2CNI_SUBPROBLEM_HH

#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dumux/freeflow/stokesncni/stokesncnimodel.hh>
#include <dumux/multidomain/couplinglocalresiduals/stokesncnicouplinglocalresidual.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>

#include "2cnistokes2p2cnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class Stokes2cniSubProblem;

//////////
// Specify the properties for the Stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cniSubProblem, 
	INHERITS_FROM(BoxStokesncni, SubDomain, TwoCNIStokesTwoPTwoCNISpatialParams));

//SET_PROP(Stokes2cniSubProblem, ModelParameterGroup)
//{ static const char *value() { return "FF"; }; };

// Set the problem property
SET_TYPE_PROP(Stokes2cniSubProblem, Problem, Dumux::Stokes2cniSubProblem<TypeTag>);

//! Use the Stokes2cniCouplingLocalResidual for the computation of the local residual in the Stokes domain
SET_TYPE_PROP(Stokes2cniSubProblem, LocalResidual, StokesncniCouplingLocalResidual<TypeTag>);

//! Set the property for the material parameters by extracting it from the material law
SET_PROP(Stokes2cniSubProblem, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
public:
    typedef typename MaterialLaw::Params type;
};

SET_TYPE_PROP(Stokes2cniSubProblem, Constraints,
		Dune::PDELab::NonoverlappingConformingDirichletConstraints);

// set the grid operator
SET_PROP(Stokes2cniSubProblem, GridOperator)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsTrafo) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef Dumux::PDELab::MultiDomainLocalOperator<TypeTag> LocalOperator;

    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
public:
    typedef Dune::PDELab::GridOperator<GridFunctionSpace,
            GridFunctionSpace, LocalOperator,
            Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
            Scalar, Scalar, Scalar,
            ConstraintsTrafo, ConstraintsTrafo,
            true> type;
};

SET_PROP(Stokes2cniSubProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) CoupledTypeTag;
    typedef typename GET_PROP_TYPE(CoupledTypeTag, FluidSystem) FluidSystem;
public:
    typedef FluidSystem type;
};

// use formulation based on mass fractions
SET_BOOL_PROP(Stokes2cniSubProblem, UseMoles, false);

// Disable gravity in the Stokes domain
SET_BOOL_PROP(Stokes2cniSubProblem, ProblemEnableGravity, false);

// switch inertia term on or off
SET_BOOL_PROP(Stokes2cniSubProblem, EnableNavierStokes, false);
}

/*!
 * \ingroup BoxStokesncniModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal two-component stokes subproblem with air flowing
 *        from the left to the right and coupling at the bottom.
 *
 * The Stokes subdomain is sized 0.25m times 0.25m. The boundary conditions
 * for the momentum balances are all set to Dirichlet, except on the right
 * boundary, where outflow conditions are set. The mass balance receives
 * outflow BCs, which are replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex receives Dirichlet BCs, to set the pressure level.
 *
 * This sub problem uses the \ref StokesncniModel. It is part of the
 * 2cnistokes2p2cni model and is combined with the 2p2cnisubproblem for
 * the Darcy domain.
 */
template <class TypeTag>
class Stokes2cniSubProblem : public StokesProblem<TypeTag>
{
    typedef Stokes2cniSubProblem<TypeTag> ThisType;
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // soil parameters for beavers & joseph
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    enum {
        // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // equation indices
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx, //!< Index of the transport equation (massfraction)
        energyEqIdx =    Indices::energyEqIdx     //!< Index of the energy equation (temperature)
    };
    enum { // primary variable indices
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        velocityZIdx = Indices::velocityZIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum { phaseIdx = Indices::phaseIdx };
    enum { numComponents = Indices::numComponents };
    enum {
        transportCompIdx = Indices::transportCompIdx, //!< water component index
        phaseCompIdx = Indices::phaseCompIdx          //!< air component index
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

    /*!
     * \brief The sub-problem for the Stokes subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
public:
    Stokes2cniSubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          spatialParams_(gridView)
    {
        try
        {
            bboxMin_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMin);
            bboxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMax);
            bboxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);
            bboxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, YMax);

            refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
            refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefPressure);
            refMassfrac_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefMassfrac);
            vxMax_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, VxMax);
            bjSlipVel_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, BeaversJosephSlipVel);

            sinusVelVar_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusVelVar);
            sinusPVar_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusPVar);
            sinusTVar_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusTVar);
            sinusXVar_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, SinusXVar);

            xMaterialInterface_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MaterialInterfaceX);
            runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, RunUpDistanceX); // first part of the interface without coupling
            initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

            alphaBJ_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, AlphaBJ);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Vtk, NameFF); }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given vertex
     *
     * \param values Stores the value of the boundary type
     * \param vertex The vertex
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition &globalPos = vertex.geometry().center();
        const Scalar time = this->timeManager().time();

        values.setAllDirichlet();

        if (onUpperBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);
            values.setNeumann(energyEqIdx);
        }

		// Left inflow boundaries should be Neumann, otherwise the
		// evaporative fluxes are much more grid dependent
        if (onLeftBoundary_(globalPos))
        {
            values.setNeumann(transportEqIdx);
            values.setNeumann(energyEqIdx);

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onRightBoundary_(globalPos))
        {
            values.setAllOutflow();

            if (onUpperBoundary_(globalPos)) // corner point
                values.setAllDirichlet();
        }

        if (onLowerBoundary_(globalPos))
        {
            values.setAllDirichlet();
            if (!onLeftBoundary_(globalPos)) // is this required?
            {
                values.setNeumann(transportEqIdx);
                values.setNeumann(energyEqIdx);
            }

            if (globalPos[0] > runUpDistanceX_-eps_ && time > initializationTime_)
            {
                values.setAllCouplingOutflow();
//                values.setCouplingInflow(energyEqIdx);
            }
        }

        // the mass balance has to be of type outflow
        // it does not get a coupling condition, since pn is a condition for stokes
        values.setOutflow(massBalanceIdx);

        // set pressure at one point, do NOT specify this
        // if the Darcy domain has a Dirichlet condition for pressure
        if (onRightBoundary_(globalPos))
        {
            if (time > initializationTime_)
                values.setDirichlet(pressureIdx, massBalanceIdx);
            else
                if (!onLowerBoundary_(globalPos) && !onUpperBoundary_(globalPos))
                    values.setDirichlet(pressureIdx, massBalanceIdx);
        }
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet vertex
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param vertex The vertex
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

//        initial_(values, globalPos);
        values = 0.0;

        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.0;
        values[pressureIdx] = refPressure()  +
                density*this->gravity()[1]*(globalPos[1] - bboxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
        values[temperatureIdx] = refTemperature();
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann intersection
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param is The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * Negative values indicate an inflow.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos =
                fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        values = 0.;

        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);
        const Scalar enthalpy =
                FluidSystem::enthalpy(fluidState, phaseIdx);
        const Scalar xVelocity = xVelocity_(globalPos);

        if (onLeftBoundary_(globalPos)
                && globalPos[1] > bboxMin_[1] && globalPos[1] < bboxMax_[1])
        {
            values[transportEqIdx] = -xVelocity*density*refMassfrac();
            values[energyEqIdx] = -xVelocity*density*enthalpy;
        }
    }

    /*!
     * \brief Evaluate the Beavers-Joseph coefficient
     *        at the center of a given intersection
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param is The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * \return Beavers-Joseph coefficient
     */
    Scalar beaversJosephCoeff(const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos =
                fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        if (onLowerBoundary_(globalPos))
            return alphaBJ_;
        else
            return 0.0;
    }

    /*!
     * \brief Returns the intrinsic permeability tensor \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar permeability(const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        return spatialParams_.intrinsicPermeability(element,
                                                    fvGeometry,
                                                    scvIdx);
    }

    // \}

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations in units of
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx) const
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }
    // \}

    /*!
     * \brief Determines if globalPos is a corner of the grid
     *
     * \param globalPos The global position
     */
    bool isCornerPoint(const GlobalPosition &globalPos)
    {
        if ((onLeftBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onLeftBoundary_(globalPos) && onUpperBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onUpperBoundary_(globalPos)))
            return true;
        else
            return false;
    }

    /*!
     * \brief Auxiliary function used for the mortar coupling, if mortar coupling,
     *        this should return true
     *
     * \param globalPos The global position
     */
    bool isInterfaceCornerPoint(const GlobalPosition &globalPos) const
    { return false; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParams &spatialParams()
    { return spatialParams_; }
    const SpatialParams &spatialParams() const
    { return spatialParams_; }

    //! \brief Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_ + diurnalVariation_(sinusPVar_); }

    //! \brief Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_+ diurnalVariation_(sinusTVar_); }

    //! \brief Returns the reference mass fraction.
    const Scalar refMassfrac() const
    { return refMassfrac_ + diurnalVariation_(sinusXVar_); }

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        FluidState fluidState;
        updateFluidStateForBC_(fluidState);

        const Scalar density =
                FluidSystem::density(fluidState, phaseIdx);

        values[velocityXIdx] = xVelocity_(globalPos);
        values[velocityYIdx] = 0.;

        values[pressureIdx] = refPressure()
                + density*this->gravity()[1]*(globalPos[1] - bboxMin_[1]);
        values[massOrMoleFracIdx] = refMassfrac();
        values[temperatureIdx] = refTemperature();
    }

    //! \brief set the profile of the inflow velocity (horizontal direction)
    const Scalar xVelocity_(const GlobalPosition &globalPos) const
    {
        const Scalar vmax = vxMax_ + hourlyVariation_(sinusVelVar_);
//        const Scalar relativeHeight = (globalPos[1]-bboxMin_[1])/height_();
        // linear profile
//        return vmax*relativeHeight + bjSlipVel_; // BJ slip velocity is added as sqrt(Kxx)
	      // parabolic profile
        return  4*vmax*(globalPos[1] - bboxMin_[1])*(bboxMax_[1] - globalPos[1])
                / (height_()*height_()) + bjSlipVel_;
//        Scalar tempVel = (globalPos[1] - bboxMin_[1])/height_();
//        return vmax*tempVel*tempVel + bjSlipVel_;
        // logarithmic profile
//        return 0.1*vmax*log((relativeHeight+1e-3)/1e-3) + bjSlipVel_;
    }

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState) const
    {
        fluidState.setTemperature(refTemperature());
        fluidState.setPressure(phaseIdx, refPressure());

        Scalar massFraction[numComponents];
        massFraction[transportCompIdx] = refMassfrac();
        massFraction[phaseCompIdx] = 1 - massFraction[transportCompIdx];

        // calculate average molar mass of the gas phase
        Scalar M1 = FluidSystem::molarMass(transportCompIdx);
        Scalar M2 = FluidSystem::molarMass(phaseCompIdx);
        Scalar X2 = massFraction[phaseCompIdx];
        Scalar massToMoleDenominator = M2 + X2*(M1 - M2);

        fluidState.setMoleFraction(phaseIdx, transportCompIdx, massFraction[transportCompIdx]*M2/massToMoleDenominator);
        fluidState.setMoleFraction(phaseIdx, phaseCompIdx, massFraction[phaseCompIdx]*M1/massToMoleDenominator);
    }


    const Scalar diurnalVariation_(const Scalar value) const
    {
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        return sin(2*M_PI*time/86400) * value;
    }


    const Scalar hourlyVariation_(const Scalar value) const
    {
        const Scalar time = this->timeManager().time();
        return sin(2*M_PI*time/3600) * value;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bboxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bboxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bboxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bboxMax_[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    const Scalar height_() const
    { return bboxMax_[1] - bboxMin_[1]; }

    // spatial parameters
    SpatialParams spatialParams_;

    static constexpr Scalar eps_ = 1e-8;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    Scalar refPressure_;
    Scalar refTemperature_;
    Scalar refMassfrac_;

    Scalar vxMax_;
    Scalar bjSlipVel_;
    Scalar alphaBJ_;

    Scalar sinusVelVar_;
    Scalar sinusPVar_;
    Scalar sinusTVar_;
    Scalar sinusXVar_;

    Scalar xMaterialInterface_;
    Scalar runUpDistanceX_;
    Scalar initializationTime_;
};
} //end namespace

#endif // DUMUX_STOKES2CNI_SUBPROBLEM_HH
