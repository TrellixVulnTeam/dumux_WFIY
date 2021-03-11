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
 * \brief A test problem for the one-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROBLEM_HH
#define DUMUX_PNM2P_PROBLEM_HH

#define CUSTOMCOMPONENTS 0

#include <dune/foamgrid/foamgrid.hh>

// base problem
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porenetwork/2pnc/model.hh>

// spatial params
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/porenetwork/common/utilities.hh>

// #if !CUSTOMCOMPONENTS
#include <dumux/material/fluidsystems/h2oair.hh>
// #else
// #include "./files/h2oairfluidsystem.hh"
// #endif

namespace Dumux
{
template <class TypeTag>
class DrainageProblem;

namespace Properties
{

namespace TTag {
#if ISOTHERMAL
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoPNC>; };
#else
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoPNC>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = Dumux::DrainageProblem<TypeTag>; };

// SET_PROP(DrainageProblem, FluidSystem)
//   {
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
// #if !CUSTOMCOMPONENTS
//     typedef Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar> > type;
// #else
//     typedef Dumux::FluidSystems::MyH2OAir<TypeTag, Scalar, Dumux::Components::SimpleH2O<Scalar> > type;
// #endif
//   };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
  {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
  };


// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };


// // SET_INT_PROP(DrainageProblem,
// //              Formulation,
// //              TwoPFormulation::pwsn);

// // Set the allowed discrepancy for the solution plausibility check
// SET_SCALAR_PROP(DrainageProblem, CheckSolutionEpsilon, 5e-9);
}

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        // grid and world dimension
        dim = GridView::dimension,
        dimworld = GridView::dimensionworld,

        // primary variable indices
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::switchIdx,
//         pnIdx = Indices::pnIdx,
//         swIdx = Indices::swIdx,


        // phase indices
        wPhaseIdx = Indices::conti0EqIdx + FluidSystem::H2OIdx,
        nPhaseIdx = Indices::conti0EqIdx + FluidSystem::AirIdx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

    enum { dofCodim =  dim };

    // phase presence
    enum { wPhaseOnly = Indices::firstPhaseOnly };

public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams), eps_(1e-7)
    {
        verbose_ = getParam<bool>("Problem.Verbose", false);
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        initialPc_ = getParam<Scalar>("Problem.InitialPc");
        finalPc_ = getParam<Scalar>("Problem.FinalPc");
        numSteps_ = getParam<int>("Problem.NumSteps");
        episodeLenght_ = getParam<Scalar>("Problem.EpisodeLength");

        swEquilibrium_.resize(3);
        swEquilibrium_ = {1.0, 1.0, 1.0};


        pcEpisopde_.resize(numSteps_);
        for (int i = 0 ; i < pcEpisopde_.size(); i++)
              pcEpisopde_[i] = initialPc_ + i*(finalPc_ - initialPc_)/numSteps_;

        std::cout << "The following global PCs are applied: " << std::endl;
        for (auto x: pcEpisopde_)
        {
            std::cout << x << std::endl;
        }

        logfile_.open("logfile_" + this->name() + ".txt"); //for the logfile
        logfile_ <<"Logfile for: " + this->name()  << std::endl;
        logfile_ << std::left << std::setw(20) << std::setfill(' ') << "Episode"
                 << std::left << std::setw(20) << std::setfill(' ') << "Time"
                 << std::left << std::setw(20) << std::setfill(' ') << "globalPc"
                 << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                 << std::endl;

        step_ = 0;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        if (vtpOutputFrequency_ < 0)
            return true;

        if (vtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return (timeStepIndex % vtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
    }


    /*!
    //  * \brief Called at the end of each time step
    //  */
    // template<class SolutionVector>
    // void postTimeStep(const GridVariables& gridVariables, const SolutionVector& sol)
    // {
    //     const auto avg = Averaging<TypeTag>::averagedValues(*this, gridVariables, sol);

    //     logfile_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << "nan"
    //                            << std::left << std::setw(20) << std::setfill(' ') << "nan"
    //                            << std::left << std::setw(20) << std::setfill(' ') <<  pcEpisopde_[step_]
    //                            << std::left << std::setw(20) << std::setfill(' ') << avg.sw
    //                            << std::endl;

    //     // store the three most recent averaged saturations
    //     std::rotate(swEquilibrium_.rbegin(), swEquilibrium_.rbegin()+1, swEquilibrium_.rend());
    //     swEquilibrium_[0]= avg.sw;

    //     // Check for steady state and end episode
    //     const Scalar change = std::abs(swEquilibrium_[0]-swEquilibrium_[2])/swEquilibrium_[0];
    //     const Scalar pc = pcEpisopde_[step_];
    //     std::cout << "global pC applied: " << pc << " / " << finalPc_ << " (step " << step_ << " of " << numSteps_ << ")" << std::endl;
    //     std::cout << "swAverage: " << swEquilibrium_[0] << " (relative shift: " << change << ")" << std::endl;
    //     if(change < 1e-8)
    //         ++step_;
    // }


#if ISOTHERMAL
    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}
#endif
     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv))
           bcTypes.setAllDirichlet();

        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();
#if !ISOTHERMAL
        bcTypes.setDirichlet(temperatureIdx);
#endif
        return bcTypes;
    }

    /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        control volume.
    *
    * \param values The dirichlet values for the primary variables
    * \param vertex The vertex (pore body) for which the condition is evaluated
    *
    * For this method, the \a values parameter stores primary variables.
    */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values.setState(Indices::bothPhases);
        // values[snIdx] = 0.1;

        // int dofIdxGlobal = this->fvGridGeometry().vertexMapper().index(vertex);
        // if(this->spatialParams().isOutletPore(dofIdxGlobal))
        // {
            // values.setState(Indices::secondPhaseOnly);
            values[Indices::switchIdx] = 0.90;
            // values[Indices::switchIdx] = 1e-4;

//


#if !ISOTHERMAL
        if(isInletPore_(scv))
            values[temperatureIdx] = 273.15 + 15;
        else
            values[temperatureIdx] = 273.15 + 10;
#endif
        return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        vertex
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * units of \f$ [ \textnormal{unit of conserved quantity} / s] \f$
     * \param vertex The vertex
     * \param volVars All volume variables for the pore
     *
     * Positive values mean that mass is created,
     * negative ones mean that it vanishes.
     */
     PrimaryVariables source(const Element &element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolume &scv) const
     {
        PrimaryVariables values(0.0);
        // const int vIdx =  this->fvGridGeometry().vertexMapper().index(vertex);
        //
        // if(this->spatialParams().isInletPore(vIdx))
        // {
        //     const Scalar pc = pcEpisopde_[step_];
        //     values[snIdx] = (volVars.pressure(nPhaseIdx) - (1e5 + pc)) * 1e8;
        // }

        return values;
    }
    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        // priVars[Indices::witchIdx] = 0.0;

        values.setState(Indices::firstPhaseOnly);

        // get global index of pore
        int dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if(isInletPore_(dofIdxGlobal))
        {
            values.setState(Indices::bothPhases);
            values[Indices::switchIdx] = 0.90;
            // values[Indices::switchIdx] = 1e-4;
        }
        else
            values[Indices::switchIdx] = 1e-8;

#if !ISOTHERMAL
        values[temperatureIdx] = 273.15 + 10;
#endif
        return values;
    }

    /*!
     * \brief Evaluate the initial invasion state of a pore throat
     *
     * Returns true for a invaded throat and false elsewise.
     */
    bool initialInvasionState(const Element &element) const
    {
        // if(this->spatialParams().isOutletThroat(element))
        //     return true;
        return false;
    }

    // \}

    /*!
     * \brief Calls gnuplot and plots the pc-S curve
     */
    void plotPcS()
    {
        std::FILE * pipe_;
        pipe_ = popen("gnuplot -persist", "w");
        std::string command = "set xrange [0:1] \n";
        command += "set xlabel 'S_w' \n";
        command += "set ylabel 'p_c' \n";
        std::string filename = "'logfile_"+ this->name() +".txt'";
        command += "plot " + filename + " using 4:3 with lines title 'p_c(S_w)'";
        fputs((command + "\n").c_str(), pipe_);
        pclose(pipe_);
    }


    const bool verbose() const
    { return verbose_; }

    bool simulationFinished() const
    { return step_ >= numSteps_ ; }

        bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

private:
    Scalar eps_;
    bool verbose_;
    int vtpOutputFrequency_;
    Scalar initialPc_;
    Scalar finalPc_;
    int numSteps_;
    Scalar episodeLenght_;
    std::vector<Scalar> pcEpisopde_;
    std::vector<Scalar> swEquilibrium_;
    std::ofstream logfile_;

    int step_;
};
} //end namespace

#endif
