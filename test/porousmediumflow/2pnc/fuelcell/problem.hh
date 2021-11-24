// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup TwoPNCTests
 * \brief Definition of a problem for water management in PEM fuel cells.
 */

#ifndef DUMUX_FUELCELL_PROBLEM_HH
#define DUMUX_FUELCELL_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

#ifdef NONISOTHERMAL
#include <dumux/material/chemistry/electrochemistry/electrochemistryni.hh>
#else
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>
#endif

namespace Dumux {

/*!
 * \ingroup TwoPNCTests
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class FuelCellProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    // Select the electrochemistry method
#ifdef NONISOTHERMAL
    using ElectroChemistry = typename Dumux::ElectroChemistryNI<Scalar, Indices, FluidSystem, GridGeometry, ElectroChemistryModel::Ochs>;
#else
    using ElectroChemistry = typename Dumux::ElectroChemistry<Scalar, Indices, FluidSystem, GridGeometry, ElectroChemistryModel::Ochs>;
#endif
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = GridGeometry::discMethod == DiscretizationMethods::box;
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    enum { dofCodim = isBox ? dim : 0 };
public:
    FuelCellProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        nTemperature_       = getParam<int>("Problem.NTemperature");
        nPressure_          = getParam<int>("Problem.NPressure");
        pressureLow_        = getParam<Scalar>("Problem.PressureLow");
        pressureHigh_       = getParam<Scalar>("Problem.PressureHigh");
        temperatureLow_     = getParam<Scalar>("Problem.TemperatureLow");
        temperatureHigh_    = getParam<Scalar>("Problem.TemperatureHigh");

        name_               = getParam<std::string>("Problem.Name");

        pO2Inlet_           = getParam<Scalar>("ElectroChemistry.pO2Inlet");

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \name Problem parameters
     */

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scv.dofPosition();

        //reaction sources from electro chemistry
        if(inReactionLayer_(globalPos))
        {
            const auto& volVars = elemVolVars[scv];
            auto currentDensity = ElectroChemistry::calculateCurrentDensity(volVars);
            ElectroChemistry::reactionSource(values, currentDensity);
        }

        return values;
    }


    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();

        if (onUpperBoundary_(globalPos)){
            bcTypes.setAllDirichlet();
        }
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        auto priVars = initial_(globalPos);

        if(onUpperBoundary_(globalPos))
        {
            Scalar pn = 1.0e5;
            priVars[Indices::pressureIdx] = pn;
            priVars[Indices::switchIdx] = 0.3;//Sw for bothPhases
            priVars[Indices::switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
#ifdef NONISOTHERMAL
            priVars[Indices::temperatureIdx] = this->spatialParams().temperatureAtPos(globalPos);
#endif
        }

        return priVars;
    }

    /*!
     * \name Volume terms
     */


    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); }


    /*!
     * \brief Adds additional VTK output data to the VTKWriter.
     *
     * Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        const auto& gridView = this->gridGeometry().gridView();
        currentDensity_.resize(gridView.size(dofCodim));
        reactionSourceH2O_.resize(gridView.size(dofCodim));
        reactionSourceO2_.resize(gridView.size(dofCodim));
        Kxx_.resize(gridView.size(dofCodim));
        Kyy_.resize(gridView.size(dofCodim));

        vtk.addField(currentDensity_, "currentDensity [A/cm^2]");
        vtk.addField(reactionSourceH2O_, "reactionSourceH2O [mol/(sm^2)]");
        vtk.addField(reactionSourceO2_, "reactionSourceO2 [mol/(sm^2)]");
        vtk.addField(Kxx_, "Kxx");
        vtk.addField(Kyy_, "Kyy");
    }

    void updateVtkFields(const SolutionVector& curSol)
    {
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto& globalPos = scv.dofPosition();
                const auto dofIdxGlobal = scv.dofIndex();

                if(inReactionLayer_(globalPos))
                {
                    //reactionSource Output
                    PrimaryVariables source;
                    auto i = ElectroChemistry::calculateCurrentDensity(volVars);
                    ElectroChemistry::reactionSource(source, i);

                    reactionSourceH2O_[dofIdxGlobal] = source[Indices::conti0EqIdx + FluidSystem::H2OIdx];
                    reactionSourceO2_[dofIdxGlobal] = source[Indices::conti0EqIdx + FluidSystem::O2Idx];

                    //Current Output in A/cm^2
                    currentDensity_[dofIdxGlobal] = i/10000;
                }
                else
                {
                    reactionSourceH2O_[dofIdxGlobal] = 0.0;
                    reactionSourceO2_[dofIdxGlobal] = 0.0;
                    currentDensity_[dofIdxGlobal] = 0.0;
                }
                Kxx_[dofIdxGlobal] = volVars.permeability()[0][0];
                Kyy_[dofIdxGlobal] = volVars.permeability()[1][1];
            }
        }
    }

private:

    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::bothPhases);

        Scalar pn = 1.0e5;
        priVars[Indices::pressureIdx] = pn;
        priVars[Indices::switchIdx] = 0.3;//Sw for bothPhases
        priVars[Indices::switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
#ifdef NONISOTHERMAL
        priVars[Indices::temperatureIdx] = this->spatialParams().temperatureAtPos(globalPos);
#endif

        return priVars;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    bool inReactionLayer_(const GlobalPosition& globalPos) const
    { return globalPos[1] < 0.1*(this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1]) + eps_; }

    static constexpr Scalar eps_ = 1e-6;
    int nTemperature_;
    int nPressure_;
    std::string name_ ;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar pO2Inlet_;
    std::vector<double> currentDensity_;
    std::vector<double> reactionSourceH2O_;
    std::vector<double> reactionSourceO2_;
    std::vector<double> Kxx_;
    std::vector<double> Kyy_;
};

} // end namespace Dumux

#endif
