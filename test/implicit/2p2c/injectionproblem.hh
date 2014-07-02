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
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dune/common/version.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

#include <dumux/implicit/2p2c/2p2cmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

#include "injectionspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(TwoPTwoC, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem, INHERITS_FROM(BoxModel, InjectionProblem));
NEW_TYPE_TAG(InjectionCCProblem, INHERITS_FROM(CCModel, InjectionProblem));
    
// Set the grid type
SET_PROP(InjectionProblem, Grid)
{
    typedef Dune::SGrid<2,2> type;
};

// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dumux::InjectionProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(InjectionProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useComplexRelations = false;
 public:
    typedef Dumux::FluidSystems::H2ON2<Scalar, useComplexRelations> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem, ProblemEnableGravity, true);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(InjectionProblem, UseMoles, true);

SET_BOOL_PROP(InjectionProblem, ImplicitEnableJacobianRecycling, true);
SET_BOOL_PROP(InjectionProblem, VtkAddVelocity, false);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one (\f$ K=10e-12\f$) for \f$ y<22m\f$ and one with a lower permeablility (\f$ K=10e-13\f$)
 * in the rest of the domain.
 *
 * A mixture of Nitrogen and Water vapor, which is composed according to the prevailing conditions (temperature, pressure)
 * enters a water-filled aquifer. This is realized with a solution-dependent Neumann boundary condition at the right boundary
 * (\f$ 5m<y<15m\f$). The aquifer is situated 2700m below sea level. The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the lower permeable aquitard.
 * 
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2c -parameterFile ./test_box2p2c.input</tt> or
 * <tt>./test_cc2p2c -parameterFile ./test_cc2p2c.input</tt>
 */
template <class TypeTag>
class InjectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        contiH2OEqIdx = Indices::contiWEqIdx,
        contiN2EqIdx = Indices::contiNEqIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem(TimeManager &timeManager,
                     const GridView &gridView)
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        : ParentType(timeManager, GridCreator::grid().leafGridView())
#else
        : ParentType(timeManager, GridCreator::grid().leafView())
#endif
    {
        try
        {
            nTemperature_       = GET_RUNTIME_PARAM(TypeTag, int, Problem.NTemperature);
            nPressure_          = GET_RUNTIME_PARAM(TypeTag, int, Problem.NPressure);
            pressureLow_        = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PressureLow);
            pressureHigh_       = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PressureHigh);
            temperatureLow_     = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.TemperatureLow);
            temperatureHigh_    = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.TemperatureHigh);
            temperature_        = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InitialTemperature);
            depthBOR_           = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.DepthBOR);
            name_               = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

        /* Alternative syntax:
         * typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
         * const Dune::ParameterTree &tree = ParameterTree::tree();
         * nTemperature_       = tree.template get<int>("Problem.NTemperature");
         *
         * + We see what we do
         * - Reporting whether it was used does not work
         * - Overwriting on command line not possible
         */


        eps_ = 1e-6;

        // initialize the tables of the fluid system
        //FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        //stateing in the console whether mole or mass fractions are used
        if(!useMoles)
        {
        	std::cout<<"problem uses mass-fractions"<<std::endl;
        }
        else
        {
        	std::cout<<"problem uses mole-fractions"<<std::endl;
        }
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageW, storageN;
        this->model().globalPhaseStorage(storageW, wPhaseIdx);
        this->model().globalPhaseStorage(storageN, nPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: wetting=[" << storageW << "]"
                     << " nonwetting=[" << storageN << "]\n";
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
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K]
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
         values = 0;

         GlobalPosition globalPos;
         if (isBox)
             globalPos = element.geometry().corner(scvIdx);
         else
             globalPos = intersection.geometry().center();

         Scalar injectedPhaseMass = 1e-3;
         Scalar moleFracW = elemVolVars[scvIdx].moleFraction(nPhaseIdx, wCompIdx);
         if (globalPos[1] < 15 && globalPos[1] > 7) {
             values[contiN2EqIdx] = -(1-moleFracW)*injectedPhaseMass/FluidSystem::molarMass(nCompIdx); //mole/(m^2*s) -> kg/(s*m^2)
             values[contiH2OEqIdx] = -moleFracW*injectedPhaseMass/FluidSystem::molarMass(wCompIdx); //mole/(m^2*s) -> kg/(s*m^2)
         }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Returns the initial phase state for a control volume.
     *
     * \param vertex The vertex
     * \param globalIdx The global index of the vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    { return Indices::wPhaseOnly; }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = FluidSystem::H2O::liquidDensity(temperature_, 1e5);

        Scalar pl = 1e5 - densityW*this->gravity()[1]*(depthBOR_ - globalPos[1]);
        Scalar moleFracLiquidN2 = pl*0.95/BinaryCoeff::H2O_N2::henry(temperature_);
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidN2;

        Scalar meanM =
            FluidSystem::molarMass(wCompIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(nCompIdx)*moleFracLiquidN2;
        if(!useMoles) //mass fraction formulation
        {
        	Scalar massFracLiquidN2 = moleFracLiquidN2*FluidSystem::molarMass(nCompIdx)/meanM;
        	values[Indices::switchIdx] = massFracLiquidN2;
        }
        else //mole-fraction formulation
        {
        	values[Indices::switchIdx] = moleFracLiquidN2;
        }
        values[Indices::pressureIdx] = pl;
    }

    Scalar temperature_;
    Scalar depthBOR_;
    Scalar eps_;

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;


};
} //end namespace

#endif
