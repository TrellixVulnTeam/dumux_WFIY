//$Id$
#ifndef DUNE_NEW_BRINECO2PROBLEM2PNI_HH
#define DUNE_NEW_BRINECO2PROBLEM2PNI_HH

//#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/phaseproperties/phaseproperties_brineco2.hh>


#include <dumux/auxiliary/timemanager.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/new_models/2pni/2pniboxmodel.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/nonlinear/new_newtonmethod.hh>
#include <dumux/nonlinear/new_newtoncontroller.hh>
#include<dumux/new_models/2pni/2pninewtoncontroller.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

namespace Dune
{
/*!
 * \todo Please doc me!
 */
template<class ScalarT>
class NewBrineCO2Problem2pni : public BasicDomain<Dune::YaspGrid<2>, ScalarT>
{
    typedef Dune::YaspGrid<2>              Grid;
    typedef BasicDomain<Grid, ScalarT>     ParentType;
    typedef NewBrineCO2Problem2pni<ScalarT>  ThisType;
    typedef TwoPNIBoxModel<ThisType>         Model;

    typedef Dune::Liq_BrineCO2                     WettingPhase;
    typedef Dune::Gas_BrineCO2                     NonwettingPhase;
    typedef Dune::HomogeneousSoil<Grid, ScalarT>  Soil;
    typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the Pw-Sn model
    typedef typename Model::TwoPNITraits          TwoPNITraits;

private:
    // some constants from the traits for convenience
    enum {
        numEq       = BoxTraits::numEq,
        pressureIdx = TwoPNITraits::pressureIdx,
        saturationIdx   = TwoPNITraits::saturationIdx,
        temperatureIdx = TwoPNITraits::temperatureIdx,

        // Grid and world dimension
        dim         = DomainTraits::dim,
        dimWorld    = DomainTraits::dimWorld,
    };

    // copy some types from the traits for convenience
    typedef typename DomainTraits::Scalar                     Scalar;
    typedef typename DomainTraits::Element                    Element;
    typedef typename DomainTraits::ElementIterator            ElementIterator;
    typedef typename DomainTraits::ReferenceElement           ReferenceElement;
    typedef typename DomainTraits::Vertex                     Vertex;
    typedef typename DomainTraits::VertexIterator             VertexIterator;
    typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
    typedef typename DomainTraits::LocalPosition              LocalPosition;
    typedef typename DomainTraits::GlobalPosition             GlobalPosition;

    typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction               SpatialFunction;
    typedef typename BoxTraits::SolutionVector                SolutionVector;
    typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

    typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>                  TimeManager;

    typedef typename Model::NewtonMethod                NewtonMethod;
    typedef Dune::TwoPNINewtonController<NewtonMethod>        NewtonController;

public:
    NewBrineCO2Problem2pni(Grid *grid,
                   Scalar dtInitial,
                   Scalar tEnd)
        : ParentType(grid),
          materialLaw_(soil_, wPhase_, nPhase_),
          timeManager_(tEnd,
                       this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("new_brineco2_2pni")
    {
        timeManager_.setStepSize(dtInitial);

        eps_    = 1e-8;
        depthBOR_ = 3238.20;

        gravity_ = 0;
        gravity_[dim - 1] = -9.81;

        wasRestarted_ = false;
    }

    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // set the initial condition
        model_.initial();

        if (!wasRestarted_) {
            // write the inital solution to disk
            writeCurrentResult_();
        }
    }

    /*!
     * \brief Called by the TimeManager in order to get a time
     *        integration on the model.
     *
     * \note timeStepSize and nextStepSize are references and may
     *       be modified by the TimeIntegration. On exit of this
     *       function 'timeStepSize' must contain the step size
     *       actually used by the time integration for the current
     *       steo, and 'nextStepSize' must contain the suggested
     *       step size for the next time step.
     */
    void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
    {
        model_.update(stepSize,
                      nextStepSize,
                      newtonMethod_,
                      newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        if (this->grid().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        // write the current result to disk
        writeCurrentResult_();

        // write restart file after every five steps
        static int dummy = 0;
        ++dummy;
        if (dummy % 1 == 0)
            serialize();
    };
    ///////////////////////////////////
    // End of simulation control stuff
    ///////////////////////////////////

    ///////////////////////////////////
    // Strings pulled by the TwoPBoxModel during the course of
    // the simulation (-> boundary conditions, initial conditions,
    // etc)
    ///////////////////////////////////
    //! Returns the current time step size in seconds
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    //! Set the time step size in seconds.
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }


    //! properties of the wetting (liquid) phase
    /*! properties of the wetting (liquid) phase
      \return    wetting phase
    */
    const WettingPhase &wettingPhase() const
    { return wPhase_; }

    //! properties of the nonwetting (liquid) phase
    /*! properties of the nonwetting (liquid) phase
      \return    nonwetting phase
    */
    const NonwettingPhase &nonwettingPhase() const
    { return nPhase_; }


    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    const Soil &soil() const
    {  return soil_; }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    Soil &soil()
    {  return soil_; }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    MaterialLaw &materialLaw ()
    { return materialLaw_; }

    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = BoundaryConditions::neumann;

        if (globalPos[0] > 99.)
            values = BoundaryConditions::dirichlet;
        values[temperatureIdx] = BoundaryConditions::dirichlet;
    }

    /////////////////////////////
    // DIRICHLET boundaries
    /////////////////////////////
    void dirichlet(SolutionVector             &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //                const LocalPosition &localPos
        //                    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        initial_(values, globalPos);
    }

    /////////////////////////////
    // NEUMANN boundaries
    /////////////////////////////
    void neumann(SolutionVector             &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0;
        // negative values for injection
        if (globalPos[1] < 10.0 && globalPos[0] < 1.0)
        	values[saturationIdx] = -0.005;
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(SolutionVector &values,
                const Element &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        values = Scalar(0.0);
    }

    //////////////////////////////

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    void initial(SolutionVector         &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        /*                const LocalPosition &localPos
                          = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);
        */

        initial_(values, globalPos);
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(SolutionVector       &values,
                  const GlobalPosition &globalPos) const
    {
           Scalar densityB = 1037.24;
           Scalar pRef = 101300.;
           Scalar TRef = 283.15;

            values[pressureIdx] = pRef - (depthBOR_ - globalPos[dim - 1]) * densityB*gravity_[1];
            values[saturationIdx] = 0.0;
            values[temperatureIdx] = TRef + (depthBOR_ - globalPos[dim - 1]) * 0.03;
    }
    const Model &model() const
    {
        return model_;
    }

    void serialize()
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.serializeBegin(this->grid(),
                           "new2pni",
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

        res.serializeEnd();
    }

    void deserialize(double t)
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.deserializeBegin(this->grid(), "new2pni", t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };

public:
    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    double depthBOR () const
    {
        return depthBOR_;
    }

    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };

private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }


    Scalar eps_;
    Scalar depthBOR_;
    GlobalPosition  gravity_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

    bool wasRestarted_;
};
} //end namespace

#endif
