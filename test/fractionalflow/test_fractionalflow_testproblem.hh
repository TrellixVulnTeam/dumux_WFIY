// $Id$
#ifndef TEST_FRACTIONALFLOW_TESTPROBLEM_HH
#define TEST_FRACTIONALFLOW_TESTPROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune {
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class GridView, class Scalar, class VC> class FractionalFlowTestProblem :
        public FractionalFlowProblem<GridView,Scalar,VC> {

    enum {dim=GridView::dimension, dimWorld=GridView::dimensionworld};
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    FractionalFlowTestProblem(VC& variables, Fluid& wettingPhase, Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const GlobalPosition LowerLeft = 0,
                const GlobalPosition UpperRight = 0) :
        FractionalFlowProblem<GridView, Scalar, VC>(variables, wettingPhase, nonWettingPhase, soil, materialLaw), Left_(LowerLeft[0]),
        Right_(UpperRight[0]), eps_(1e-8)
    {}

    virtual Scalar sourcePress  (const GlobalPosition& globalPos, const Element& element,
                                 const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
                                                   const Element& element, const LocalPosition& localPos) const {
        if ((globalPos[0] < eps_))
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                         const LocalPosition& localPos) const
    {
//        if (globalPos[0] > (Right_ - eps_) || globalPos[0] < eps_)
            if (globalPos[0] < eps_)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
                          const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return 2e5;
        // all other boundaries
        return 2e5;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[0] > Right_ - eps_)
            return 3e-7;
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos, Scalar factor) const {
        if (globalPos[0] > Right_ - eps_)
            return factor;
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                    const LocalPosition& localPos) const
    {
            return 0.0;
    }

private:
    Scalar Left_;
    Scalar Right_;

    Scalar eps_;
};
}

#endif
