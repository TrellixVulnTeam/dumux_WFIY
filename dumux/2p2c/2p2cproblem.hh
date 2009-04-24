// $Id$

#ifndef DUNE_2P2CPROBLEM_HH
#define DUNE_2P2CPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/property_baseclasses.hh>
/**
 * @file
 * @brief  Base class for defining an instance of the 2P2C problem
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{
//! base class that defines the parameters of a diffusion equation
/*! An interface for defining parameters for the stationary diffusion equation
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation.
 *
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar>
class TwoPTwoCProblem {
    enum {dim=Grid::dimension, numEq=2};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

public:
    //    //! evaluate diffusion tensor
    //    /*! Evaluate the diffusion tensor at given location
    //      @param[in]  globalPos    position in global coordinates
    //      @param[in]  element    entity of codim 0
    //      @param[in]  localPos   position in reference element of element
    //      @param[out] D    diffusion tensor to be filled
    //     */
    //    virtual const FieldMatrix<Scalar,dim,dim>& K (const FieldVector<Scalar,dim>& globalPos) = 0;

    //! evaluate source term
    /*! evaluate source term at given location
      @param[in]  globalPos    position in global coordinates
      @param[in]  element    entity of codim 0
      @param[in]  localPos   position in reference element of element
      \return     value of source term
    */
    virtual FieldVector<Scalar,numEq> q (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                         const FieldVector<Scalar,dim>& localPos) const = 0;

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition type given by enum in this class
    */
    //    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<Scalar,n>& globalPos, const Element& element,
    //            const IntersectionIterator& intersectionIt,
    //            const FieldVector<Scalar,n>& localPos) const = 0;

    virtual FieldVector<BoundaryConditions::Flags, numEq>bctype(
                                                                const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                                                const IntersectionIterator& intersectionIt,
                                                                const FieldVector<Scalar,dim>& localPos) const = 0;

    //! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
    /*! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
      @param[in]  globalPos    position in global coordinates
      \return     index of the primary variable
    */

    virtual void dirichletIndex(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<Scalar,dim>& localPos, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
    }

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,numEq>dirichlet(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,numEq> neumann(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                         const IntersectionIterator& intersectionIt,
                                         const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate initial condition at given position
    /*! evaluate initial boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,numEq> initial(const FieldVector<Scalar,dim>& globalPos,
                                              const Element& element, const FieldVector<Scalar,dim>& localPos) const = 0;

    //! initiate phase state at given position
    /*! initiate phase state at given position
      @param[in]  globalPos    position in global coordinates
      \return     initial phase state
    */
    virtual int initialPhaseState(const FieldVector<Scalar,dim>& globalPos,
                                  const Element& element, const FieldVector<Scalar,dim>& localPos) const = 0;

    //    virtual double porosity(const FieldVector<Scalar,dim>& globalPos, const Element& element,
    //            const FieldVector<Scalar,dim>& localPos) const = 0;

    virtual FieldVector<Scalar,dim> gravity() const= 0;

    virtual double depthBOR () const = 0;

    //! properties of the wetting (liquid) phase
    /*! properties of the wetting (liquid) phase
      \return    wetting phase
    */
    virtual Liquid_GL& wettingPhase () const
    {
        return wettingPhase_;
    }

    //! properties of the nonwetting (liquid) phase
    /*! properties of the nonwetting (liquid) phase
      \return    nonwetting phase
    */
    virtual Gas_GL& nonwettingPhase () const
    {
        return nonwettingPhase_;
    }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    virtual Matrix2p<Grid, Scalar>& soil () const
    {
        return soil_;
    }

    //! object for multicomponent calculations
    /*! object for multicomponent calculations including mass fractions,
     * mole fractions and some basic laws
     \return    multicomponent object
    */
    virtual MultiComp& multicomp () const
    {
        return multicomp_;
    }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    virtual TwoPhaseRelations<Grid, Scalar>& materialLaw () const
    {
        return materialLaw_;
    }

    //element-wise return of the values of an Exact solution
    virtual Scalar uExOutVertex(int &ElementIndex, int VariableIndex) const {
        DUNE_THROW(NotImplemented, "Ex(akt) Solution");
        return 0;
    }

    //updates an exact/analytic solution
    virtual void updateExSol(double &dt,
                             BlockVector<FieldVector<Scalar, numEq> > &approxSol) {
        DUNE_THROW(NotImplemented, "Ex(akt) Solution");
        return;
    }


    TwoPTwoCProblem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, Scalar>& soil,
                    MultiComp& multicomp = *(new CWaterAir),
                    TwoPhaseRelations<Grid,Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>),
                    const bool exsol = false)
        : exsolution(exsol), wettingPhase_(liq), nonwettingPhase_(gas), soil_(soil),
          multicomp_(multicomp), materialLaw_(materialLaw)
    {    }

    //! always define virtual destructor in abstract base class
    virtual ~TwoPTwoCProblem () {}

    const bool exsolution;

protected:
    Liquid_GL& wettingPhase_;
    Gas_GL& nonwettingPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    MultiComp& multicomp_;
    TwoPhaseRelations<Grid, Scalar>& materialLaw_;
};

}
#endif
