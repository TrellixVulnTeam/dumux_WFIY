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
 * \brief Base class for fully-implicit models
 */
#ifndef DUMUX_IMPLICIT_MODEL_HH
#define DUMUX_IMPLICIT_MODEL_HH

#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/istl/bvector.hh>

#include "implicitproperties.hh"
#include <dumux/common/valgrind.hh>
#include <dumux/parallel/vertexhandles.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The base class for the vertex centered finite volume
 *        discretization scheme.
 */
template<class TypeTag>
class ImplicitModel
{
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, LocalJacobian) LocalJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) LocalResidual;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    // copying a model is not a good idea
    ImplicitModel(const ImplicitModel &);

public:
    /*!
     * \brief The constructor.
     */
    ImplicitModel()
    : problemPtr_(0)
    {
        enableHints_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Implicit, EnableHints);
    }

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &problem)
    {
        problemPtr_ = &problem;

        updateBoundaryIndices_();

        int numDofs = asImp_().numDofs();
        uCur_.resize(numDofs);
        uPrev_.resize(numDofs);
        if (isBox)
            boxVolume_.resize(numDofs);

        localJacobian_.init(problem_());
        jacAsm_ = Dune::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());

        asImp_().applyInitialSolution_();

        // resize the hint vectors
        if (isBox && enableHints_) {
            int numVertices = gridView_().size(dim);
            curHints_.resize(numVertices);
            prevHints_.resize(numVertices);
            hintsUsable_.resize(numVertices);
            std::fill(hintsUsable_.begin(),
                      hintsUsable_.end(),
                      false);
        }

        // also set the solution of the "previous" time step to the
        // initial solution.
        uPrev_ = uCur_;
    }

    void setHints(const Element &element,
                  ElementVolumeVariables &prevVolVars,
                  ElementVolumeVariables &curVolVars) const
    {
        if (!isBox || !enableHints_)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int n = element.subEntities(dim);
#else
        int n = element.template count<dim>();
#endif
        prevVolVars.resize(n);
        curVolVars.resize(n);
        for (int i = 0; i < n; ++i) {
            int globalIdx = vertexMapper().map(element, i, dim);

            if (!hintsUsable_[globalIdx]) {
                curVolVars[i].setHint(NULL);
                prevVolVars[i].setHint(NULL);
            }
            else {
                curVolVars[i].setHint(&curHints_[globalIdx]);
                prevVolVars[i].setHint(&prevHints_[globalIdx]);
            }
        }
    }

    void setHints(const Element &element,
                  ElementVolumeVariables &curVolVars) const
    {
        if (!isBox || !enableHints_)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int n = element.subEntities(dim);
#else
        int n = element.template count<dim>();
#endif
        curVolVars.resize(n);
        for (int i = 0; i < n; ++i) {
            int globalIdx = vertexMapper().map(element, i, dim);

            if (!hintsUsable_[globalIdx])
                curVolVars[i].setHint(NULL);
            else
                curVolVars[i].setHint(&curHints_[globalIdx]);
        }
    }

    void updatePrevHints()
    {
        if (!isBox || !enableHints_)
            return;

        prevHints_ = curHints_;
    }

    void updateCurHints(const Element &element,
                        const ElementVolumeVariables &elemVolVars) const
    {
        if (!isBox || !enableHints_)
            return;

        for (unsigned int i = 0; i < elemVolVars.size(); ++i) {
            int globalIdx = vertexMapper().map(element, i, dim);
            curHints_[globalIdx] = elemVolVars[i];
            if (!hintsUsable_[globalIdx])
                prevHints_[globalIdx] = elemVolVars[i];
            hintsUsable_[globalIdx] = true;
        }
    }


    /*!
     * \brief Compute the global residual for an arbitrary solution
     *        vector.
     *
     * \param residual Stores the result
     * \param u The solution for which the residual ought to be calculated
     */
    Scalar globalResidual(SolutionVector &residual,
                          const SolutionVector &u)
    {
        SolutionVector tmp(curSol());
        curSol() = u;
        Scalar res = globalResidual(residual);
        curSol() = tmp;
        return res;
    }

    /*!
     * \brief Compute the global residual for the current solution
     *        vector.
     *
     * \param residual Stores the result
     */
    Scalar globalResidual(SolutionVector &residual)
    {
        residual = 0;

        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            localResidual().eval(*eIt);

            if (isBox)
            {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                for (int i = 0; i < eIt->subEntities(dim); ++i) {
#else
                for (int i = 0; i < eIt->template count<dim>(); ++i) {
#endif
                    int globalI = vertexMapper().map(*eIt, i, dim);
                    residual[globalI] += localResidual().residual(i);
                }
            }
            else
            {
                int globalI = elementMapper().map(*eIt);
                residual[globalI] = localResidual().residual(0);
            }
        }

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView_().comm().size() > 1)
            result2 = gridView_().comm().sum(result2);

        // add up the residuals on the process borders
        if (isBox && gridView_().comm().size() > 1) {
            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumHandle(residual, vertexMapper());
            gridView_().communicate(sumHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }

        return std::sqrt(result2);
    }

    /*!
     * \brief Compute the integral over the domain of the storage
     *        terms of all conservation quantities.
     *
     * \param storage Stores the result
     */
    void globalStorage(PrimaryVariables &storage)
    {
        storage = 0;

        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if(eIt->partitionType() == Dune::InteriorEntity)
            {
                localResidual().evalStorage(*eIt);

                if (isBox)
                {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                    for (int i = 0; i < eIt->subEntities(dim); ++i)
#else
                    for (int i = 0; i < eIt->template count<dim>(); ++i)
#endif
                    {
                        storage += localResidual().storageTerm()[i];
                    }
                }
                else
                {
                    storage += localResidual().storageTerm()[0];
                }
            }
        }

        if (gridView_().comm().size() > 1)
            storage = gridView_().comm().sum(storage);
    }

    /*!
     * \brief Returns the volume \f$\mathrm{[m^3]}\f$ of a given control volume.
     *
     * \param globalIdx The global index of the control volume's
     *                  associated vertex
     */
    Scalar boxVolume(const int globalIdx) const
    {
        if (isBox)
        {
            return boxVolume_[globalIdx][0];
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "requested box volume for cell-centered model");
        }
    }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    const SolutionVector &curSol() const
    { return uCur_; }

    /*!
     * \brief Reference to the current solution as a block vector.
     */
    SolutionVector &curSol()
    { return uCur_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    const SolutionVector &prevSol() const
    { return uPrev_; }

    /*!
     * \brief Reference to the previous solution as a block vector.
     */
    SolutionVector &prevSol()
    { return uPrev_; }

    /*!
     * \brief Returns the operator assembler for the global jacobian of
     *        the problem.
     */
    JacobianAssembler &jacobianAssembler()
    { return *jacAsm_; }

    /*!
     * \copydoc jacobianAssembler()
     */
    const JacobianAssembler &jacobianAssembler() const
    { return *jacAsm_; }

    /*!
     * \brief Returns the local jacobian which calculates the local
     *        stiffness matrix for an arbitrary element.
     *
     * The local stiffness matrices of the element are used by
     * the jacobian assembler to produce a global linerization of the
     * problem.
     */
    LocalJacobian &localJacobian()
    { return localJacobian_; }
    /*!
     * \copydoc localJacobian()
     */
    const LocalJacobian &localJacobian() const
    { return localJacobian_; }

    /*!
     * \brief Returns the local residual function.
     */
    LocalResidual &localResidual()
    { return localJacobian().localResidual(); }
    /*!
     * \copydoc localResidual()
     */
    const LocalResidual &localResidual() const
    { return localJacobian().localResidual(); }

    /*!
     * \brief Returns the relative error between two vectors of
     *        primary variables.
     *
     * \param dofIdx The global index of the control volume's
     *               associated degree of freedom
     * \param priVars1 The first vector of primary variables
     * \param priVars2 The second vector of primary variables
     */
    Scalar relativeErrorDof(const int dofIdx,
                            const PrimaryVariables &priVars1,
                            const PrimaryVariables &priVars2)
    {
        Scalar result = 0.0;
        for (int j = 0; j < numEq; ++j) {
            Scalar eqErr = std::abs(priVars1[j] - priVars2[j]);
            eqErr /= std::max<Scalar>(1.0, std::abs(priVars1[j] + priVars2[j])/2);
            
            result = std::max(result, eqErr);
        }
        return result;
    }
    
    /*!
     * \brief Try to progress the model to the next timestep.
     *
     * \param solver The non-linear solver
     * \param controller The controller which specifies the behaviour
     *                   of the non-linear solver
     */
    bool update(NewtonMethod &solver,
                NewtonController &controller)
    {
#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i)
            Valgrind::CheckDefined(curSol()[i]);
#endif // HAVE_VALGRIND

        asImp_().updateBegin();

        int converged = solver.execute(controller);
         
        if (this->gridView_().comm().size() > 1)
        {
            converged = this->gridView_().comm().min(converged);
        }
        if (converged) {
            asImp_().updateSuccessful();
        }
        else
            asImp_().updateFailed();

#if HAVE_VALGRIND
        for (size_t i = 0; i < curSol().size(); ++i) {
            Valgrind::CheckDefined(curSol()[i]);
        }
#endif // HAVE_VALGRIND

        return converged;
    }

    /*!
     * \brief Check the plausibility of the current solution
     *
     *        This has to be done by the actual model, it knows
     *        best, what (ranges of) variables to check.
     *        This is primarily a hook
     *        which the actual model can overload.
     */
    void checkPlausibility() const
    { }

    /*!
     * \brief Called by the update() method before it tries to
     *        apply the newton method. This is primarily a hook
     *        which the actual model can overload.
     */
    void updateBegin()
    { }


    /*!
     * \brief Called by the update() method if it was
     *        successful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateSuccessful()
    { }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primarily a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        uCur_ = uPrev_;
        if (isBox)
            curHints_ = prevHints_;

        jacAsm_->reassembleAll();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and
     *        the result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        // make the current solution the previous one.
        uPrev_ = uCur_;
        if (isBox)
            prevHints_ = curHints_;

        updatePrevHints();
    }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        if (isBox)
            res.template serializeEntities<dim>(asImp_(), this->gridView_()); 
        else
            res.template serializeEntities<0>(asImp_(), this->gridView_());
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        if (isBox)
            res.template deserializeEntities<dim>(asImp_(), this->gridView_());
        else
            res.template deserializeEntities<0>(asImp_(), this->gridView_());
        
        prevSol() = curSol();
    }

    /*!
     * \brief Write the current solution for a vertex to a restart
     *        file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void serializeEntity(std::ostream &outstream,
                         const Entity &entity)
    {
        int dofIdx = dofMapper().map(entity);
        
        // write phase state
        if (!outstream.good()) {
            DUNE_THROW(Dune::IOError,
                       "Could not serialize vertex "
                       << dofIdx);
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            outstream << curSol()[dofIdx][eqIdx] << " ";
        }
    }

    /*!
     * \brief Reads the current solution variables for a vertex from a
     *        restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param entity The entity which's data should be
     *               serialized, i.e. a vertex for the box method
     *               and an element for the cell-centered method
     */
    template <class Entity>
    void deserializeEntity(std::istream &instream,
                           const Entity &entity)
    {
        int dofIdx = dofMapper().map(entity);

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                DUNE_THROW(Dune::IOError,
                           "Could not deserialize vertex "
                           << dofIdx);
            instream >> curSol()[dofIdx][eqIdx];
        }
    }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDofs() const
    {
        if (isBox)
            return gridView_().size(dim); 
        else
            return gridView_().size(0); 
    }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are
     *        defined to indices.
     *
     * Is the box method is used, this means a mapper 
     * for vertices, if the cell centered method is used,
     * this means a mapper for elements.
     */
    template <class T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox), VertexMapper>::type &dofMapper() const
    {
        return problem_().vertexMapper();
    }
    template <class T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox), ElementMapper>::type &dofMapper() const
    {
        return problem_().elementMapper();
    }

    /*!
     * \brief Mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return problem_().vertexMapper(); }

    /*!
     * \brief Mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return problem_().elementMapper(); }

    /*!
     * \brief Resets the Jacobian matrix assembler, so that the
     *        boundary types can be altered.
     */
    void resetJacobianAssembler ()
    {
        jacAsm_.template reset<JacobianAssembler>(0);
        jacAsm_ = Dune::make_shared<JacobianAssembler>();
        jacAsm_->init(problem_());
    }

    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    void updatePVWeights(const Element &element,
                         const ElementVolumeVariables &volVars) const
    { }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param writer  The VTK multi writer object on which the fields should be added.
     * \param u       The solution function
     * \param deltaU  The delta of the solution function before and after the Newton update
     */
    template <class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
                                 const SolutionVector &u,
                                 const SolutionVector &deltaU)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        SolutionVector residual(u);
        asImp_().globalResidual(residual, u);

        // create the required scalar fields
        unsigned numDofs = asImp_().numDofs();

        // global defect of the two auxiliary equations
        ScalarField* def[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            x[eqIdx] = writer.allocateManagedBuffer(numDofs);
            delta[eqIdx] = writer.allocateManagedBuffer(numDofs);
            def[eqIdx] = writer.allocateManagedBuffer(numDofs);
        }

        for (unsigned int globalIdx = 0; globalIdx < u.size(); globalIdx++)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) 
            {
                (*x[eqIdx])[globalIdx] = u[globalIdx][eqIdx];
                (*delta[eqIdx])[globalIdx] = - deltaU[globalIdx][eqIdx];
                (*def[eqIdx])[globalIdx] = residual[globalIdx][eqIdx];
            }
        }
        
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            std::ostringstream oss;
            oss.str(""); oss << "x_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*x[eqIdx], oss.str());
            else
                writer.attachCellData(*x[eqIdx], oss.str());
            oss.str(""); oss << "delta_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*delta[eqIdx], oss.str());
            else
                writer.attachCellData(*delta[eqIdx], oss.str());
            oss.str(""); oss << "defect_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*def[eqIdx], oss.str());
            else
                writer.attachCellData(*def[eqIdx], oss.str());
        }

        asImp_().addOutputVtkFields(u, writer);
    }

    /*!
     * \brief Add the quantities of a time step which ought to be written to disk.
     *
     * This should be overwritten by the actual model if any secondary
     * variables should be written out. Read: This should _always_ be
     * overwritten by well behaved models!
     *
     * \tparam MultiWriter The type of the VTK multi writer
     *
     * \param sol The global vector of primary variable values.
     * \param writer The VTK multi writer where the fields should be added.
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numDofs = asImp_().numDofs();

        // global defect of the two auxiliary equations
        ScalarField* x[numEq];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            x[eqIdx] = writer.allocateManagedBuffer(numDofs);
        }

        for (int globalIdx = 0; globalIdx < sol.size(); globalIdx++)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                (*x[eqIdx])[globalIdx] = sol[globalIdx][eqIdx];
            }
        }

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            std::ostringstream oss;
            oss << "primaryVar_" << eqIdx;
            if (isBox)
                writer.attachVertexData(*x[eqIdx], oss.str());
            else
                writer.attachCellData(*x[eqIdx], oss.str());
        }
    }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView() const
    { return problem_().gridView(); }

    /*!
     * \brief Returns true if the entity indicated by 'globalIdx' 
     * is located on / touches the grid's boundary.
     *
     * \param globalIdx The global index of the entity
     */
    bool onBoundary(const int globalIdx) const
    { return boundaryIndices_[globalIdx]; }

    /*!
     * \brief Returns true if a vertex is located on the grid's
     *        boundary.
     *
     * \param element A DUNE Codim<0> entity which contains the control
     *             volume's associated vertex.
     * \param vIdx The local vertex index inside element
     */
    bool onBoundary(const Element &element, const int vIdx) const
    {
        if (isBox)
            return onBoundary(vertexMapper().map(element, vIdx, dim));
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "requested for cell-centered model");            
    }

    
    /*!
     * \brief Returns true if the control volume touches
     *        the grid's boundary.
     *
     * \param element A DUNE Codim<0> entity coinciding with the control
     *             volume.
     */
    bool onBoundary(const Element &element) const
    {
        if (!isBox)
            return onBoundary(elementMapper().map(element));
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "requested for box model");
    }
    
    /*!
     * \brief Fill the fluid state according to the primary variables. 
     * 
     * Taking the information from the primary variables, 
     * the fluid state is filled with every information that is 
     * necessary to evaluate the model's local residual. 
     * 
     * \param priVars The primary variables of the model.
     * \param problem The problem at hand. 
     * \param element The current element. 
     * \param fvGeometry The finite volume element geometry.
     * \param scvIdx The index of the subcontrol volume. 
     * \param fluidState The fluid state to fill. 
     */
    template <class FluidState>
    static void completeFluidState(const PrimaryVariables& priVars,
                                   const Problem& problem,
                                   const Element& element,
                                   const FVElementGeometry& fvGeometry,
                                   const int scvIdx,
                                   FluidState& fluidState)
    {
        VolumeVariables::completeFluidState(priVars, problem, element,
                                            fvGeometry, scvIdx, fluidState);
    }
protected:
    /*!
     * \brief A reference to the problem on which the model is applied.
     */
    Problem &problem_()
    { return *problemPtr_; }
    /*!
     * \copydoc problem_()
     */
    const Problem &problem_() const
    { return *problemPtr_; }

    /*!
     * \brief Reference to the grid view of the spatial domain.
     */
    const GridView &gridView_() const
    { return problem_().gridView(); }

    /*!
     * \brief Reference to the local residal object
     */
    LocalResidual &localResidual_()
    { return localJacobian_.localResidual(); }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     */
    void applyInitialSolution_()
    {
        // first set the whole domain to zero
        uCur_ = Scalar(0.0);
        boxVolume_ = Scalar(0.0);

        FVElementGeometry fvGeometry;

        // iterate through leaf grid and evaluate initial
        // condition at the center of each sub control volume
        //
        // TODO: the initial condition needs to be unique for
        // each vertex. we should think about the API...
        ElementIterator eIt = gridView_().template begin<0>();
        const ElementIterator &eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            // deal with the current element
            fvGeometry.update(gridView_(), *eIt);

            // loop over all element vertices, i.e. sub control volumes
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; scvIdx++)
            {
                // get the global index of the degree of freedom
                int globalIdx = dofMapper().map(*eIt, scvIdx, dofCodim);

                // let the problem do the dirty work of nailing down
                // the initial solution.
                PrimaryVariables initPriVars;
                Valgrind::SetUndefined(initPriVars);
                problem_().initial(initPriVars,
                                   *eIt,
                                   fvGeometry,
                                   scvIdx);
                Valgrind::CheckDefined(initPriVars);

                if (isBox)
                {
                    // add up the initial values of all sub-control
                    // volumes. If the initial values disagree for
                    // different sub control volumes, the initial value
                    // will be the arithmetic mean.
                    initPriVars *= fvGeometry.subContVol[scvIdx].volume;
                    boxVolume_[globalIdx] += fvGeometry.subContVol[scvIdx].volume;
                }
                
                uCur_[globalIdx] += initPriVars;
                Valgrind::CheckDefined(uCur_[globalIdx]);
            }
        }

        // add up the primary variables and the volumes of the boxes
        // which cross process borders
        if (isBox && gridView_().comm().size() > 1) {
            VertexHandleSum<Dune::FieldVector<Scalar, 1>,
                Dune::BlockVector<Dune::FieldVector<Scalar, 1> >,
                VertexMapper> sumVolumeHandle(boxVolume_, vertexMapper());
            gridView_().communicate(sumVolumeHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);

            VertexHandleSum<PrimaryVariables, SolutionVector, VertexMapper>
                sumPVHandle(uCur_, vertexMapper());
            gridView_().communicate(sumPVHandle,
                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }

        if (isBox)
        {
            // divide all primary variables by the volume of their boxes
            for (unsigned int i = 0; i < uCur_.size(); ++i) {
                uCur_[i] /= boxVolume(i);
            }
        }
    }

    /*!
     * \brief Find all indices of boundary vertices (box) / elements (cell centered).
     */
    void updateBoundaryIndices_()
    {
        boundaryIndices_.resize(numDofs());
        std::fill(boundaryIndices_.begin(), boundaryIndices_.end(), false);

        ElementIterator eIt = gridView_().template begin<0>();
        ElementIterator eEndIt = gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            Dune::GeometryType geomType = eIt->geometry().type();
            const ReferenceElement &refElement = ReferenceElements::general(geomType);

            IntersectionIterator isIt = gridView_().ibegin(*eIt);
            IntersectionIterator isEndIt = gridView_().iend(*eIt);
            for (; isIt != isEndIt; ++isIt) {
                if (isIt->boundary()) {
                    if (isBox)
                    {
                        // add all vertices on the intersection to the set of
                        // boundary vertices
                        int faceIdx = isIt->indexInInside();
                        int numFaceVerts = refElement.size(faceIdx, 1, dim);
                        for (int faceVertIdx = 0;
                             faceVertIdx < numFaceVerts;
                             ++faceVertIdx)
                        {
                            int elemVertIdx = refElement.subEntity(faceIdx,
                                                                   1,
                                                                   faceVertIdx,
                                                                   dim);
                            int globalVertIdx = vertexMapper().map(*eIt, elemVertIdx, dim);
                            boundaryIndices_[globalVertIdx] = true;
                        }
                    }
                    else 
                    {
                        int globalIdx = elementMapper().map(*eIt);
                        boundaryIndices_[globalIdx] = true;
                    }
                }
            }
        }
    }

    // the hint cache for the previous and the current volume
    // variables
    mutable std::vector<bool> hintsUsable_;
    mutable std::vector<VolumeVariables> curHints_;
    mutable std::vector<VolumeVariables> prevHints_;

    // the problem we want to solve. defines the constitutive
    // relations, matxerial laws, etc.
    Problem *problemPtr_;

    // calculates the local jacobian matrix for a given element
    LocalJacobian localJacobian_;
    // Linearizes the problem at the current time step using the
    // local jacobian
    Dune::shared_ptr<JacobianAssembler> jacAsm_;

    // the set of all indices of vertices on the boundary
    std::vector<bool> boundaryIndices_;

    // cur is the current iterative solution, prev the converged
    // solution of the previous time step
    SolutionVector uCur_;
    SolutionVector uPrev_;

    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > boxVolume_;

private:
    /*!
     * \brief Returns whether messages should be printed
     */
    bool verbose_() const
    { return gridView_().comm().rank() == 0; }

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    bool enableHints_;
};
} // end namespace Dumux

#include "implicitpropertydefaults.hh"

#endif
