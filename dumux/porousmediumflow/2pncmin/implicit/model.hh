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
* \brief Adaption of the fully implicit box scheme to the two-phase n-component flow model.
*/

#ifndef DUMUX_2PNCMIN_MODEL_HH
#define DUMUX_2PNCMIN_MODEL_HH

#include "properties.hh"
#include "indices.hh"
#include "primaryvariableswitch.hh"
#include "localresidual.hh"

#include <dumux/porousmediumflow/2pnc/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCMinModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model with additional solid/mineral phases.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme (this is not done for 2pnc approach yet, however possible) as
 * spatial and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 *Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 * The primary variable of the solid phases is the volume fraction \f$\phi_\lambda = \frac{V_\lambda}{V_{total}}\f$.
 *
 * The source an sink terms link the mass balances of the n-transported component to the solid phases.
 * The porosity \f$\phi\f$ is updated according to the reduction of the initial (or solid-phase-free porous medium) porosity \f$\phi_0\f$
 * by the accumulated volume fractions of the solid phases:
 * \f$ \phi = \phi_0 - \sum (\phi_\lambda)\f$
 * Additionally, the permeability is updated depending on the current porosity.
 */

template<class TypeTag>
class TwoPNCMinModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    // the parent class needs to access the variable switch
    friend typename GET_PROP_TYPE(TypeTag, BaseModel);

    using ThisType = Dumux::TwoPNCMinModel<TypeTag>;
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseModel);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numSPhases = GET_PROP_VALUE(TypeTag, NumSPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using CoordScalar = typename GridView::ctype;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:

    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param problem The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem& problem)
    {
        ParentType::init(problem);

        // register standardized vtk output fields
        auto& vtkOutputModule = problem.vtkOutputModule();
        vtkOutputModule.addSecondaryVariable("Sn", [](const VolumeVariables& v){ return v.saturation(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("Sw", [](const VolumeVariables& v){ return v.saturation(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pn", [](const VolumeVariables& v){ return v.pressure(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pw", [](const VolumeVariables& v){ return v.pressure(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("pc", [](const VolumeVariables& v){ return v.capillaryPressure(); });
        vtkOutputModule.addSecondaryVariable("rhoW", [](const VolumeVariables& v){ return v.density(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("rhoN", [](const VolumeVariables& v){ return v.density(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("mobW", [](const VolumeVariables& v){ return v.mobility(wPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("mobN", [](const VolumeVariables& v){ return v.mobility(nPhaseIdx); });
        vtkOutputModule.addSecondaryVariable("porosity", [](const VolumeVariables& v){ return v.porosity(); });
        vtkOutputModule.addSecondaryVariable("temperature", [](const VolumeVariables& v){ return v.temperature(); });

        for (int sPhaseIdx = 0; sPhaseIdx < numSPhases; ++sPhaseIdx)
            vtkOutputModule.addSecondaryVariable("precipitateVolumeFraction_" + FluidSystem::phaseName(numPhases + sPhaseIdx),
                                                 [sPhaseIdx](const VolumeVariables& v)
                                                 { return v.precipitateVolumeFraction(numPhases + sPhaseIdx); });

        vtkOutputModule.addSecondaryVariable("Kxx",
                                             [this](const VolumeVariables& v){ return this->perm_(v.permeability())[0][0]; });
        if (dim >= 2)
            vtkOutputModule.addSecondaryVariable("Kyy",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[1][1]; });
        if (dim >= 3)
            vtkOutputModule.addSecondaryVariable("Kzz",
                                                 [this](const VolumeVariables& v){ return this->perm_(v.permeability())[2][2]; });

        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                vtkOutputModule.addSecondaryVariable("x" + FluidSystem::componentName(j) + FluidSystem::phaseName(i),
                                                     [i,j](const VolumeVariables& v){ return v.moleFraction(i,j); });

        for (int j = 0; j < numComponents; ++j)
            vtkOutputModule.addSecondaryVariable("m^w_" + FluidSystem::componentName(j),
                                                 [j](const VolumeVariables& v){ return v.molarity(wPhaseIdx,j); });
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& phasePresence = outputModule.createScalarField("phase presence", dofCodim);
        for (std::size_t i = 0; i < phasePresence.size(); ++i)
            phasePresence[i] = priVarSwitch().phasePresence(i);
    }

    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // \todo resize volvars vector if grid was adapted

        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());

        // update the secondary variables if global caching is enabled
        // \note we only updated if phase presence changed as the volume variables
        //       are already updated once by the switch
        for (const auto& element : elements(this->problem_().gridView()))
        {
            // make sure FVElementGeometry & vol vars are bound to the element
            auto fvGeometry = localView(this->globalFvGeometry());
            fvGeometry.bindElement(element);

            if (switchFlag_)
            {
                for (auto&& scv : scvs(fvGeometry))
                {
                    auto dofIdxGlobal = scv.dofIndex();
                    if (priVarSwitch_().wasSwitched(dofIdxGlobal))
                    {
                        this->nonConstCurGlobalVolVars().volVars(scv.index()).update(this->curSol()[dofIdxGlobal],
                                                                                     this->problem_(),
                                                                                     element,
                                                                                     scv);
                    }
                }
            }

            // handle the boundary volume variables
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // if we are not on a boundary, skip the rest
                if (!scvf.boundary())
                    continue;

                // check if boundary is a pure dirichlet boundary
                const auto bcTypes = this->problem_().boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    const auto insideScvIdx = scvf.insideScvIdx();
                    const auto& insideScv = fvGeometry.scv(insideScvIdx);
                    const auto dirichletPriVars = this->problem_().dirichlet(element, scvf);

                    this->nonConstCurGlobalVolVars().volVars(scvf.outsideScvIdx()).update(dirichletPriVars, this->problem_(), element, insideScv);
                }
            }
        }
    }

    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    template<typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalVolumeVariablesCache), void>::type
    newtonEndStep()
    {
        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol());
    }

    /*!
     * \brief Called by the update() method if applying the Newton
     *        method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        switchFlag_ = false;
        priVarSwitch_().resetPhasePresence();
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();

        // update the phase state
        priVarSwitch_().updateOldPhasePresence();
        switchFlag_ = false;
    }

    /*!
     * \brief Returns true if the primary variables were switched for
     *        at least one dof after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one entity for the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);

        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << dofIdxGlobal);

        outStream << priVarSwitch().phasePresence(dofIdxGlobal) << " ";
    }

    /*!
     * \brief Reads the current solution from a restart file.
     *
     * \param inStream The input stream of one entity from the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);

        // read phase presence
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good())
            DUNE_THROW(Dune::IOError, "Could not deserialize entity " << dofIdxGlobal);

        int phasePresence;
        inStream >> phasePresence;

        priVarSwitch_().setPhasePresence(dofIdxGlobal, phasePresence);
        priVarSwitch_().setOldPhasePresence(dofIdxGlobal, phasePresence);
    }

    const Dumux::TwoPNCMinPrimaryVariableSwitch<TypeTag>& priVarSwitch() const
    { return switch_; }

private:

    Dumux::TwoPNCMinPrimaryVariableSwitch<TypeTag>& priVarSwitch_()
    { return switch_; }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        ParentType::applyInitialSolution_();

        // initialize the primary variable switch
        priVarSwitch_().init(this->problem_());
    }

    //! the class handling the primary variable switch
    Dumux::TwoPNCMinPrimaryVariableSwitch<TypeTag> switch_;
    bool switchFlag_;

    Tensor perm_(Scalar perm) const
    {
        Tensor K(0.0);

        for(int i=0; i<dim; i++)
            K[i][i] = perm;

       return K;
    }

    const Tensor& perm_(const Tensor& perm) const
    {
       return perm;
    }
};

} // end namespace Dumux

#include "propertydefaults.hh"

#endif
