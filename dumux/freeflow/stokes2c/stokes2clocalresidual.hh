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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the compositional Stokes box model.
 *
 */
#ifndef DUMUX_STOKES2C_LOCAL_RESIDUAL_HH
#define DUMUX_STOKES2C_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokes/stokeslocalresidual.hh>

#include <dumux/freeflow/stokes2c/stokes2cvolumevariables.hh>
#include <dumux/freeflow/stokes2c/stokes2cfluxvariables.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the compositional Stokes box model. This is derived
 *        from the Stokes box model.
 */
template<class TypeTag>
class Stokes2cLocalResidual : public StokesLocalResidual<TypeTag>
{
    typedef StokesLocalResidual<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { transportEqIdx = Indices::transportEqIdx }; //!< Index of the transport equation
    enum { phaseIdx = Indices::phaseIdx }; //!< Index of the considered phase (only of interest when using two-phase fluidsystems)

    // component indices
    enum { phaseCompIdx = Indices::phaseCompIdx };          //!< Index of the main component of the fluid phase
    enum { transportCompIdx = Indices::transportCompIdx };  //!< Index of the minor component of the fluid phase

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief Evaluate the stored amount of quantities additional to the Stokes model
     *        (transport equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param storage The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage, const int scvIdx, const bool usePrevSol) const
    {
        // compute the storage term for the transport equation
        ParentType::computeStorage(storage, scvIdx, usePrevSol);

        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const ElementVolumeVariables &elemVolVars = usePrevSol ? this->prevVolVars_() : this->curVolVars_();
        const VolumeVariables &volVars = elemVolVars[scvIdx];

        // compute the storage of the component
        storage[transportEqIdx] =
            volVars.density() *
            volVars.fluidState().massFraction(phaseIdx, transportCompIdx);

        Valgrind::CheckDefined(volVars.density());
        Valgrind::CheckDefined(volVars.fluidState().massFraction(phaseIdx, transportCompIdx));
    }

    /*!
     * \brief Evaluates the advective component (mass) flux
     * over a face of a sub-control volume and writes the result in
     * the flux vector.
     *
     * This method is called by compute flux (base class).
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // call computation of the advective fluxes of the stokes model
        // (momentum and mass fluxes)
        ParentType::computeAdvectiveFlux(flux, fluxVars);

        // vertex data of the upstream and the downstream vertices
        const VolumeVariables &up = this->curVolVars_(fluxVars.upstreamIdx());
        const VolumeVariables &dn = this->curVolVars_(fluxVars.downstreamIdx());

        Scalar tmp = fluxVars.normalVelocity();

        if (this->massUpwindWeight_ > 0.0)
            tmp *=  this->massUpwindWeight_ *         // upwind data
                up.density() * up.fluidState().massFraction(phaseIdx, transportCompIdx);
        if (this->massUpwindWeight_ < 1.0)
            tmp += (1.0 - this->massUpwindWeight_) *     // rest
                dn.density() * dn.fluidState().massFraction(phaseIdx, transportCompIdx);

        flux[transportEqIdx] += tmp;
        Valgrind::CheckDefined(flux[transportEqIdx]);
    }

    /*!
     * \brief Adds the diffusive component flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the SCV face or boundary face
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxVars);

        // diffusive component flux
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            flux[transportEqIdx] -=
                fluxVars.moleFractionGrad()[dimIdx] *
                fluxVars.face().normal[dimIdx] *
                (fluxVars.diffusionCoeff() + fluxVars.eddyDiffusivity()) *
                fluxVars.molarDensity() *
                FluidSystem::molarMass(transportCompIdx);

        Valgrind::CheckDefined(flux[transportEqIdx]);
    }
};

}

#endif
