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
 *        using the non-isothermal compositional n-compontent Stokes box model.
 *
 */
#ifndef DUMUX_STOKESNCNI_LOCAL_RESIDUAL_HH
#define DUMUX_STOKESNCNI_LOCAL_RESIDUAL_HH

#include <dumux/freeflow/stokesnc/stokesnclocalresidual.hh>

#include "stokesncnivolumevariables.hh"
#include "stokesncnifluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup BoxStokesncniModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal compositional n-component Stokes box model. This class is derived
 *        from the stokesnc box local residual and adds the energy balance equation.
 */
template<class TypeTag>
class StokesncniLocalResidual : public StokesncLocalResidual<TypeTag>
{
    typedef StokesncLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension }; // dimensions
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) }; // number of equations
    enum { numComponents = Indices::numComponents }; // number of components
    enum { energyEqIdx = Indices::energyEqIdx}; // equation indices

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief Evaluate the amount the additional quantities to the stokesnc model
     *        (energy equation).
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
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

        // compute the storage of energy
        storage[energyEqIdx] =
            volVars.density() *
            volVars.internalEnergy();
    }

    /*!
     * \brief Evaluates the convective energy flux
     *        over a face of a sub-control volume and writes the result in
     *        the flux vector. This method is called by computeFlux in the base class.
     *
     * \param flux The vector for the fluxes over the SCV/boundary face
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

        tmp *=  this->massUpwindWeight_ *         // upwind data
            up.density() * up.enthalpy() +
            (1.0 - this->massUpwindWeight_) *     // rest
            dn.density() * dn.enthalpy();

        flux[energyEqIdx] += tmp;
        Valgrind::CheckDefined(flux[energyEqIdx]);
    }

    /*!
     * \brief Evaluates the diffusive component energy flux
     *        over the face of a sub-control volume.
     *
     * \param flux The vector for the fluxes over the SCV face
     * \param fluxVars The flux variables at the current SCV face
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const FluxVariables &fluxVars) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, fluxVars);

        // conductive energy flux
        computeConductiveFlux(flux, fluxVars);

        // diffusive component energy flux
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            for (int compIdx=0; compIdx<numComponents; compIdx++)
            {
                flux[energyEqIdx] -= fluxVars.moleFractionGrad(compIdx)[dimIdx]
                                     * fluxVars.face().normal[dimIdx]
                                     *(fluxVars.diffusionCoeff(compIdx) + fluxVars.eddyDiffusivity())
                                     * fluxVars.molarDensity()
                                     * FluidSystem::molarMass(compIdx) // Multiplied by molarMass [kg/mol] to convert from [mol/m^3 s] to [kg/m^3 s]
                                     * fluxVars.componentEnthalpy(compIdx);
            }
        }
        Valgrind::CheckDefined(flux[energyEqIdx]);
    }

    /*!
     * \brief Evaluates the conductive energy flux
     *        over the face of a sub-control volume.
     *
     * \param flux The vector for the fluxes over the SCV face
     * \param fluxVars The flux variables at the current SCV face
     */
    void computeConductiveFlux(PrimaryVariables &flux,
                               const FluxVariables &fluxVars) const
    {
        // diffusive heat flux
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            flux[energyEqIdx] -=
                fluxVars.temperatureGrad()[dimIdx] *
                fluxVars.face().normal[dimIdx] *
                (fluxVars.thermalConductivity() + fluxVars.thermalEddyConductivity());
        }
        Valgrind::CheckDefined(flux[energyEqIdx]);
    }
};

} // end namespace

#endif
