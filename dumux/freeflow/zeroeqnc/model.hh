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
 * \brief Base class for all models which use the compositional ZeroEq box model.
 */
#ifndef DUMUX_ZEROEQNC_MODEL_HH
#define DUMUX_ZEROEQNC_MODEL_HH

#include "indices.hh"
#include "properties.hh"
#include <dumux/freeflow/stokesnc/model.hh>
#include <dumux/freeflow/zeroeq/model.hh>

namespace Dumux
{
/*!
 * \ingroup BoxZeroEqncModel
 * \brief Adaption of the box scheme to the compositional ZeroEq model.
 *
 * This model implements an single-phase isothermal compositional free flow
 * solving the mass and the momentum balance. For the momentum balance
 * the Reynolds-averaged Navier-Stokes (RANS) equation with zero equation
 * (algebraic) turbulence model is used.
 *
 * Mass balance:
 * \f[
 *  \frac{\partial \varrho_\textrm{g}}{\partial t}
 *  + \text{div} \left( \varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g} \right)
 *  - q_\textrm{g} = 0
 * \f]
 *
 * Momentum Balance:
 * \f[
 *   \frac{\partial \left(\varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g}\right)}{\partial t}
 *   + \text{div} \left(
 *     \varrho_\textrm{g} {\boldsymbol{v}_\textrm{g} {\boldsymbol{v}}_\textrm{g}}
 *     - \left[ \mu_\textrm{g} + \mu_\textrm{g,t} \right]
 *       \left( \textbf{grad}\, \boldsymbol{v}_\textrm{g}
 *              + \textbf{grad}\, \boldsymbol{v}_\textrm{g}^T \right)
 *   \right)
 *   + \left(p_\textrm{g} {\bf {I}} \right)
 *   - \varrho_\textrm{g} {\bf g} = 0
 * \f]
 *
 * Component mass balance equations:
 * \f[
 *  \frac{\partial \left(\varrho_\textrm{g} X_\textrm{g}^\kappa\right)}{\partial t}
 *  + \text{div} \left( \varrho_\textrm{g} {\boldsymbol{v}}_\textrm{g} X_\textrm{g}^\kappa
 *  - \left[ D^\kappa_\textrm{g} + D^\kappa_\textrm{g,t} \right]
 *    \varrho_\textrm{g} \frac{M^\kappa}{M_\textrm{g}} \textbf{grad}\, x_\textrm{g}^\kappa \right)
 *  - q_\textrm{g}^\kappa = 0
 * \f]
 * Please note that, even though it is n-component model, the diffusive
 * fluxes are still calculated with binary diffusion.
 *
 * This is discretized by a fully-coupled vertex-centered finite volume
 * (box) scheme in space and by the implicit Euler method in time.
 */
template<class TypeTag>
class ZeroEqncModel : public ZeroEqModel<TypeTag>
{
    typedef ZeroEqModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        intervals = GET_PROP_VALUE(TypeTag, NumberOfIntervals),
        walls = (GET_PROP_VALUE(TypeTag, BBoxMinIsWall) ? 1 : 0)
                + (GET_PROP_VALUE(TypeTag, BBoxMaxIsWall) ? 1 : 0)
    };
    enum { transportCompIdx = Indices::transportCompIdx,
           numComponents = Indices::numComponents };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx) };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;


public:
    ZeroEqncModel()
        :   flowNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal))
        ,   wallNormal_(GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal))
    {
        eps_ = 1e-6;

        // check whether sand grain roughness may be used
        if (GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMinSandGrainRoughness) > 0
            || GET_PARAM_FROM_GROUP(TypeTag, Scalar, ZeroEq, BBoxMaxSandGrainRoughness) > 0)
        {
            Dune::dwarn << "warning: surface roughness will not be used for eddy diffusivity models." << std::endl;
        }
    }

    //! \copydoc ImplicitModel::addOutputVtkFields
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > VectorField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_().size(dim);
        unsigned numElements = this->gridView_().size(0);
        ScalarField &pN = *writer.allocateManagedBuffer(numVertices);
        ScalarField &delP = *writer.allocateManagedBuffer(numVertices);
        ScalarField *moleFraction[numComponents];
        ScalarField *massFraction[numComponents];
        ScalarField &rho = *writer.allocateManagedBuffer(numVertices);
        VectorField &velocity = *writer.template allocateManagedBuffer<Scalar, dim> (numVertices);
        ScalarField &mu = *writer.allocateManagedBuffer(numVertices);
        ScalarField &D = *writer.allocateManagedBuffer(numVertices);
        ScalarField &mut = *writer.allocateManagedBuffer(numElements);
        ScalarField &nut = *writer.allocateManagedBuffer(numElements);
        ScalarField &lmix = *writer.allocateManagedBuffer(numElements);
        ScalarField &uPlus = *writer.allocateManagedBuffer(numElements);
        ScalarField &yPlus = *writer.allocateManagedBuffer(numElements);
        ScalarField &Dt = *writer.allocateManagedBuffer(numElements);
        ScalarField &rank = *writer.allocateManagedBuffer(numElements);
        for (int i = 0; i < numComponents; ++i)
            moleFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);
        for (int i = 0; i < numComponents; ++i)
            massFraction[i] = writer.template allocateManagedBuffer<Scalar, 1>(numVertices);

        // write volume values to .vtu and .csv
        std::ofstream volVarsFile("volVarsData.csv", std::ios_base::out);
        asImp_().writeVolVarsHeader(volVarsFile);
        volVarsFile << std::endl;

        FVElementGeometry fvGeometry;
        ElementBoundaryTypes elemBcTypes;
        VolumeVariables volVars;

        for (const auto& element : elements(this->gridView_()))
        {
            int idx = this->elementMapper().index(element);
            rank[idx] = this->gridView_().comm().rank();

            fvGeometry.update(this->gridView_(), element);
            elemBcTypes.update(this->problem_(), element, fvGeometry);

            int numLocalVerts = element.template subEntities(dim);
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int vIdxGlobal = this->vertexMapper().subIndex(element, i, dim);
                volVars.update(sol[vIdxGlobal],
                               this->problem_(),
                               element,
                               fvGeometry,
                               i,
                               false);

                pN[vIdxGlobal] = volVars.pressure();
                delP[vIdxGlobal] = volVars.pressure() - 1e5;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                {
                    (*moleFraction[compIdx])[vIdxGlobal]= volVars.moleFraction(compIdx);
                    (*massFraction[compIdx])[vIdxGlobal]= volVars.massFraction(compIdx);
                    Valgrind::CheckDefined((*moleFraction[compIdx])[vIdxGlobal]);
                    Valgrind::CheckDefined((*massFraction[compIdx])[vIdxGlobal]);
                }
                rho[vIdxGlobal] = volVars.density();
                mu[vIdxGlobal] = volVars.dynamicViscosity();
                D[vIdxGlobal] = volVars.diffusionCoeff(transportCompIdx);
                velocity[vIdxGlobal] = volVars.velocity();

                asImp_().writeVolVarsData(volVarsFile, volVars);
                volVarsFile << std::endl;
            }
        }
        writer.attachVertexData(pN, "P");
        writer.attachVertexData(delP, "delP");

        for (int j = 0; j < numComponents; ++j)
        {
            std::ostringstream moleFrac, massFrac;
            moleFrac << "x_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*moleFraction[j], moleFrac.str().c_str());

            massFrac << "X_" << FluidSystem::phaseName(phaseIdx)
                     << "^" << FluidSystem::componentName(j);
            writer.attachVertexData(*massFraction[j], massFrac.str().c_str());
        }
        volVarsFile.close();

        writer.attachVertexData(rho, "rho");
        writer.attachVertexData(mu, "mu");
        writer.attachVertexData(D, "D");
        writer.attachVertexData(velocity, "v", dim);

        // ensure that the actual values are given out
        asImp_().updateWallProperties();

        // write flux values to .vtu and .csv
        std::ofstream fluxFile("fluxVarsData.csv", std::ios_base::out);
        asImp_().writeFluxVarsHeader(fluxFile);
        fluxFile << std::endl;

        for (const auto& element : elements(this->gridView_()))
        {
            fvGeometry.update(this->gridView_(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false);

            unsigned int numFluxVars = 0;
            Scalar sumDynamicEddyViscosity = 0.0;
            Scalar sumKinematicEddyViscosity = 0.0;
            Scalar sumMixingLength = 0.0;
            Scalar sumUPlus = 0.0;
            Scalar sumYPlus = 0.0;
            Scalar sumEddyDiffusivity = 0.0;

            for (const auto& intersection : intersections(this->gridView_(), element))
            {
                int fIdx = intersection.indexInInside();

                FluxVariables fluxVars;
                fluxVars.update(this->problem_(),
                                element,
                                fvGeometry,
                                fIdx,
                                elemVolVars,
                                false);

                asImp_().writeFluxVarsData(fluxFile, fluxVars);
                fluxFile << std::endl;

                sumDynamicEddyViscosity += fluxVars.dynamicEddyViscosity();
                sumKinematicEddyViscosity += fluxVars.kinematicEddyViscosity();
                sumMixingLength += fluxVars.mixingLength();
                sumUPlus += fluxVars.uPlus();
                sumYPlus += fluxVars.yPlusRough();
                sumEddyDiffusivity += fluxVars.eddyDiffusivity();
                numFluxVars += 1;
            }

            int eIdxGlobal = this->elementMapper().index(element);
            mut[eIdxGlobal] = sumDynamicEddyViscosity / numFluxVars;
            nut[eIdxGlobal] = sumKinematicEddyViscosity / numFluxVars;
            lmix[eIdxGlobal] = sumMixingLength / numFluxVars;
            uPlus[eIdxGlobal] = sumUPlus / numFluxVars;
            yPlus[eIdxGlobal] = sumYPlus / numFluxVars;
            Dt[eIdxGlobal] = sumEddyDiffusivity / numFluxVars;
        }
        fluxFile.close();

        writer.attachCellData(mut, "mu_t");
        writer.attachCellData(nut, "nu_t");
        writer.attachCellData(lmix, "l_mix");
        writer.attachCellData(uPlus, "u^+");
        writer.attachCellData(yPlus, "y^+");
        writer.attachCellData(Dt, "D_t");
    }

    //! \copydoc ZeroEqModel::writeFluxVarsHeader
    void writeFluxVarsHeader(std::ofstream &stream)
    {
        ParentType::writeFluxVarsHeader(stream);
        stream << "," << "eddyDiffusivity";
    }

    //! \copydoc ZeroEqModel::writeFluxVarsData
    void writeFluxVarsData(std::ofstream &stream, const FluxVariables &fluxVars)
    {
        ParentType::writeFluxVarsData(stream, fluxVars);
        stream << "," << fluxVars.eddyDiffusivity();
    }

    /*!
     * \name Wall properties
     */
    // \{

    //! \copydoc ZeroEqModel::resetWallProperties
    void resetWallProperties()
    {
        ParentType::resetWallProperties();

        for (int wallIdx = 0; wallIdx < walls; ++wallIdx)
            for (int posIdx = 0; posIdx < intervals; ++posIdx)
            {
                this->wall[wallIdx].maxMassFraction[posIdx] = 0.0;
                this->wall[wallIdx].maxMoleFraction[posIdx] = 0.0;
            }
    }

    //! \copydoc ZeroEqModel::calculateMaxFluxVars
    void calculateMaxFluxVars(const FluxVariables &fluxVars, const GlobalPosition &globalPos)
    {
        ParentType::calculateMaxFluxVars(fluxVars, globalPos);

        int posIdx = this->getPosIdx(globalPos);
        int wallIdx = this->getWallIdx(globalPos, posIdx);
        // mass fraction
        if (this->wall[wallIdx].maxMassFraction[posIdx] < fluxVars.massFraction(transportCompIdx))
            for (int wIdx = 0; wIdx < walls; ++wIdx)
                this->wall[wIdx].maxMassFraction[posIdx] = fluxVars.massFraction(transportCompIdx);
        // mole fraction
        if (this->wall[wallIdx].maxMoleFraction[posIdx] < fluxVars.moleFraction(transportCompIdx))
            for (int wIdx = 0; wIdx < walls; ++wIdx)
                this->wall[wIdx].maxMoleFraction[posIdx] = fluxVars.moleFraction(transportCompIdx);
    }

    //! \copydoc ZeroEqModel::doInterpolationFluxValues
    const void doInterpolationFluxValues(const int wallIdx, const int posIdx, const int prevIdx, const int nextIdx)
    {
        ParentType::doInterpolationFluxValues(wallIdx, posIdx, prevIdx, nextIdx);

        this->wall[wallIdx].maxMassFraction[posIdx] = this->interpolation(posIdx, prevIdx, this->wall[wallIdx].maxMassFraction[prevIdx], nextIdx, this->wall[wallIdx].maxMassFraction[nextIdx]);
        this->wall[wallIdx].maxMoleFraction[posIdx] = this->interpolation(posIdx, prevIdx, this->wall[wallIdx].maxMoleFraction[prevIdx], nextIdx, this->wall[wallIdx].maxMoleFraction[nextIdx]);
    }

    // \} // wall properties

    /*!
     * \brief Returns the name of used the eddy diffusivity model.
     */
    std::string eddyDiffusivityModelName() const
    {
        switch (GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyDiffusivityModel))
        {
            case EddyDiffusivityIndices::noEddyDiffusivityModel: // 0
                return "noEddyDiffusivityModel";
                break;
            case EddyDiffusivityIndices::reynoldsAnalogy: // 1
                return "reynoldsAnalogy";
                break;
            case EddyDiffusivityIndices::modifiedVanDriest: // 2
                return "modifiedVanDriest";
                break;
            case EddyDiffusivityIndices::deissler: // 3
                return "deissler";
                break;
            case EddyDiffusivityIndices::meier: // 4
                return "meier";
                break;
            case EddyDiffusivityIndices::exponential: // 5
                return "exponential";
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "This eddy diffusivity model is not implemented.");
        }
    }

private:
    Scalar eps_;
    const int flowNormal_;
    const int wallNormal_;

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#include "propertydefaults.hh"

#endif // DUMUX_ZEROEQNC_MODEL_HH
