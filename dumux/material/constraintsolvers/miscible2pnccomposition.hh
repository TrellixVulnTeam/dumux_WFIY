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
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 */
#ifndef DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH
#define DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1. After calling the solve() method the following quantities
 * are calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molarities of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
 * - if the setInternalEnergy parameter is true, also specific enthalpies and internal energies of *all* phases
 */
template <class Scalar, class FluidSystem>
class miscible2pncComposition
{
    static const int numPhases = FluidSystem::numPhases;
    static const int numComponents = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numPhases;


public:
    /*!
     * \brief Computes the composition of all phases of a M-phase,
     *        N-component fluid system assuming that all M phases are
     *        present
     *
     * It constraint solver assumes the following quantities to be set:
     *
     * - temperatures of *all* phases
     * - saturations of *all* phases
     * - pressures of *all* phases
     *
     * It also assumes that the mole/mass fractions of all phases sum up
     * to 1. After calling the solve() method the following quantities
     * are calculated in addition:
     *
     * - temperature of *all* phases
     * - density, molar density, molar volume of *all* phases
     * - composition in mole and mass fractions and molarities of *all* phases
     * - mean molar masses of *all* phases
     * - fugacity coefficients of *all* components in *all* phases
     * - if the setViscosity parameter is true, also dynamic viscosities of *all* phases
     * - if the setInternalEnergy parameter is true, also specific enthalpies and internal energies of *all* phases
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      const int knownPhaseIdx,
                      bool setViscosity,
                      bool setInternalEnergy)
    {
#ifndef NDEBUG
        // currently this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids. TODO: relax this
        // (requires solving a non-linear system of equations, i.e. using
        // newton method.)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));
        }
#endif

        //get the known mole fractions from the fluidState
        Dune::FieldVector<Scalar, numComponents-numMajorComponents> xKnown(0.0);
        for (int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
        	xKnown[knownCompIdx] = fluidState.moleFraction(knownPhaseIdx, knownCompIdx + numMajorComponents);
        }

        // compute all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
            }
        }


        // create the linear system of equations which defines the
        // mole fractions
        Dune::FieldMatrix<Scalar, numComponents*numPhases, numComponents*numPhases> M(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> x(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> b(0.0);

        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int rowIdx = numComponents*(numPhases - 1) + phaseIdx;
            b[rowIdx] = 1.0;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int colIdx = phaseIdx*numComponents + compIdx;

                M[rowIdx][colIdx] = 1.0;
            }
        }
        //set the additional equations for the numComponents-numMajorComponents
        //Components, of which the molefractions are known, set to molefraction(knownCompIdx)=xKnown
        for(int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
        	int rowIdx = numComponents + numPhases + knownCompIdx;
        	int colIdx = knownPhaseIdx*numComponents + knownCompIdx + numMajorComponents;
        	M[rowIdx][colIdx] = 1.0;
        	b[rowIdx] = xKnown[knownCompIdx];
        }

        // assemble the equations expressing the fact that the
        // fugacities of each component is equal in all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            Scalar entryCol1 =
                fluidState.fugacityCoefficient(/*phaseIdx=*/0, compIdx)
                * fluidState.pressure(/*phaseIdx=*/0);
            int col1Idx = compIdx;

            for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
                int rowIdx = (phaseIdx - 1)*numComponents + compIdx;
                int col2Idx = phaseIdx*numComponents + compIdx;

                Scalar entryCol2 =
                    fluidState.fugacityCoefficient(phaseIdx, compIdx)
                    * fluidState.pressure(phaseIdx);

                M[rowIdx][col1Idx] = entryCol1;
                M[rowIdx][col2Idx] = -entryCol2;
            }
        }

        //preconditioning of M to reduce condition number
        //prevents matrix meeting dune's singularity criteria
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
        	if (compIdx < numMajorComponents)
        	{
        		for (int colIdx = 0; colIdx < numPhases*numComponents; colIdx++)
        		{
        			//Multiply row of main component (Raoult's Law) with 10e-5 (order of magn. of pressure)
        			M[compIdx][colIdx] *= 10e-5;
        		}
        	} else {
        		for (int colIdx = 0; colIdx < numPhases*numComponents; colIdx++)
        		{
        			//Multiply row of sec. components (Henry's Law) with 10e-9 (order of magn. of Henry constant)
        			M[compIdx][colIdx] *= 10e-9;
        		}
        	}
        }

        // solve for all mole fractions
        M.solve(x, b);
        
        // set all mole fractions and the the additional quantities in
        // the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int rowIdx = phaseIdx*numComponents + compIdx;
                fluidState.setMoleFraction(phaseIdx, compIdx, x[rowIdx]);
            }
            paramCache.updateComposition(fluidState, phaseIdx);

            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            if (setViscosity) {
                value = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                fluidState.setViscosity(phaseIdx, value);
            }

            if (setInternalEnergy) {
                value = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                fluidState.setEnthalpy(phaseIdx, value);
            }
        }
    }
};

} // end namespace Dumux

#endif
