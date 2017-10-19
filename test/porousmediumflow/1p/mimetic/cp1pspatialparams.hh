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
 * \brief The spatial parameters for the LensProblem which uses the
 *        two-phase fully implicit model
 */
#ifndef DUMUX_MIMETIC_CPOneP_SPATIAL_PARAMS_HH
#define DUMUX_MIMETIC_CPOneP_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class CPOnePSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(CPOnePSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(CPOnePSpatialParams, SpatialParams, CPOnePSpatialParams<TypeTag>);

}
/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the CPTwoPProblem which uses the
 *        two-phase fully implicit model
 */
template<class TypeTag>
class CPOnePSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

public:
    // export permeability type
    using PermeabilityType = Tensor;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    CPOnePSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        homogeneous_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Homogeneous);

        const std::vector<int>& globalCell = GridCreator::grid().globalCell();

        if (GridCreator::deck().hasKeyword("PORO")) {
            std::cout << "Found PORO..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PORO").getRawDoubleData();
            porosity_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                if (homogeneous_)
                    porosity_[i] = 0.2;
                else
                    porosity_[i] = eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck().hasKeyword("PERMX")) {
            std::cout << "Found PERMX..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PERMX").getRawDoubleData();
            permX_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permX_[i] = 9.86923e-16;
                else
                    permX_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck().hasKeyword("PERMZ")) {
            std::cout << "Found PERMZ..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PERMZ").getRawDoubleData();
            permZ_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permZ_[i] = 9.86923e-16;
                else
                    permZ_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }
        else {
            std::cout << "PERMZ not found, set equal to PERMX..." << std::endl;
            permZ_ = permX_;
        }
    }

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param globalPos The global position
     */
    PermeabilityType  permeability(const Element& element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        Tensor K(0);
        K[0][0] = K[1][1] = permX_[eIdx];
        K[2][2] = permZ_[eIdx];

        return K;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return porosity_[eIdx];
    }


private:
    std::vector<Scalar> porosity_;
    std::vector<Scalar> permX_;
    std::vector<Scalar> permZ_;
    static constexpr Scalar eps_ = 1.5e-7;
    bool homogeneous_;
};

} // end namespace Dumux

#endif
