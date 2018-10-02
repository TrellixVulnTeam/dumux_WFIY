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
 * \ingroup OnePTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *        The fracture is a 2D network embedded in 3D.
 */

#ifndef DUMUX_ONEP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_ONEP_FRACTURE_TEST_PROBLEM_HH

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>

#include "fracturespatialparams.hh"

namespace Dumux {

template <class TypeTag>
class FractureProblem;

namespace Properties {
NEW_TYPE_TAG(FractureTypeTag, INHERITS_FROM(OneP));
NEW_TYPE_TAG(FractureBoxTypeTag, INHERITS_FROM(BoxModel, FractureTypeTag));
NEW_TYPE_TAG(FractureCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, FractureTypeTag));
NEW_TYPE_TAG(FractureCCMpfaTypeTag, INHERITS_FROM(CCMpfaModel, FractureTypeTag));

//! Enable caching (more memory, but faster runtime)
SET_BOOL_PROP(FractureTypeTag, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(FractureTypeTag, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(FractureTypeTag, EnableGridFluxVariablesCache, true);

//! The grid type
#if HAVE_DUNE_FOAMGRID
SET_TYPE_PROP(FractureTypeTag, Grid, Dune::FoamGrid<2, 3>);
#endif

// Set the problem property
SET_TYPE_PROP(FractureTypeTag, Problem, Dumux::FractureProblem<TypeTag>);

// the fluid system
SET_PROP(FractureTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};
// Set the spatial parameters
SET_PROP(FractureTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FractureSpatialParams<FVGridGeometry, Scalar>;
};

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *        The fracture is a 2D network embedded in 3D.
 *
 * This problem uses the \ref OnePModel.
 */
template <class TypeTag>
class FractureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    enum { dimWorld = GridView::dimensionworld };

    enum {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

public:
    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    FractureProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    // void postTimeStep()
    // {
    //     // Calculate storage terms
    //     PrimaryVariables storage;
    //     this->model().globalStorage(storage);
    //
    //     // Write mass balance information for rank 0
    //     if (this->gridView().comm().rank() == 0) {
    //         std::cout<<"Storage: " << storage << std::endl;
    //     }
    // }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();
        const auto& gg = this->fvGridGeometry();
        if (globalPos[0] > gg.bBoxMax()[0] - eps_ || globalPos[0] < gg.bBoxMin()[0] + eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        const auto& gg = this->fvGridGeometry();
        values[pressureIdx] = 1.0e5*(globalPos[0] - gg.bBoxMin()[0])/(gg.bBoxMax()[0] - gg.bBoxMin()[0]) + 1.0e5;
        return values;
    }
    // \}

private:
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} //end namespace Dumux

#endif
