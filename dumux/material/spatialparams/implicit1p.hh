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
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
#ifndef DUMUX_IMPLICIT_SPATIAL_PARAMS_ONE_P_HH
#define DUMUX_IMPLICIT_SPATIAL_PARAMS_ONE_P_HH

#include <dumux/common/propertysystem.hh>
#include <dumux/common/math.hh>

#include <dumux/implicit/properties.hh>

#include <dune/common/fmatrix.hh>

namespace Dumux {
// forward declaration of property tags
namespace Properties {
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(SpatialParamsForchCoeff);
}

/*!
 * \ingroup SpatialParameters
 */


/**
 * \brief The base class for spatial parameters of one-phase problems
 * using a fully implicit discretization method.
 */
template<class TypeTag>
class ImplicitSpatialParamsOneP
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Implementation = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<CoordScalar,dimWorld>;

public:
    ImplicitSpatialParamsOneP(const Problem& problem, const GridView &gridView)
    : problemPtr_(&problem)
    {}

    /*!
     * \brief Averages a diffusion coefficient (Scalar).
     * \return the averaged diffusion coefficient
     * \param T1 diffusion coefficient of the first scv
     * \param T2 diffusion coefficient of the second scv
     * \param normal The unit normal vector of the face on which to average
     */
    Scalar meanDiffusionTensor(const Scalar T1,
                               const Scalar T2,
                               const GlobalPosition& normal) const
    { return harmonicMean(T1, T2); }

    /*!
     * \brief Averages a diffusion Tensor.
     * \return the averaged diffusion tensor
     * \param T1 diffusion tensor of the first scv
     * \param T2 diffusion tensor of the second scv
     * \param normal The unit normal vector of the face on which to average
     */
    DimWorldMatrix meanDiffusionTensor(const DimWorldMatrix& T1,
                                       const DimWorldMatrix& T2,
                                       const GlobalPosition& normal) const
    {
        // determine nT*k*n
        GlobalPosition tmp;
        GlobalPosition tmp2;
        T1.mv(normal, tmp);
        T2.mv(normal, tmp2);
        Scalar alpha1 = tmp*normal;
        Scalar alpha2 = tmp2*normal;

        Scalar alphaHarmonic = harmonicMean(alpha1, alpha2);
        Scalar alphaAverage = 0.5*(alpha1 + alpha2);

        DimWorldMatrix T(0.0);
        for (int i = 0; i < dimWorld; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
            {
                T[i][j] += 0.5*(T1[i][j] + T2[i][j]);
                if (i == j)
                    T[i][j] += alphaHarmonic - alphaAverage;
            }
        }

        return T;
    }

    /*!
     * \brief Function for defining the permeability.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    auto permeability(const Element& element,
                      const SubControlVolume& scv,
                      const ElementSolutionVector& elemSol) const
    {
        return asImp_().permeabilityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the permeability.
     *
     * \return permeability
     * \param globalPos The position of the center of the scv
     */
    DimWorldMatrix permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a intrinsicPermeabilityAtPos() method.");
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolutionVector& elemSol) const
    {
        return asImp_().porosityAtPos(scv.center());
    }

    /*!
     * \brief Function for defining the porosity.
     *
     * \return porosity
     * \param globalPos The position of the center of the scv
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a porosityAtPos() method.");
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidHeatCapacity(const Element &element,
                             const SubControlVolume& scv,
                             const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidHeatCapacityAtPos(scv.center());
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidHeatCapacityAtPos() method.");
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidDensity(const Element &element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidDensityAtPos(scv.center());
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidDensityAtPos() method.");
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv,
                                    const ElementSolutionVector& elemSol) const
    {
        return asImp_().solidThermalConductivityAtPos(scv.center());
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a solidThermalConductivityAtPos() method.");
    }

    /*!
     * \brief Function for defining the Beavers-Joseph coefficient for multidomain
     *        problems\f$\mathrm{[-]}\f$.
     *
     * \return Beavers-Joseph coefficient \f$\mathrm{[-]}\f$
     * \param globalPos The global position
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide a beaversJosephCoeffAtPos() method.");
    }

    /*!
     * \brief Apply the Forchheimer coefficient for inertial forces
     *        calculation.
     *
     *        Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90 \cite ward1964 .
     *        Actually the Forchheimer coefficient is also a function of the dimensions of the
     *        porous medium. Taking it as a constant is only a first approximation
     *        (Nield, Bejan, Convection in porous media, 2006, p. 10 \cite nield2006 )
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    Scalar forchCoeff(const SubControlVolume &scv) const
    {
        static Scalar forchCoeff = GET_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, ForchCoeff);
        return forchCoeff;
    }

    const Problem& problem()
    {
        return *problemPtr_;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Problem *problemPtr_;
};

} // namespace Dumux

#endif
