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
 * \brief Definition of the spatial parameters for the kuevette problem
 *        which uses the isothermal two-phase two component
 *        fully implicit model.
 */
#ifndef DUMUX_INFILTRATION_SPATIAL_PARAMETERS_HH
#define DUMUX_INFILTRATION_SPATIAL_PARAMETERS_HH

#include <dumux/implicit/3p3c/3p3cindices.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangen3pparams.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class InfiltrationSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(InfiltrationSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(InfiltrationSpatialParams, SpatialParams, Dumux::InfiltrationSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(InfiltrationSpatialParams, MaterialLaw)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    // define the material law
    typedef ParkerVanGen3P<Scalar> type;
};
}

/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class TypeTag>
class InfiltrationSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
	dim = GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar,dim> DimVector;


    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;


public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InfiltrationSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        // intrinsic permeabilities
        fineK_ = 1.e-11;
        coarseK_ = 1.e-11;

        // porosities
        Porosity_ = 0.40;

        // residual saturations
        MaterialParams_.setSwr(0.12);
        MaterialParams_.setSwrx(0.12);
        MaterialParams_.setSnr(0.07);
        MaterialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        MaterialParams_.setVgAlpha(0.0005);
        MaterialParams_.setVgn(4.);
        MaterialParams_.setKrRegardsSnr(false);

        // parameters for adsorption
        MaterialParams_.setKdNAPL(0.);
        MaterialParams_.setRhoBulk(1500.);
    }

    ~InfiltrationSpatialParams()
    {}


    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSolution The global solution vector
     */
    void update(const SolutionVector &globalSolution)
    {
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        //const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        // if (isFineMaterial_(globalPos))
        //     return finePorosity_;
        // else
        //     return coarsePorosity_;
        return Porosity_;
    }


    /*!
     * \brief return the parameter object for the material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               int scvIdx) const
    {
        //const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
        //if (isFineMaterial_(globalPos))
        //return fineMaterialParams_;
        //else
        //return coarseMaterialParams_;
        return MaterialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    double heatCapacity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        int scvIdx) const
    {
        return
            850. // specific heat capacity [J / (kg K)]
            * 2650. // density of sand [kg/m^3]
            * (1 - porosity(element, fvGeometry, scvIdx));
    }

    /*!
     * \brief Calculate the heat flux \f$[W/m^2]\f$ through the
     *        rock matrix based on the temperature gradient \f$[K / m]\f$
     *
     * This is only required for non-isothermal models.
     *
     * \param heatFlux The resulting heat flux vector
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables
     * \param tempGrad The temperature gradient
     * \param element The current finite element
     * \param fvGeometry The finite volume geometry of the current element
     * \param scvfIdx The local index of the sub-control volume face where
     *                    the matrix heat flux should be calculated
     */
    void matrixHeatFlux(DimVector &heatFlux,
                        const FluxVariables &fluxVars,
                        const ElementVolumeVariables &elemVolVars,
                        const DimVector &tempGrad,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        int scvfIdx) const
    {
        static const Scalar lDry = 0.35;
        static const Scalar lsw1 = 1.8;
        static const Scalar lsn1 = 0.65;

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fluxVars.face().i;
        const int j = fluxVars.face().j;
        Scalar sw = std::max(0.0, (elemVolVars[i].saturation(wPhaseIdx) +
                                   elemVolVars[j].saturation(wPhaseIdx)) / 2);
        Scalar sn = std::max(0.0, (elemVolVars[i].saturation(nPhaseIdx) +
                                   elemVolVars[j].saturation(nPhaseIdx)) / 2);

        // the heat conductivity of the matrix. in general this is a
        // tensorial value, but we assume isotropic heat conductivity.
        Scalar heatCond = lDry + sqrt(sw) * (lsw1-lDry) + sqrt(sn) * (lsn1-lDry);

        // the matrix heat flux is the negative temperature gradient
        // times the heat conductivity.
        heatFlux = tempGrad;
        heatFlux *= -heatCond;
    }

    const MaterialLawParams& materialLawParams() const
    {
        return MaterialParams_;
    }
private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return
            70. <= globalPos[0] && globalPos[0] <= 85. &&
            7.0 <= globalPos[1] && globalPos[1] <= 7.50;
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar Porosity_;

    MaterialLawParams MaterialParams_;
};

}

#endif
