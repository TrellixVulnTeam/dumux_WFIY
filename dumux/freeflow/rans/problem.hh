// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSProblem
 */
#ifndef DUMUX_RANS_PROBLEM_HH
#define DUMUX_RANS_PROBLEM_HH

#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/deprecated.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/walldistance.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include "model.hh"

namespace Dumux {

//! forward declare
template<class TypeTag, TurbulenceModel turbulenceModel>
class RANSProblemImpl;

//! the turbulence-model-specfic RANS problem
template<class TypeTag>
using RANSProblem = RANSProblemImpl<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::turbulenceModel()>;

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes problem base class.
 *
 * This implements some base functionality for RANS models.
 * Especially vectors containing all wall-relevant properties, which are accessed
 * by the volumevariables.
 */
template<class TypeTag>
class RANSProblemBase : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr auto dim = GridView::dimension;
    static constexpr int numCorners = SubControlVolumeFace::numCornersPerFace;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    struct WallElementInformation
    {
        // store the element indicies for all elements with an intersection on the wall
        unsigned int wallElementIdx;
        // for each wall element, store the faces normal axis
        unsigned int wallFaceNormalAxis;
        // for each wall element, store the location of the face center and each corner.
        GlobalPosition wallFaceCenter;
        std::array<GlobalPosition, numCorners> wallFaceCorners;
    };

public:

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblemBase(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {
        if ( !(hasParamInGroup(this->paramGroup(), "RANS.IsFlatWallBounded")))
        {
            std::cout << "The parameter \"Rans.IsFlatWallBounded\" is not specified. \n"
                    << " -- Based on the grid and the boundary conditions specified by the user,"
                    << " this parameter is set to be "<< std::boolalpha << isFlatWallBounded() << "\n";
        }

        // update size and initial values of the global vectors
        wallDistance_.resize(this->gridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());
        neighborIdx_.resize(this->gridGeometry().elementMapper().size());
        velocity_.resize(this->gridGeometry().elementMapper().size(), DimVector(0.0));
        velocityGradients_.resize(this->gridGeometry().elementMapper().size(), DimMatrix(0.0));
        stressTensorScalarProduct_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        vorticityTensorScalarProduct_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        flowDirectionAxis_.resize(this->gridGeometry().elementMapper().size(), fixedFlowDirectionAxis_);
        storedViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedDensity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
    }

    /*!
     * \brief Update the static (solution independent) relations to the walls and neighbors
     */
    void updateStaticWallProperties()
    {
        std::cout << "Update static wall properties. ";
        calledUpdateStaticWallProperties = true;

        checkForWalls_();
        findWallDistances_();
        findNeighborIndices_();
    }

    /*!
     * \brief Update the dynamic (solution dependent) turbulence parameters
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        std::cout << "Update dynamic wall properties." << std::endl;
        if (!calledUpdateStaticWallProperties)
            DUNE_THROW(Dune::InvalidStateException,
                       "You have to call updateStaticWallProperties() once before you call updateDynamicWallProperties().");

        calculateCCVelocities_(curSol);
        calculateCCVelocityGradients_();
        calculateMaxMinVelocities_();
        calculateStressTensor_();
        calculateVorticityTensor_();
        storeViscosities_(curSol);
    }

    /*!
     * \brief Returns whether a wall function should be used at a given face
     *
     * \param element The element.
     * \param scvf The sub control volume face.
     * \param eqIdx The equation index.
     */
    bool useWallFunction(const Element& element,
                         const SubControlVolumeFace& scvf,
                         const int& eqIdx) const
    { return false; }

    /*!
     * \brief Returns an additional wall function momentum flux
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    FacePrimaryVariables wallFunction(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFaceVariables& elemFaceVars,
                                      const SubControlVolumeFace& scvf,
                                      const SubControlVolumeFace& lateralBoundaryFace) const
    { return FacePrimaryVariables(0.0); }

    /*!
     * \brief  Returns an additional wall function flux for cell-centered quantities
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    { return CellCenterPrimaryVariables(0.0); }

    /*!
     * \brief Returns whether a given sub control volume face is on a wall
     *
     * \param scvf The sub control volume face.
     */
    [[deprecated("The isOnWall and IsOnWallAtPos functions will be removed after release 3.5. "
                 "Please use the Rans specific boundarytypes, and mark wall boundaries with the setWall() function.")]]
    bool isOnWall(const SubControlVolumeFace& scvf) const
    {
        return asImp_().isOnWallAtPos(scvf.center());
    }

    bool isFlatWallBounded() const
    {
        static const bool hasAlignedWalls = hasAlignedWalls_();
        return hasAlignedWalls;
    }

    /*!
     * \brief Returns the Karman constant
     */
    const Scalar karmanConstant() const
    { return 0.41; }

    //! \brief Returns the \f$ \beta_{\omega} \f$ constant
    const Scalar betaOmega() const
    {
        return 0.0708;
    }

    /*!
     * \brief Return the turbulent Prandtl number \f$ [-] \f$ which is used to convert
     *        the eddy viscosity to an eddy thermal conductivity
     */
    Scalar turbulentPrandtlNumber() const
    {
        static const Scalar turbulentPrandtlNumber
            = getParamFromGroup<Scalar>(this->paramGroup(), "RANS.TurbulentPrandtlNumber", 1.0);
        return turbulentPrandtlNumber;
    }

    /*!
     * \brief Return the turbulent Schmidt number \f$ [-] \f$ which is used to convert
     *        the eddy viscosity to an eddy diffusivity
     */
    Scalar turbulentSchmidtNumber() const
    {
        static const Scalar turbulentSchmidtNumber
            = getParamFromGroup<Scalar>(this->paramGroup(), "RANS.TurbulentSchmidtNumber", 1.0);
        return turbulentSchmidtNumber;
    }

    int wallNormalAxis(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a wallNormalAxis "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return wallNormalAxis_[elementIdx];
    }

    int flowDirectionAxis(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a flowDirectionAxis "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return flowDirectionAxis_[elementIdx];
    }

    unsigned int wallElementIndex(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a wallElementIndex "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return wallElementIdx_[elementIdx];

    }

    Scalar wallDistance(const int elementIdx) const
    { return wallDistance_[elementIdx]; }

    GlobalPosition cellCenter(const int elementIdx) const
    {
        const auto& element = this->gridGeometry().element(elementIdx);
        return element.geometry().center();
    }

    unsigned int neighborIndex(const int elementIdx, const int dimIdx, const int sideIdx) const
    { return neighborIdx_[elementIdx][dimIdx][sideIdx];}

    DimVector ccVelocityVector(const int elementIdx) const
    { return velocity_[elementIdx]; }

    Scalar ccVelocity(const int elementIdx, const int dimIdx) const
    { return velocity_[elementIdx][dimIdx]; }

    DimVector velocityMaximum(const int elementIdx) const
    { return velocityMaximum_[elementIdx]; }

    DimVector velocityMinimum(const int elementIdx) const
    { return velocityMinimum_[elementIdx]; }

    DimMatrix velocityGradientTensor(const int elementIdx) const
    { return velocityGradients_[elementIdx]; }

    Scalar velocityGradient(const int elementIdx, const int dimIdx, const int velIdx) const
    { return velocityGradients_[elementIdx][dimIdx][velIdx]; }

    Scalar stressTensorScalarProduct(const int elementIdx) const
    { return stressTensorScalarProduct_[elementIdx]; }

    Scalar vorticityTensorScalarProduct(const int elementIdx) const
    { return vorticityTensorScalarProduct_[elementIdx]; }

    Scalar storedViscosity(const int elementIdx) const
    { return storedViscosity_[elementIdx]; }

    Scalar storedDensity(const int elementIdx) const
    { return storedDensity_[elementIdx]; }

    Scalar kinematicViscosity(const int elementIdx) const
    { return storedViscosity(elementIdx) / storedDensity(elementIdx); }

    bool calledUpdateStaticWallProperties = false;

private:

    bool hasAlignedWalls_() const
    {
        if ( hasParamInGroup(this->paramGroup(), "RANS.IsFlatWallBounded"))
        {
            static const bool isFlatWallBounded = getParamFromGroup<bool>(this->paramGroup(), "RANS.IsFlatWallBounded");
            return isFlatWallBounded;
        }

        std::vector<int> wallFaceAxis;
        wallFaceAxis.reserve(this->gridGeometry().numBoundaryScvf());

        const auto gridView = this->gridGeometry().gridView();
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // Remove this check after release 3.5. IsOnWall Interface is deprecated
                if constexpr (Deprecated::hasIsOnWall<Implementation, GlobalPosition>())
                {
                    // Remove this part
                    if (!scvf.boundary() && asImp_().isOnWall(scvf)) // only search for walls at a global boundary
                        wallFaceAxis.push_back(scvf.directionIndex());
                }
                else
                {
                    // Keep this part
                    if (!scvf.boundary() && asImp_().boundaryTypes(element, scvf).hasWall())  // only search for walls at a global boundary
                        wallFaceAxis.push_back(scvf.directionIndex());
                }
            }
        }

        // Returns if all wall directions are the same
        return std::all_of(wallFaceAxis.begin(), wallFaceAxis.end(), [firstDir=wallFaceAxis[0]](auto dir){ return (dir == firstDir);} ) ;
    }

    void checkForWalls_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
                if constexpr (!Deprecated::hasIsOnWall<Implementation, GlobalPosition>())
                    if (asImp_().boundaryTypes(element, scvf).hasWall())
                        return;
        }

        // If reached, no walls were found, throw exception. Remove check after 3.5
        if constexpr (!Deprecated::hasIsOnWall<Implementation, GlobalPosition>())
            DUNE_THROW(Dune::InvalidStateException, "No walls are are specified with the setWall() function");
    }

    /*!
     * \brief Use the boundary search algorithm to find the shortest distance to a wall for each element
     *
     *  Also store the wall element's index, and its direction in the case of flat wall bounded problems
     */
    void findWallDistances_()
    {
        // Remove this check after release 3.5. IsOnWall Interface is deprecated
        if constexpr (Deprecated::hasIsOnWall<Implementation, GlobalPosition>())
        {
            WallDistance wallInformation(this->gridGeometry(), WallDistance<GridGeometry>::atElementCenters,
                [this] (const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
                { return asImp_().isOnWall(scvf); });
            wallDistance_ = wallInformation.wallDistance();
            storeWallElementAndDirectionIndex_(wallInformation.wallData());
        }
        else
        {
            WallDistance wallInformation(this->gridGeometry(), WallDistance<GridGeometry>::atElementCenters,
                [this] (const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
                { return asImp_().boundaryTypes(fvGeometry.element(), scvf).hasWall(); });
            wallDistance_ = wallInformation.wallDistance();
            storeWallElementAndDirectionIndex_(wallInformation.wallData());
        }
    }

    template <class WallData>
    void storeWallElementAndDirectionIndex_(const WallData& wallData)
    {
        // The wall Direction Index is used for flat quadrilateral channel problems only
        if (!(GridGeometry::discMethod == DiscretizationMethods::staggered))
            DUNE_THROW(Dune::NotImplemented, "The wall direction Index can only be calculated for quadrilateral structured grids");

        // If isFlatWallBounded, the corresonding wall element is stored for each element
        if (isFlatWallBounded())
        {
            wallNormalAxis_.resize(wallData.size());
            wallElementIdx_.resize(wallData.size());

            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                wallElementIdx_[elementIdx] = wallData[elementIdx].eIdx;
                if ( ! (hasParam("RANS.WallNormalAxis")) )
                {
                    GlobalPosition wallOuterNormal = wallData[elementIdx].scvfOuterNormal;
                    if constexpr (dim == 2) // 2D
                        wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : 1;
                    else // 3D
                        wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : ((wallOuterNormal[1] == 1) ? 1 : 2);
                }
                else
                    wallNormalAxis_[elementIdx] = fixedWallNormalAxis_;
            }
        }
    }

    /*!
     * \brief Store all direct neighbor indicies for each element
     */
    void findNeighborIndices_()
    {
        // search for neighbor Idxs
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborIdx_[elementIdx][dimIdx][0] = elementIdx;
                neighborIdx_[elementIdx][dimIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(this->gridGeometry().gridView(), element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborIdx = this->gridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter(elementIdx)[dimIdx] - cellCenter(neighborIdx)[dimIdx]) > 1e-8)
                    {
                        if (cellCenter(elementIdx)[dimIdx] > cellCenter(neighborIdx)[dimIdx])
                            neighborIdx_[elementIdx][dimIdx][0] = neighborIdx;

                        if (cellCenter(elementIdx)[dimIdx] < cellCenter(neighborIdx)[dimIdx])
                            neighborIdx_[elementIdx][dimIdx][1] = neighborIdx;
                    }
                }
            }
        }
    }

    void calculateCCVelocities_(const SolutionVector& curSol)
    {
        auto fvGeometry = localView(this->gridGeometry());
        // calculate cell-center-averaged velocities
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            // calculate velocities
            DimVector velocityTemp(0.0);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const int dofIdxFace = scvf.dofIndex();
                const auto numericalSolutionFace = curSol[GridGeometry::faceIdx()][dofIdxFace][Indices::velocity(scvf.directionIndex())];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementIdx][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }
    }


    void calculateCCVelocityGradients_()
    {
        using std::abs;

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    const unsigned int neighborIndex0 = neighborIndex(elementIdx, dimIdx, 0);
                    const unsigned int neighborIndex1 = neighborIndex(elementIdx, dimIdx, 1);

                    velocityGradients_[elementIdx][velIdx][dimIdx]
                        = (ccVelocity(neighborIndex1, velIdx) - ccVelocity(neighborIndex0, velIdx))
                        / (cellCenter(neighborIndex1)[dimIdx] - cellCenter(neighborIndex0)[dimIdx]);

                    if (abs(cellCenter(neighborIndex1)[dimIdx] - cellCenter(neighborIndex0)[dimIdx]) < 1e-8)
                        velocityGradients_[elementIdx][velIdx][dimIdx] = 0.0;
                }
            }

            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // adapt calculations for Dirichlet condition
                unsigned int scvfNormDim = scvf.directionIndex();
                if (scvf.boundary())
                {
                    for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                    {
                        if (!asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::velocity(velIdx)))
                            continue;

                        Scalar dirichletVelocity = asImp_().dirichlet(element, scvf)[Indices::velocity(velIdx)];

                        unsigned int neighborIdx = neighborIndex(elementIdx, scvfNormDim, 0);
                        if (scvf.center()[scvfNormDim] < cellCenter(elementIdx)[scvfNormDim])
                            neighborIdx = neighborIndex(elementIdx, scvfNormDim, 1);

                        velocityGradients_[elementIdx][velIdx][scvfNormDim]
                            = (ccVelocity(neighborIdx, velIdx) - dirichletVelocity)
                              / (cellCenter(neighborIdx)[scvfNormDim] - scvf.center()[scvfNormDim]);
                    }
                }

                // Calculate the BJS-velocity by accounting for all sub faces.
                std::vector<int> bjsNumFaces(dim, 0);
                std::vector<unsigned int> bjsNeighbor(dim, 0);
                DimVector bjsVelocityAverage(0.0);
                DimVector normalNormCoordinate(0.0);
                unsigned int velIdx = Indices::velocity(scvfNormDim);
                const int numSubFaces = scvf.pairData().size();
                for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                {
                    const auto& lateralFace = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData()[localSubFaceIdx].localLateralFaceIdx);

                    // adapt calculations for Beavers-Joseph-Saffman condition
                    unsigned int normalNormDim = lateralFace.directionIndex();
                    if (lateralFace.boundary() && (asImp_().boundaryTypes(element, lateralFace).isBeaversJoseph(Indices::velocity(velIdx))))
                    {
                        unsigned int neighborIdx = neighborIndex(elementIdx, normalNormDim, 0);
                        if (lateralFace.center()[normalNormDim] < cellCenter(elementIdx)[normalNormDim])
                            neighborIdx = neighborIndex(elementIdx, normalNormDim, 1);

                        const SubControlVolume& scv = fvGeometry.scv(scvf.insideScvIdx());
                        bjsVelocityAverage[normalNormDim] += ParentType::beaversJosephVelocity(element, scv, scvf, lateralFace, ccVelocity(elementIdx, velIdx), 0.0);
                        if (bjsNumFaces[normalNormDim] > 0 && neighborIdx != bjsNeighbor[normalNormDim])
                            DUNE_THROW(Dune::InvalidStateException, "Two different neighborIdx should not occur");
                        bjsNeighbor[normalNormDim] = neighborIdx;
                        normalNormCoordinate[normalNormDim] = lateralFace.center()[normalNormDim];
                        bjsNumFaces[normalNormDim]++;
                    }
                }
                for (unsigned dirIdx = 0; dirIdx < dim; ++dirIdx)
                {
                    if (bjsNumFaces[dirIdx] == 0)
                        continue;

                    unsigned int neighborIdx = bjsNeighbor[dirIdx];
                    bjsVelocityAverage[dirIdx] /= bjsNumFaces[dirIdx];

                    velocityGradients_[elementIdx][velIdx][dirIdx]
                        = (ccVelocity(neighborIdx, velIdx) - bjsVelocityAverage[dirIdx])
                        / (cellCenter(neighborIdx)[dirIdx] - normalNormCoordinate[dirIdx]);

                }
            }
        }
    }

    void calculateMaxMinVelocities_()
    {
        using std::abs;
        if (isFlatWallBounded())
        {
            // If the parameter isFlatWallBounded is set to true,
            // the maximum/minimum velocities are calculated along a profile perpendicular to the corresponding wall face.

            // re-initialize min and max values
            velocityMaximum_.assign(this->gridGeometry().elementMapper().size(), DimVector(1e-16));
            velocityMinimum_.assign(this->gridGeometry().elementMapper().size(), DimVector(std::numeric_limits<Scalar>::max()));

            // For each profile perpendicular to the channel wall, find the max and minimum velocities
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                Scalar maxVelocity = 0.0;
                const unsigned int wallElementIdx = wallElementIndex(elementIdx);

                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(ccVelocity(elementIdx, dimIdx)) > abs(velocityMaximum_[wallElementIdx][dimIdx]))
                        velocityMaximum_[wallElementIdx][dimIdx] = ccVelocity(elementIdx, dimIdx);

                    if (abs(ccVelocity(elementIdx, dimIdx)) < abs(velocityMinimum_[wallElementIdx][dimIdx]))
                        velocityMinimum_[wallElementIdx][dimIdx] = ccVelocity(elementIdx, dimIdx);

                    // Set the flow direction axis as the direction of the max velocity
                    if ((hasParam("RANS.FlowDirectionAxis") != 1) && (maxVelocity) < abs(ccVelocity(elementIdx, dimIdx)))
                    {
                        maxVelocity = abs(ccVelocity(elementIdx, dimIdx));
                        flowDirectionAxis_[elementIdx] = dimIdx;
                    }
                }
            }
        }
        else
        {
            // If the parameter isFlatWallBounded is set to false, or not set,
            // the maximum/minimum velocities are calculated as a global max/min throughout the domain.

            DimVector maxVelocity(0.0);
            DimVector minVelocity(std::numeric_limits<Scalar>::max());
            // Find the max and minimum velocities in the full domain
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(ccVelocity(elementIdx, dimIdx)) > abs(maxVelocity[dimIdx]))
                        maxVelocity[dimIdx] = ccVelocity(elementIdx, dimIdx);

                    if (abs(ccVelocity(elementIdx, dimIdx)) < abs(minVelocity[dimIdx]))
                        minVelocity[dimIdx] = ccVelocity(elementIdx, dimIdx);
                }
            }
            velocityMaximum_.assign(this->gridGeometry().elementMapper().size(), maxVelocity);
            velocityMinimum_.assign(this->gridGeometry().elementMapper().size(), minVelocity);
        }
    }

    void calculateStressTensor_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> stressTensor(0.0);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensor[dimIdx][velIdx] = 0.5 * velocityGradient(elementIdx, dimIdx, velIdx)
                                                 + 0.5 * velocityGradient(elementIdx, velIdx, dimIdx);
              }
            }
            stressTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensorScalarProduct_[elementIdx] += stressTensor[dimIdx][velIdx] * stressTensor[dimIdx][velIdx];
                }
            }
        }
    }

    void calculateVorticityTensor_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> vorticityTensor(0.0);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    vorticityTensor[dimIdx][velIdx] = 0.5 * velocityGradient(elementIdx, dimIdx, velIdx)
                                                    - 0.5 * velocityGradient(elementIdx, velIdx, dimIdx);
              }
            }
            vorticityTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    vorticityTensorScalarProduct_[elementIdx] += vorticityTensor[dimIdx][velIdx] * vorticityTensor[dimIdx][velIdx];
                }
            }
        }
    }

    void storeViscosities_(const SolutionVector& curSol)
    {
        // calculate or call all secondary variables
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));

                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDensity_[elementIdx] = volVars.density();
                storedViscosity_[elementIdx] = volVars.viscosity();
            }
        }
    }

    const int fixedFlowDirectionAxis_ = getParam<int>("RANS.FlowDirectionAxis", 0);
    const int fixedWallNormalAxis_ = getParam<int>("RANS.WallNormalAxis", 1);

    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowDirectionAxis_;
    std::vector<Scalar> wallDistance_;
    std::vector<unsigned int> wallElementIdx_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIdx_;

    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;

    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;

    std::vector<Scalar> storedDensity_;
    std::vector<Scalar> storedViscosity_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
