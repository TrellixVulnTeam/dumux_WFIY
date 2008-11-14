/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief The required state classes for the lenhard problem
 */
#ifndef DUMUX_LENHARDSTATE_HH
#define DUMUX_LENHARDSTATE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/new_material/statehelpermacros.hh>
#include <dumux/new_material/parkerlenhardstate.hh>
#include <dumux/new_material/parkerlenhard.hh>

#include <vector>

namespace Dune
{
namespace Lenhard
{
    // parameters which are globally constant
    template <class ScalarT>
    class LenhardGlobalState
    {
    public:
        typedef ScalarT                          Scalar;
        typedef Dune::FieldVector<Scalar, 1>     NVector;
        typedef Dune::FieldMatrix<Scalar, 1, 1>  NxNMatrix;

        // Density of the wetting phase
        PROPERTY(Scalar, densityW, setDensityW);

        // Density of the non-wetting phase
        PROPERTY(Scalar, densityN, setDensityN);

        // the gravity vector
        PROPERTY(NVector, gravity, setGravity);

        // the viscosity of the wetting fluid
        PROPERTY(Scalar, viscosityW, setViscosityW);

        // the viscosity of the non-wetting fluid
        PROPERTY(Scalar, viscosityN, setViscosityN);
    };

    // parameters which are constant for a porous medium
    template <class CapPressureParamsT>
    class LenhardMediumState
    {
    public:
        typedef CapPressureParamsT                  CapPressureParams;
        typedef typename CapPressureParams::Scalar  Scalar;
        typedef Dune::FieldVector<Scalar, 1>        NVector;
        typedef Dune::FieldMatrix<Scalar, 1, 1>     NxNMatrix;
        typedef LenhardGlobalState<Scalar>          GlobalState;


        LenhardMediumState()
            {
                _Swr = _Snr = 0;
                _porosity = 0;
                _globalState = NULL;
            }

        // Low level hysteresis model parameters (i.e. parameters
        // of the underlying model for the main imbibition and
        // main drainage curves)
        MUTABLE_PROPERTY(CapPressureParams, micParams, setMicParams);
        MUTABLE_PROPERTY(CapPressureParams, mdcParams, setMdcParams);

        // Twophase parameters
        PARAMETER(Scalar, Swr);
        void setSwr(Scalar Swr)
            {
                _Swr = Swr;
                _Snre = _Snr/(1 - _Swr);
            }
        PARAMETER(Scalar, Snr);
        void setSnr(Scalar Snr){
            _Snr = Snr;
            _Snre = _Snr/(1 - _Swr);
        }
        PARAMETER(Scalar, Snre);

        // absolute permeability
        PROPERTY(NxNMatrix, permeability, setPermeability);
        void setPermeability(Scalar permeability)
            {
                _permeability[0][0] = permeability;
            }

        // the absolute porosity of the medium
        PROPERTY(Scalar, porosity, setPorosity);

        // pointer to the global properties
        PTR_PROPERTY(GlobalState, globalState, setGlobalState);
    };


#if !defined USE_NODE_PARAMETERS
    // node dependend parameters for cell centered parameter storage
    template <class MediumStateT>
    class LenhardNodeState
    {
    public:
        LenhardNodeState()
            {
                _neighbours = NULL;
            };

        ~LenhardNodeState()
            {
                delete _neighbours;
            };

        bool isOnInterface() const
            { return _neighbours != NULL; }

        void addNeighbourCellIdx(int cellIdx)
            {
                if (_neighbours == NULL) {
                    // allocate the array for the neighbouring cell
                    // indices.  in order not having to reallocate the
                    // array every time a new index is added, we
                    // reserve 4 slots from the beginning.
                    _neighbours = new _NeighbourIdxArray;
                    _neighbours->reserve(4);
                }

                // make sure the cellIdx is not already in the array
                for (unsigned i = 0; i < _neighbours->size(); ++i)
                    if ((*_neighbours)[i] == cellIdx)
                        return;

                // append the cell index
                _neighbours->push_back(cellIdx);
            }

        int numNeighbourCells() const
            {
                return _neighbours == NULL
                    ? 0
                    : _neighbours->size();
            }

        int neighbourCellIdx(int localIdx) const
            {
                assert(_neighbours != NULL);
                return _neighbours->at(localIdx);
            }

    private:
        typedef std::vector<int> _NeighbourIdxArray;
        _NeighbourIdxArray *_neighbours;
    };
#else // USE_NODE_PARAMETERS
    template <class MediumStateT>
    class LenhardCellState
    {
    public:
        typedef MediumStateT                              MediumState;

        typedef typename MediumState::Scalar              Scalar;
        typedef typename MediumState::NVector             NVector;
        typedef typename MediumState::NxNMatrix           NxNMatrix;

        // pointer to the medium uniform properties
        const MediumState *mediumState() const
            { return _mediumState; }
        void setMediumState(const MediumState *ms)
            { _mediumState = ms; }

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), NxNMatrix, permeability);
        PROXY_PARAMETER(mediumState(), Scalar, porosity);
        //////////////////////////////////

    private:
        const MediumState *_mediumState;
    };
#endif


    // cell/node dependent parameters
#if defined USE_NODE_PARAMETERS
#define MAIN_STATE_CLASS LenhardNodeState
#else
#define MAIN_STATE_CLASS LenhardCellState
#endif

    template <class MediumStateT>
    class MAIN_STATE_CLASS
    {
    public:
        typedef MediumStateT                              MediumState;
        typedef typename MediumState::GlobalState         GlobalState;

        typedef typename MediumState::Scalar              Scalar;
        typedef typename MediumState::NVector             NVector;
        typedef typename MediumState::NxNMatrix           NxNMatrix;
        typedef typename MediumState::CapPressureParams   CapPressureParams;
        typedef Dune::PLScanningCurve<Scalar>           ScanningCurve;

        MAIN_STATE_CLASS()
        {
            _init();
        };

        MAIN_STATE_CLASS(const MAIN_STATE_CLASS &s)
        {
            // we don't really copy the other cell state, but only
            // use the same medium and global parameters.  (The sole
            // purpose of the copy constructor is to make STL
            // containers work properly.)
            _init(s._mediumState);
        }

        MAIN_STATE_CLASS(const MediumState *mediumState)
        {
            _init(mediumState);
        }


        ~MAIN_STATE_CLASS()
            { }


        // pointer to the medium uniform properties
        const MediumState *mediumState() const
            { return _mediumState; }
        void setMediumState(const MediumState *ms)
            { _mediumState = ms; }

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), NxNMatrix, permeability);
        PROXY_PARAMETER(mediumState(), Scalar, porosity);
        //////////////////////////////////

        //////////////////////////////////
        // medium uniform parameters (implemented as proxy parameters)
        PROXY_PARAMETER(mediumState(), CapPressureParams, micParams);
        PROXY_PARAMETER(mediumState(), CapPressureParams, mdcParams);
        PROXY_PARAMETER(mediumState(), Scalar, Swr);
        PROXY_PARAMETER(mediumState(), Scalar, Snr);
        PROXY_PARAMETER(mediumState(), Scalar, Snre);
        //////////////////////////////////

        //////////////////////////////////
        // cell specific properties
        MUTABLE_PTR_PROPERTY(ScanningCurve, mdc, setMdc); // main drainage scanning curve
        MUTABLE_PTR_PROPERTY(ScanningCurve, pisc, setPisc); // primary imbibition scanning curve
        MUTABLE_PTR_PROPERTY(ScanningCurve, csc, setCsc); // current scanning curve
        PROPERTY(Scalar, Snrei, setSnrei); // current effective residual non-wetting saturation
        //////////////////////////////////

    private:
        void _init(const MediumState *mediumState = NULL)
            {
                _mediumState = mediumState;

                _Snrei = 0;
                _mdc = new ScanningCurve();
                _pisc = _csc = NULL;
            }

        const MediumState *_mediumState;
    };
}
}

#endif
