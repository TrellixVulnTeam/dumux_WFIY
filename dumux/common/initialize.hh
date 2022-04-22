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
 * \ingroup Common
 * \brief The initialize function to be called before using Dumux
 */
#ifndef DUMUX_COMMON_INITIALIZE_HH
#define DUMUX_COMMON_INITIALIZE_HH

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_TBB
#include <cstdlib>
#include <string>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/global_control.h>

#ifndef DOXYGEN
namespace Dumux::Detail {

class TBBGlobalControl
{
public:
    static oneapi::tbb::global_control& instance(int& argc, char* argv[])
    {
        int maxNumThreads = oneapi::tbb::info::default_concurrency();
        if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
            maxNumThreads = std::max(1, std::stoi(std::string{ dumuxNumThreads }));

        static oneapi::tbb::global_control global_limit(
            oneapi::tbb::global_control::max_allowed_parallelism, maxNumThreads
        );

        return global_limit;
    }
};

} // namespace Dumux::Detail
#endif // DOXYGEN

#endif // HAVE_TBB


#if HAVE_OPENMP
#include <omp.h>
#endif // HAVE_OPENMP


#if HAVE_KOKKOS
#include <Kokkos_Core.hpp>

#ifndef DOXYGEN
namespace Dumux::Detail {

class KokkosScopeGuard
{
public:
    static Kokkos::ScopeGuard& instance(int& argc, char* argv[])
    {
        static Kokkos::ScopeGuard guard(argc, argv);
        return guard;
    }
};

} // namespace Dumux::Detail
#endif // DOXYGEN

#endif // HAVE_KOKKOS

namespace Dumux {

void initialize(int& argc, char* argv[])
{
    // initialize MPI if available
    // otherwise this will create a sequential (fake) helper
    Dune::MPIHelper::instance(argc, argv);

#if HAVE_TBB
    // initialize TBB and keep global control alive
    Detail::TBBGlobalControl::instance(argc, argv);
#endif

#if HAVE_OPENMP
    if (const char* dumuxNumThreads = std::getenv("DUMUX_NUM_THREADS"))
        omp_set_num_threads(
            std::max(1, std::stoi(std::string{ dumuxNumThreads }))
        );
#endif

#if HAVE_KOKKOS
    // initialize Kokkos (command line / environmental variable interface)
    Detail::KokkosScopeGuard::instance(argc, argv);
#endif
}

} // end namespace Dumux

#endif
