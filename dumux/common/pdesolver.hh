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
 * \brief Defines a high-level interface for a PDESolver
 */
#ifndef DUMUX_COMMON_PDESOLVER_HH
#define DUMUX_COMMON_PDESOLVER_HH

#include <memory>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/type_traits.hh>

#include <dumux/common/timeloop.hh>

// forward declare
namespace Dune {
template <class FirstRow, class ... Args>
class MultiTypeBlockMatrix;
} // end namespace Dune

namespace Dumux {
namespace Detail {

template<class Assembler>
using AssemblerVariablesType = typename Assembler::Variables;

template<class Assembler>
inline constexpr bool exportsVariables = Dune::Std::is_detected_v<AssemblerVariablesType, Assembler>;

template<class A, bool exports = exportsVariables<A>> struct VariablesChooser;
template<class A> struct VariablesChooser<A, true> { using Type = AssemblerVariablesType<A>; };
template<class A> struct VariablesChooser<A, false> { using Type = typename A::ResidualType; };

template<class Assembler>
using AssemblerVariables = typename VariablesChooser<Assembler>::Type;

} // end namespace Detail

/*!
 * \ingroup Common
 * \brief A high-level interface for a PDESolver
 *
 * A PDESolver is constructed with an assembler and a linear solver
 * and has a method solve that linearizes (if not already linear), assembles, solves and updates
 * given an initial solution producing a new solution.
 *
 * \tparam A Assembler for linearized system of the PDE
 * \tparam LS Linear system solver
 */
template<class A, class LS>
class PDESolver
{
    using Scalar = typename A::Scalar;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! export the assembler and linear solver types
    using Assembler = A;
    using LinearSolver = LS;

    //! export the type of variables that represent a numerical solution
    using Variables = Detail::AssemblerVariables<Assembler>;

    /*!
     * \brief Constructor
     * \param assembler pointer to the assembler of the linear system
     * \param linearSolver pointer to the solver of the resulting linear system
     */
    PDESolver(std::shared_ptr<Assembler> assembler,
              std::shared_ptr<Assembler> assemblerRef,
              std::shared_ptr<LinearSolver> linearSolver)
    : assembler_(assembler)
    , assemblerRef_(assemblerRef)
    , linearSolver_(linearSolver)
    {}

    virtual ~PDESolver() = default;

    /*!
     * \brief Solve the given PDE system (usually assemble + solve linear system + update)
     * \param vars instance of the `Variables` class representing a numerical
     *             solution, defining primary and possibly secondary variables
     *             and information on the time level.
     */
    virtual void solve(Variables& vars) = 0;

    /*!
     * \brief Solve the given PDE system with time step control
     * \note This is used for solvers that are allowed to e.g. automatically reduce the
     *       time step if the solve was not successful
     * \param vars instance of the `Variables` class representing a numerical solution
     * \param timeLoop a reference to the current time loop
     */
    virtual void solve(Variables& vars, TimeLoop& timeLoop)
    {
        // per default we just forward to the method without time step control
        solve(vars);
    }

    /*!
     * \brief Access the assembler
     */
    const Assembler& assembler() const
    { return *assembler_; }

    /*!
     * \brief Access the assembler
     */
    Assembler& assembler()
    { return *assembler_; }

    /*!
     * \brief Access the assemblerRef
     */
    const Assembler& assemblerRef() const
    { return *assemblerRef_; }

    /*!
     * \brief Access the assembler
     */
    Assembler& assemblerRef()
    { return *assemblerRef_; }

    /*!
     * \brief Access the linear solver
     */
    const LinearSolver& linearSolver() const
    { return *linearSolver_; }

protected:

    /*!
     * \brief Access the linear solver
     */
    LinearSolver& linearSolver()
    { return *linearSolver_; }

    /*!
     * \brief Helper function to assure the MultiTypeBlockMatrix's sub-blocks have the correct sizes.
     */
    template <class FirstRow, class ... Args>
    bool checkSizesOfSubMatrices(const Dune::MultiTypeBlockMatrix<FirstRow, Args...>& matrix) const
    {
        bool matrixHasCorrectSize = true;
        using namespace Dune::Hybrid;
        forEach(std::make_index_sequence<Dune::MultiTypeBlockMatrix<FirstRow, Args...>::N()>(), [&](const auto i)
        {
            const auto& row = matrix[i];
            const auto numRowsLeftMostBlock = row[Dune::index_constant<0>{}].N();
            forEach(row, [&](const auto& subBlock)
            {
                if (subBlock.N() != numRowsLeftMostBlock)
                    matrixHasCorrectSize = false;
            });
        });
        return matrixHasCorrectSize;
    }

private:
    std::shared_ptr<Assembler> assembler_;
    std::shared_ptr<Assembler> assemblerRef_;
    std::shared_ptr<LinearSolver> linearSolver_;
};

} // namespace Dumux

#endif
