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
 * \copydoc Dumux::Variables
 */
#ifndef DUMUX_VARIABLES_HH
#define DUMUX_VARIABLES_HH

#include <type_traits>

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/timestepping/timelevel.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Discretization
 * \brief Class that represents the variables of a model.
 *        We assume that models are formulated on the basis of primary and
 *        possibly secondary variables, where the latter may non-linearly
 *        depend on the former. Variables in Dumux represent the state of
 *        a numerical solution of a model, consisting of all primary/secondary
 *        variables and, if the a transient problem is modeled, of time information.
 *
 *        This class defines the interface that is expected of variable classes,
 *        and it provides the implementation for models that do not require storing
 *        any additional information besides the primary variables and (optionally)
 *        time.
 * \tparam X The type used for solution vectors, i.e. all primary variables.
 */
template<class X>
class Variables
{
    template<class T, bool indexable>
    struct ScalarExtractorHelper;

    template<class T>
    struct ScalarExtractorHelper<T, false>
    { using Type = T; };

    template<class T>
    struct ScalarExtractorHelper<T, true>
    {
    private:
        using ValueType = std::decay_t<decltype(std::declval<T>()[0])>;
        static constexpr bool indexable = IsIndexable<ValueType>::value;
    public:
        using Type = typename ScalarExtractorHelper<ValueType, indexable>::Type;
    };

public:
    //! export the type of solution vector
    using SolutionVector = X;

    //! export the underlying scalar type
    using Scalar = typename ScalarExtractorHelper<X, IsIndexable<X>::value>::Type;

    //! export the time representation
    using TimeLevel = Dumux::Experimental::TimeLevel<Scalar>;

    //! Default constructor
    explicit Variables() : x_(), t_(0.0) {}

    //! Construction from a solution
    explicit Variables(const SolutionVector& x,
                       const TimeLevel& t = TimeLevel{0.0})
    : x_(x), t_(t)
    {}

    //! Construction from a solution
    explicit Variables(SolutionVector&& x,
                       const TimeLevel& t = TimeLevel{0.0})
    : x_(std::move(x)), t_(t)
    {}

    //! Construction from initializer lambda
    template<class Initializer,
              std::enable_if_t<(std::is_invocable_r_v<void, Initializer, X&>), int> = 0>
    explicit Variables(const Initializer& initializeSolution,
                       const TimeLevel& timeLevel = TimeLevel{0.0})
    : t_(timeLevel)
    {
        initializeSolution(x_);
    }

    //! Return the time level
    const TimeLevel& timeLevel() const
    { return t_; }

    //! Return reference to the solution
    const SolutionVector& dofs() const { return x_; }

    //! Non-const access still required for privar switch (TODO: Remove dependency)
    SolutionVector& dofs() { return x_; }

    //! Update the state to a new solution
    void update(const SolutionVector& x)
    {  x_ = x; }

    //! Update the time level only
    void updateTime(const TimeLevel& t)
    { t_ = t; }

    //! Update the state to a new solution & time level
    void update(const SolutionVector& x,
                const TimeLevel& t)
    {
        x_ = x;
        t_ = t;
    }

private:
    SolutionVector x_;
    TimeLevel t_;
};

} // end namespace Dumux::Experimental

#endif
