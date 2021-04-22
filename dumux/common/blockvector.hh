// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \brief docme
 */
#ifndef DUMUX_COMMON_BLOCKVECTOR_HH
#define DUMUX_COMMON_BLOCKVECTOR_HH

// // forward declaration
// namespace Dune {
//
// template<class... Args>
// class MultiTypeBlockVector;
//
// } // end namespace Dune

namespace Dumux::Istl {

template<class BVType>
class BlockVector
{
public:
    using field_type = typename BVType::field_type;
    using block_type = typename BVType::block_type;
    using allocator_type = typename BVType::allocator_type;
    using size_type = typename BVType::size_type;
    using Iterator = typename BVType::Iterator;
    using ConstIterator = typename BVType::ConstIterator;

    BlockVector() = default;

    //! make vector with _n components
    explicit BlockVector (size_type n) : blockVector_(n)
    {}

    /** \brief Construct from a std::initializer_list */
    BlockVector (const std::initializer_list<block_type>& l) : blockVector_(l)
    {}

    template<typename S>
    BlockVector (size_type n, S capacity) : blockVector_(n, capacity)
    {}

    block_type& operator [](int i)
    { return blockVector_[i]; }

    const block_type& operator[] (int i) const
    { return blockVector_[i]; }

    BlockVector& operator= (const field_type& k)
    {
        blockVector_ = k;
        return *this;
    }

    BlockVector& operator*= (const field_type& k)
    {
        blockVector_ *= k;
        return *this;
    }

    BlockVector& operator/= (const field_type& k)
    {
        blockVector_ /= k;
        return *this;
    }

    BlockVector& operator+= (const BlockVector& y)
    {
        blockVector_ += y.blockVector_;
        return *this;
    }

    BlockVector& operator-= (const BlockVector& y)
    {
        blockVector_ -= y.blockVector_;
        return *this;
    }

    BlockVector operator* (const BlockVector& y)
    {
        return blockVector_ * y.blockVector_;
    }

    BlockVector& axpy(const field_type& a, const BlockVector& y)
    {
        blockVector_.axpy(a, y.blockVector_);
        return *this;
    }

    BlockVector& dot(const BlockVector& y)
    {
        blockVector_.dot(y.blockVector_);
        return *this;
    }

    auto one_norm () const
    { return blockVector_.one_norm(); }

    auto two_norm () const
    { return blockVector_.two_norm(); }

    auto two_norm2 () const
    { return blockVector_.two_norm2(); }

    auto infinity_norm () const
    { return blockVector_.infinity_norm(); }

    auto begin() const
    { return blockVector_.begin(); }

    auto end() const
    { return blockVector_.end(); }

    void reserve(size_type capacity)
    { blockVector_.reserve(capacity); }

    size_type capacity() const
    { return blockVector_.capacity(); }

    size_type size() const
    { return blockVector_.size(); }

    void resize(size_type size)
    { blockVector_.resize(size); }

    BVType& native()
    { return blockVector_; }

    const BVType& native() const
    { return blockVector_; }

private:
    BVType blockVector_;
};

template<class BVType, class State>
class BlockVectorWithState
{
    struct BlockVectorView
    {
        BlockVectorView(BVType& sol, State& state)
        {
            solution_ = &sol;
            state_ = &state;
        }

       auto& operator= (const BVType& other)
       {
           *solution_ = other;
           return *this;
       }

       void setState(int s)
       { *state_ = s; }

       operator BVType()
       { return *solution_; }

       State state() const
       { return *state_; }

   private:
       BVType* solution_;
       State* state_;
    };

    struct ConstBlockVectorView
    {
        ConstBlockVectorView(const BVType& sol, const State& state) : solution_(&sol), state_(&state)
        {}

       operator BVType() const
       { return *solution_; }

       State state() const
       { return *state_; }

   private:
       const BVType* solution_;
       const  State* state_;
    };

public:
    using field_type = typename BVType::field_type;
    using block_type = typename BVType::block_type;
    using allocator_type = typename BVType::allocator_type;
    using size_type = typename BVType::size_type;
    using Iterator = typename BVType::Iterator;
    using ConstIterator = typename BVType::ConstIterator;

    BlockVectorWithState() = default;

    //! make vector with _n components
    explicit BlockVectorWithState (size_type n) : blockVector_(n)
    {
        states_.resize(blockVector_.size());
    }

    /** \brief Construct from a std::initializer_list */
    BlockVectorWithState (const std::initializer_list<block_type>& l) : blockVector_(l)
    {
        states_.resize(blockVector_.size());
    }

    template<typename S>
    BlockVectorWithState (size_type n, S capacity) : blockVector_(n, capacity)
    {
        states_.resize(blockVector_.size());
    }

    block_type& operator [](int i)
    {
        return BlockVectorView(blockVector_[dofIdx], states_[dofIdx]);
    }

    const block_type& operator[] (int i) const
    {
        return ConstBlockVectorView(blockVector_[dofIdx], states_[dofIdx]);
    }

    BlockVectorWithState& operator= (const field_type& k)
    {
        blockVector_ = k;
        return *this;
    }

    BlockVectorWithState& operator*= (const field_type& k)
    {
        blockVector_ *= k;
        return *this;
    }

    BlockVectorWithState& operator/= (const field_type& k)
    {
        blockVector_ /= k;
        return *this;
    }

    BlockVectorWithState& operator+= (const BlockVectorWithState& y)
    {
        blockVector_ += y.blockVector_;
        return *this;
    }

    BlockVectorWithState& operator-= (const BlockVectorWithState& y)
    {
        blockVector_ -= y.blockVector_;
        return *this;
    }

    BlockVectorWithState operator* (const BlockVectorWithState& y)
    {
        return blockVector_ * y.blockVector_;
    }

    BlockVectorWithState& axpy(const field_type& a, const BlockVectorWithState& y)
    {
        blockVector_.axpy(a, y.blockVector_);
        return *this;
    }

    BlockVectorWithState& dot(const BlockVectorWithState& y)
    {
        blockVector_.dot(y.blockVector_);
        return *this;
    }

    auto one_norm () const
    { return blockVector_.one_norm(); }

    auto two_norm () const
    { return blockVector_.two_norm(); }

    auto two_norm2 () const
    { return blockVector_.two_norm2(); }

    auto infinity_norm () const
    { return blockVector_.infinity_norm(); }

    auto begin() const
    { return blockVector_.begin(); }

    auto end() const
    { return blockVector_.end(); }

    void reserve(size_type capacity)
    {
        blockVector_.reserve(capacity);
        states_.reserve(capacity);
    }

    size_type capacity() const
    { return blockVector_.capacity(); }

    size_type size() const
    { return blockVector_.size(); }

    void resize(size_type size)
    {
        blockVector_.resize(size);
        states_.resize(size);
    }

    BVType& native()
    { return blockVector_; }

    const BVType& native() const
    { return blockVector_; }

private:
    BVType blockVector_;
    States states_;
};

} // end namespace Dumux::Istl

#endif
