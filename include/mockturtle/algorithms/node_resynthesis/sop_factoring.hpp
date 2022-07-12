/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file sop_factoring.hpp
  \brief Resynthesis with SOP factoring

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operators.hpp>
#include <kitty/isop.hpp>

#include "../../utils/stopwatch.hpp"

namespace mockturtle
{

/*! \brief Parameters for sop_factoring function. */
struct sop_factoring_params
{
  /*! \brief Factoring is also tried for the negated TT. */
  bool try_both_polarities{true};
};

/*! \brief Resynthesis function based on SOP factoring.
 *
 * This resynthesis function can be passed to ``node_resynthesis``,
 * ``cut_rewriting``, and ``refactoring``. The method converts a
 * given truth table in an ISOP, then factors the ISOP, and 
 * returns the factored form.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      aig_network aig = ...;

      sop_factoring<aig_network> resyn;
      refactoring( aig, resyn );
   \endverbatim
 *
 */
template<class Ntk>
class sop_factoring
{
public:
  using signal = typename Ntk::signal;
  using sop_t = std::vector<uint64_t>;

public:
  explicit sop_factoring( sop_factoring_params const& ps = {} )
      : _ps( ps ) {}

public:
  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& dest, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    assert( function.num_vars() <= 31 );

    /* derive ISOP */
    bool negated;
    auto cubes = get_isop( function, negated );

    if ( cubes.size() == 0 )
    {
      /* constant 0 */
      fn( dest.get_constant( negated ) );
      return;
    }
    else if ( cubes.size() == 1 && cubes[0]._mask == 0 )
    {
      /* constant 1 */
      fn( dest.get_constant( !negated ) );
      return;
    }

    /* create literal form of SOP */
    sop_t sop = cubes_to_sop( cubes, function.num_vars() );

    /* derive the factored form */
    signal f = gen_factor_rec( dest, {begin, end}, sop, 2 * function.num_vars() );

    fn( negated ? !f : f );

    // auto [and_terms, num_and_gates] = create_function( dest, function, inputs );
    // const auto num_gates = num_and_gates + ( and_terms.empty() ? 0u : static_cast<uint32_t>( and_terms.size() ) - 1u );
    // const auto cand = balanced_tree( dest, and_terms, false );

    // if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
    // {
    //   callback( cand, num_gates );
    // }
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& dest, kitty::dynamic_truth_table const& function, kitty::dynamic_truth_table const& dc, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    assert( function.num_vars() <= 31 );

    /* derive ISOP */
    bool negated;
    auto cubes = get_isop_dc( function, dc, negated );

    if ( cubes.size() == 0 )
    {
      /* constant 0 */
      fn( dest.get_constant( negated ) );
      return;
    }
    else if ( cubes.size() == 1 && cubes[0]._mask == 0 )
    {
      /* constant 1 */
      fn( dest.get_constant( !negated ) );
      return;
    }

    /* create literal form of SOP */
    sop_t sop = cubes_to_sop( cubes, function.num_vars() );

    /* derive the factored form */
    signal f = gen_factor_rec( dest, {begin, end}, sop, 2 * function.num_vars() );

    fn( negated ? !f : f );

    // auto [and_terms, num_and_gates] = create_function( dest, function, inputs );
    // const auto num_gates = num_and_gates + ( and_terms.empty() ? 0u : static_cast<uint32_t>( and_terms.size() ) - 1u );
    // const auto cand = balanced_tree( dest, and_terms, false );

    // if ( cand.level < best_level || ( cand.level == best_level && num_gates < best_cost ) )
    // {
    //   callback( cand, num_gates );
    // }
  }

private:
  std::vector<kitty::cube> get_isop( kitty::dynamic_truth_table const& function, bool& negated ) const
  {
    std::vector<kitty::cube> cubes = kitty::isop( function );

    if ( _ps.try_both_polarities )
    {
      std::vector<kitty::cube> n_cubes = kitty::isop( ~function );

      if ( n_cubes.size() < cubes.size() )
      {
        negated = true;
        return n_cubes; 
      }
      else if ( n_cubes.size() == cubes.size() )
      {
        uint32_t n_lit = 0;
        uint32_t lit = 0;
        for ( auto const& c : n_cubes ) { n_lit += c.num_literals(); }
        for ( auto const& c : cubes ) { lit += c.num_literals(); }

        if ( n_lit < lit )
        {
          negated = true;
          return n_cubes; 
        }
      }
    }

    negated = false;
    return cubes;
  }

  std::vector<kitty::cube> get_isop_dc( kitty::dynamic_truth_table const& function, kitty::dynamic_truth_table const& dc, bool& negated ) const
  {
    std::vector<kitty::cube> cubes;
    kitty::detail::isop_rec( function, function | dc, function.num_vars(), cubes );

    if ( _ps.try_both_polarities )
    {
      std::vector<kitty::cube> n_cubes;
      kitty::detail::isop_rec( ~function, ~function | dc, function.num_vars(), n_cubes );

      if ( n_cubes.size() < cubes.size() )
      {
        negated = true;
        return n_cubes; 
      }
      else if ( n_cubes.size() == cubes.size() )
      {
        uint32_t n_lit = 0;
        uint32_t lit = 0;
        for ( auto const& c : n_cubes ) { n_lit += c.num_literals(); }
        for ( auto const& c : cubes ) { lit += c.num_literals(); }

        if ( n_lit < lit )
        {
          negated = true;
          return n_cubes; 
        }
      }
    }

    negated = false;
    return cubes;
  }

  sop_t cubes_to_sop( std::vector<kitty::cube> const& cubes, uint32_t const num_vars ) const
  {
    sop_t sop( cubes.size() );

    auto it = sop.begin();
    /* represent literals instead of variables ex. a a' b b' */
    /* bit 63 is reserved, up to 31 varibles are supported */
    for ( auto const& c : cubes )
    {
      uint64_t& product = *it++;
      for ( auto i = 0; i < num_vars; ++i )
      {
        if ( c.get_mask( i ) )
          product |= static_cast<uint64_t>( 1 ) << ( 2*i + static_cast<unsigned>( c.get_bit( i ) ) );
      }
    }

    /* TODO: separating positive and negative literals instead of interleaving has better construction performance */
    return sop;
  }

#pragma region SOP factoring
  signal gen_factor_rec( Ntk& ntk, std::vector<signal> const& children, sop_t& sop, uint32_t const num_lit ) const
  {
    sop_t divisor, quotient, reminder;

    assert( sop.size() );

    /* compute the divisor */
    if ( !quick_divisor( sop, divisor, num_lit ) )
    {
      /* generate trivial sop circuit */
      return gen_andor_circuit_rec( ntk, children, sop.begin(), sop.end(), num_lit );
    }

    /* divide the SOP by the divisor */
    divide( sop, divisor, quotient, reminder );

    assert( quotient.size() > 0 );

    if ( quotient.size() == 1 )
    {
      return lit_factor_rec( ntk, children, sop, quotient[0], num_lit );
    }

    make_cube_free( quotient );

    /* divide the SOP by the quotient */
    divide( sop, quotient, divisor, reminder );

    if ( is_cube_free( divisor ) )
    {
      signal div_s = gen_factor_rec( ntk, children, divisor, num_lit );
      signal quot_s = gen_factor_rec( ntk, children, quotient, num_lit );

      /* build (D)*(Q) + R */
      signal dq_and = ntk.create_and( div_s, quot_s );

      if ( reminder.size() )
      {
        signal rem_s = gen_factor_rec( ntk, children, reminder, num_lit );
        return ntk.create_or( dq_and, rem_s );
      }

      return dq_and;
    }

    /* get the common cube */
    uint64_t cube = UINT64_MAX;
    for ( auto const& c : divisor )
    {
      cube &= c;
    }

    return lit_factor_rec( ntk, children, sop, cube, num_lit );
  }

  signal lit_factor_rec( Ntk& ntk, std::vector<signal> const& children, sop_t const& sop, uint64_t const c_sop, uint32_t const num_lit ) const
  {
    sop_t divisor, quotient, reminder;

    /* extract the best literal */
    best_literal( sop, divisor, c_sop, num_lit );

    /* divide SOP by the literal */
    divide_by_cube( sop, divisor, quotient, reminder );

    /* factor the divisor: cube */
    signal div_s = gen_and_circuit_rec( ntk, children, divisor[0], 0, num_lit );

    /* factor the quotient */
    signal quot_s = gen_factor_rec( ntk, children, quotient, num_lit );

    /* build l*Q + R */
    signal dq_and = ntk.create_and( div_s, quot_s );

    /* factor the reminder */
    if ( reminder.size() != 0 )
    {
      signal rem_s = gen_factor_rec( ntk, children, reminder, num_lit );
      return ntk.create_or( dq_and, rem_s );
    }
    
    return dq_and;
  }

  bool quick_divisor( sop_t const& sop, sop_t& res, uint32_t const num_lit ) const
  {
    if ( sop.size() <= 1 )
      return false;

    /* each literal appears no more than once */
    if ( literals_occurrences( sop, num_lit ) < 0 )
      return false;

    /* one level 0-kernel */
    res = sop;
    one_level_zero_kernel_rec( res, num_lit );

    assert( res.size() );
    return true;
  }

  void divide( sop_t& divident, sop_t const& divisor, sop_t& quotient, sop_t& reminder ) const
  {
    /* divisor contains a single cube */
    if ( divisor.size() == 1 )
    {
      divide_by_cube( divident, divisor, quotient, reminder );
      return;
    }

    quotient.clear();
    reminder.clear();

    /* perform division */
    for ( auto i = 0; i < divident.size(); ++i )
    {
      auto const c = divident[i];

      /* cube has been already covered */
      if ( cube_has_lit( c, 63 ) )
        continue;

      uint32_t div_i;
      for ( div_i = 0u; div_i < divisor.size(); ++div_i )
      {
        if ( ( c & divisor[div_i] ) == divisor[div_i] )
          break;
      }

      /* cube is not divisible -> reminder */
      if ( div_i >= divisor.size() )
        continue;

      /* extract quotient */
      uint64_t c_quotient = c & ~divisor[div_i];

      /* find if c_quotient can be obtained for all the divisors */
      bool found = true;
      for ( auto const& div : divisor )
      {
        if ( div == divisor[div_i] )
          continue;

        found = false;
        for ( auto const& c2 : divident )
        {
          /* cube has been already covered */
          if ( cube_has_lit( c2, 63 ) )
            continue;

          if ( ( ( c2 & div ) == div ) && ( c_quotient == ( c2 & ~div ) ) )
          {
            found = true;
            break;
          }
        }

        if ( !found )
          break;
      }

      if ( !found )
        continue;
      
      /* valid divisor, select covered cubes */
      quotient.push_back( c_quotient );

      divident[i] |= static_cast<uint64_t>( 1 ) << 63;
      for ( auto const& div : divisor )
      {
        if ( div == divisor[div_i] )
          continue;

        for ( auto& c2 : divident )
        {
          /* cube has been already covered */
          if ( cube_has_lit( c2, 63 ) )
            continue;

          if ( ( ( c2 & div ) == div ) && ( c_quotient == ( c2 & ~div ) ) )
          {
            c2 |= static_cast<uint64_t>( 1 ) << 63;
            break;
          }
        }
      }
    }

    /* add remainder */
    for ( auto& c : divident )
    {
      if ( !cube_has_lit( c, 63 ) )
      {
        reminder.push_back( c );
      }
      else
      {
        /* unmark */
        c &= ~( static_cast<uint64_t>( 1 ) << 63 );
      }
    }
  }

  int64_t literals_occurrences( sop_t const& sop, uint32_t const num_lit ) const
  {
    /* find the first literal which occurs more than once */
    for ( uint64_t lit = 0; lit < num_lit; ++lit )
    {
      unsigned occurrences = 0;
      for ( auto const& cube : sop )
      {
        if ( cube_has_lit( cube, lit ) )
          ++occurrences;

        if ( occurrences > 1 )
          return lit;
      }
    }

    /* each literal appears once */
    return -1;
  }

  void one_level_zero_kernel_rec( sop_t& sop, uint32_t const num_lit ) const
  {
    /* find least occurring leteral which occurs more than once. TODO: test other metrics */
    int64_t min_lit = least_occurrent_literal( sop, num_lit );

    if ( min_lit == -1 )
      return;
    
    divide_by_literal( sop, static_cast<uint64_t>( min_lit ) );
    make_cube_free( sop );

    one_level_zero_kernel_rec( sop, num_lit );
  }

  int64_t least_occurrent_literal( sop_t const& sop, uint32_t const num_lit ) const
  {
    int64_t min_lit = -1;
    uint32_t min_occurrences = UINT32_MAX;

    /* find the first literal which occurs more than once */
    for ( uint64_t lit = 0; lit < num_lit; ++lit )
    {
      uint32_t occurrences = 0;
      for ( auto const& cube : sop )
      {
        if ( cube_has_lit( cube, lit ) )
          ++occurrences;
      }

      if ( occurrences > 1 && occurrences < min_occurrences )
      {
        min_lit = static_cast<int64_t>( lit );
        min_occurrences = occurrences;
      }
    }

    return min_lit;
  }

  int64_t most_occurrent_literal_masked( sop_t const& sop, uint64_t const cube, uint32_t const num_lit ) const
  {
    int64_t max_lit = -1;
    uint32_t max_occurrences = 1;

    for ( uint64_t lit = 0; lit < num_lit; ++lit )
    {
      if ( !cube_has_lit( cube, lit ) )
        continue;

      uint32_t occurrences = 0;
      for ( auto const& c : sop )
      {
        if ( cube_has_lit( c, lit ) )
          ++occurrences;
      }

      if ( occurrences > max_occurrences )
      {
        max_lit = static_cast<int64_t>( lit );
        max_occurrences = occurrences;
      }
    }

    return max_lit;
  }

  void best_literal( sop_t const& sop, sop_t& result, uint64_t const cube, uint32_t const num_lit ) const
  {
    int64_t max_lit = most_occurrent_literal_masked( sop, cube, num_lit );
    assert( max_lit >= 0 );

    result.push_back( static_cast<uint64_t>( 1 ) << max_lit );
  }

  void divide_by_literal( sop_t& sop, uint64_t const lit ) const
  {
    uint32_t p = 0;
    for ( auto i = 0; i < sop.size(); ++i )
    {
      if ( cube_has_lit( sop[i], lit ) )
      {
        sop[p++] = sop[i] & ( ~( static_cast<uint64_t>( 1 ) << lit ) );
      }
    }

    sop.resize( p );
  }

  void divide_by_cube( sop_t const& divident, sop_t const& divisor, sop_t& quotient, sop_t& reminder ) const
  {
    assert( divisor.size() == 1 );

    quotient.clear();
    reminder.clear();

    uint32_t p = 0;
    for ( auto const& c : divident )
    {
      if ( ( c & divisor[0] ) > 0 )
      {
        quotient.push_back( c & ( ~divisor[0] ) );
      }
      else
      {
        reminder.push_back( c );
      }
    }
  }

  void make_cube_free( sop_t& sop ) const
  {
    /* find common cube */
    uint64_t mask = UINT64_MAX;
    for ( auto const& c : sop )
    {
      mask &= c;
    }

    if ( mask == 0 )
      return;
    
    /* make cube free */
    for ( auto& c : sop )
    {
      c &= ~mask;
    }
  }

  bool is_cube_free( sop_t const& sop ) const
  {
    /* find common cube */
    uint64_t mask = UINT64_MAX;
    for ( auto const& c : sop )
    {
      mask &= c;
    }

    return mask == 0;
  }
#pragma endregion

#pragma region Circuit generation from SOP
  signal gen_and_circuit_rec( Ntk& ntk, std::vector<signal> const& children, uint64_t const cube, uint32_t const begin, uint32_t const end ) const
  {
    /* count set literals */
    uint32_t num_lit = 0;
    uint32_t lit = begin;
    uint32_t i;
    for ( i = begin; i < end; ++i )
    {
      if ( cube_has_lit( cube, i ) )
      {
        ++num_lit;
        lit = i;
      }
    }

    assert( num_lit > 0 );

    if ( num_lit == 1 )
    {
      /* return the coprresponding signal with the correct polarity */
      if ( lit % 2 == 1 )
        return children[lit / 2];
      else
        return !children[lit / 2];
    }

    /* find splitting point */
    uint32_t count_lit = 0;
    for ( i = begin; i < end; ++i )
    {
      if ( cube_has_lit( cube, i ) )
      {
        if ( count_lit >= num_lit / 2 )
          break;

        ++count_lit;
      }
    }

    signal tree1 = gen_and_circuit_rec( ntk, children, cube, begin, i );
    signal tree2 = gen_and_circuit_rec( ntk, children, cube, i, end );

    return ntk.create_and( tree1, tree2 );
  }

  signal gen_andor_circuit_rec( Ntk& ntk, std::vector<signal> const& children, sop_t::const_iterator const& begin, sop_t::const_iterator const& end, uint32_t const num_lit ) const
  {
    auto num_prod = std::distance( begin, end );

    assert( num_prod > 0 );

    if ( num_prod == 1 )
      return gen_and_circuit_rec( ntk, children, *begin, 0, num_lit );

    /* create or tree */
    signal tree1 = gen_andor_circuit_rec( ntk, children, begin, begin + num_prod / 2, num_lit );
    signal tree2 = gen_andor_circuit_rec( ntk, children, begin + num_prod / 2, end, num_lit );

    return ntk.create_or( tree1, tree2 );
  }
#pragma endregion

//   std::pair<arrival_time_queue<Ntk>, uint32_t> create_function( Ntk& dest, kitty::dynamic_truth_table const& func, std::vector<arrival_time_pair<Ntk>> const& arrival_times ) const
//   {
//     const auto sop = create_sop_form( func );

//     stopwatch<> t_tree( time_tree_balancing );
//     arrival_time_queue<Ntk> and_terms;
//     uint32_t num_and_gates{};
//     for ( auto const& cube : sop )
//     {
//       arrival_time_queue<Ntk> product_queue;
//       for ( auto i = 0u; i < func.num_vars(); ++i )
//       {
//         if ( cube.get_mask( i ) )
//         {
//           const auto [f, l] = arrival_times[i];
//           product_queue.push( {cube.get_bit( i ) ? f : dest.create_not( f ), l} );
//         }
//       }
//       if ( !product_queue.empty() )
//       {
//         num_and_gates += static_cast<uint32_t>( product_queue.size() ) - 1u;
//       }
//       and_terms.push( balanced_tree( dest, product_queue ) );
//     }
//     return {and_terms, num_and_gates};
//   }

//   arrival_time_pair<Ntk> balanced_tree( Ntk& dest, arrival_time_queue<Ntk>& queue, bool _and = true ) const
//   {
//     if ( queue.empty() )
//     {
//       return {dest.get_constant( true ), 0u};
//     }

//     while ( queue.size() > 1u )
//     {
//       auto [s1, l1] = queue.top();
//       queue.pop();
//       auto [s2, l2] = queue.top();
//       queue.pop();
//       const auto s = _and ? dest.create_and( s1, s2 ) : dest.create_or( s1, s2 );
//       const auto l = std::max( l1, l2 ) + 1;
//       queue.push( {s, l} );
//     }
//     return queue.top();
//   }

//   std::vector<kitty::cube> create_sop_form( kitty::dynamic_truth_table const& func ) const
//   {
//     stopwatch<> t( time_sop );
//     if ( auto it = sop_hash_.find( func ); it != sop_hash_.end() )
//     {
//       sop_cache_hits++;
//       return it->second;
//     }
//     else
//     {
//       sop_cache_misses++;
//       return sop_hash_[func] = kitty::isop( func ); // TODO generalize
//     }
//   }
private:
  static inline bool cube_has_lit( uint64_t cube, uint64_t lit ) { return ( cube & ( static_cast<uint64_t>( 1 ) << lit ) ) > 0; }

private:
//   mutable std::unordered_map<kitty::dynamic_truth_table, std::vector<kitty::cube>, kitty::hash<kitty::dynamic_truth_table>> sop_hash_;
  sop_factoring_params const& _ps;

public:
  mutable uint32_t sop_cache_hits{};
  mutable uint32_t sop_cache_misses{};

  mutable stopwatch<>::duration time_factoring{};
};

} // namespace mockturtle