/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file decomposition_utils.hpp
  \brief Utilities for functional decomposition

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <kitty/traits.hpp>

namespace mockturtle
{

namespace detail
{

inline void print_perm( std::array<uint32_t, 16>& perm, uint32_t n )
{
  std::cout << "[";
  for ( uint32_t i = 0; i < n; ++i )
  {
    std::cout << perm.at( i ) << " ";
  }
  std::cout << "]\n";
}

template<typename TT, typename Fn, typename = std::enable_if_t<kitty::is_complete_truth_table<TT>::value>>
inline std::pair<TT, std::vector<uint32_t>> enumerate_iset_combinations( TT tt, uint32_t k, Fn&& fn )
{
  uint32_t n = tt.num_vars();

  /* works up to 16 input truth tables */
  assert( n <= 16 );

  /* special case */
  if ( n <= k )
  {
    std::vector<uint32_t> res_perm( n );
    std::iota( res_perm.begin(), res_perm.end(), 0u );
    return { tt, res_perm };
  }

  /* select k */
  k = std::min( k, n - k );

  /* init permutation array */
  std::array<uint32_t, 16> perm, best_perm;
  std::iota( perm.begin(), perm.begin() + n, 0u );
  best_perm = perm;

  /* TT with best cost */
  TT best_tt = tt;
  uint32_t best_cost = UINT32_MAX;

  /* enumerate combinations */
  if ( k == 1 )
  {
    uint32_t cost = fn( tt );
    if ( cost < best_cost )
    {
      best_tt = tt;
      best_cost = cost;
      best_perm = perm;
    }

    for ( uint32_t i = 2; i <= n; ++i )
    {
      std::swap( perm[n - 1], perm[n - i] );
      kitty::swap_inplace( tt, n - 1, n - i );

      uint32_t cost = fn( tt );
      if ( cost < best_cost )
      {
        best_tt = tt;
        best_cost = cost;
        best_perm = perm;
      }
    }
  }
  else if ( k == 2 )
  {
    for ( uint32_t i = 0; i < n - 1; ++i)
    {
      uint32_t cost = fn( tt );
      if ( cost < best_cost )
      {
        best_tt = tt;
        best_cost = cost;
        best_perm = perm;
      }

      for ( uint32_t j = 3; j <= n - i; ++j )
      {
        std::swap( perm[n - 2], perm[n - j] );
        kitty::swap_inplace( tt, n - 2, n - j );

        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }
      }

      std::swap( perm[n - 1], perm[i] );
      kitty::swap_inplace( tt, n - 1, i );
    }
  }
  else if ( k == 3 )
  {
    for ( uint32_t i = 0; i < n - 2; ++i )
    {
      for ( uint32_t j = i; j < n - 2; ++j )
      {
        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }
        print_perm( perm, n );
 
        for ( uint32_t k = 4; k <= n - j; ++k )
        {
          std::swap( perm[n - 3], perm[n - k] );
          kitty::swap_inplace( tt, n - 3, n - k );

          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best_tt = tt;
            best_cost = cost;
            best_perm = perm;
          }

          print_perm( perm, n );
        }

        std::swap( perm[n - 2], perm[j] );
        kitty::swap_inplace( tt, n - 2, j );
      }

      std::swap( perm[n - 1], perm[i] );
      kitty::swap_inplace( tt, n - 1, i );
    }
  }

  std::vector<uint32_t> res_perm( n );
  std::copy( best_perm.begin(), best_perm.begin() + n, res_perm.begin() );
  
  return { tt, res_perm };
}

} // namespace detail

template<typename TT, typename = std::enable_if_t<kitty::is_complete_truth_table<TT>::value>>
inline uint32_t acd_column_multiplicity( TT tt, uint32_t free_set_size )
{
  uint32_t multiplicity_set = 0;

  /* supports up to 64 values of free set (16 for |FS| == 2)*/
  assert( free_set_size <= 2 );

  /* extract iset functions */
  if ( free_set_size == 1 )
  {
    auto it = std::begin( tt );
    for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
    {
      for ( auto j = 0; j < 32; ++j )
      {
        multiplicity_set |= UINT32_C( 1 ) << ( *it & 0x3 );
        *it >>= 2;
      }
      ++it;
    }
  }
  else /* free set size 2 */
  {
    auto it = std::begin( tt );
    for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
    {
      for ( auto j = 0; j < 16; ++j )
      {
        multiplicity_set |= UINT32_C( 1 ) << ( *it & 0xF );
        *it >>= 4;
      }
      ++it;
    }
  }

  return __builtin_popcount( multiplicity_set );
}

template<typename TT, typename = std::enable_if_t<kitty::is_complete_truth_table<TT>::value>>
inline void acd_enumerate_combinations( TT tt, uint32_t free_set_size )
{
  auto evaluate_fn = [&] ( TT const& tt ) { return acd_column_multiplicity( tt, free_set_size ); };

  detail::enumerate_iset_combinations( tt, free_set_size, [&] ( TT const& tt ) { return 0; } );
}

} // namespace mockturtle
