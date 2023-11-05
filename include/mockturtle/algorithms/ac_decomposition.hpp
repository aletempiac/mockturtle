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
  \file ac_decomposition.hpp
  \brief Ashenhurst-Curtis decomposition

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <kitty/operations.hpp>
#include <kitty/print.hpp>
#include <kitty/traits.hpp>

namespace mockturtle
{

/*! \brief Parameters for ac_decomposition */
struct ac_decomposition_params
{
  /*! \brief LUT size for decomposition. */
  uint32_t lut_size{ 6 };
};

namespace detail
{

template<typename TT, typename = std::enable_if_t<kitty::is_complete_truth_table<TT>::value>>
class ac_decomposition_impl
{
public:
  ac_decomposition_impl( TT const& tt, uint32_t num_vars, ac_decomposition_params const& ps )
      : tt_start( tt )
      , num_vars( num_vars )
      , ps( ps )
      , permutations( num_vars )
  {}

  // /*! \brief Runs ACD using late arriving variables */
  // void run( std::vector<uint32_t> late_arriving )
  // {

  // }

  // /*! \brief Runs ACD trying different bound sets and free sets */
  // void run()
  // {

  // }

  /*! \brief Runs ACD trying different bound sets */
  uint32_t run_offset( uint32_t free_set_size, uint32_t offset )
  {
    auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, free_set_size ); };

    auto [tt_p, perm] = enumerate_iset_combinations_offset( free_set_size, offset, evaluate_fn, true );

    return column_multiplicity( tt_p, free_set_size );
  }

  /*! \brief Runs ACD trying different bound sets */
  uint32_t run( uint32_t free_set_size )
  {
    auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, free_set_size ); };

    auto [tt_p, perm] = enumerate_iset_combinations( free_set_size, evaluate_fn, false );

    return column_multiplicity( tt_p, free_set_size );
  }

  /*! \brief Runs ACD without trying different bound sets */
  uint32_t run_no_permutations( uint32_t free_set_size )
  {
    uint32_t multiplicity = column_multiplicity( tt_start, free_set_size );

    return multiplicity;
  }
  
private:
  uint32_t column_multiplicity( TT tt, uint32_t free_set_size )
  {
    uint32_t multiplicity_set[4] = { 0u, 0u, 0u, 0u };
    uint32_t multiplicity = 0;

    /* supports up to 64 values of free set (256 for |FS| == 3)*/
    assert( free_set_size <= 3 );

    /* extract iset functions */
    if ( free_set_size == 1 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
      {
        for ( auto j = 0; j < 32; ++j )
        {
          multiplicity_set[0] |= UINT32_C( 1 ) << ( *it & 0x3 );
          *it >>= 2;
        }
        ++it;
      }
    }
    else if ( free_set_size == 2 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
      {
        for ( auto j = 0; j < 16; ++j )
        {
          multiplicity_set[0] |= UINT32_C( 1 ) << ( *it & 0xF );
          *it >>= 4;
        }
        ++it;
      }
    }
    else /* free set size 3 */
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
      {
        for ( auto j = 0; j < 8; ++j )
        {
          multiplicity_set[( *it >> 6 ) & 0x3] |= UINT32_C( 1 ) << ( *it & 0x3F );
          *it >>= 8;
        }
        ++it;
      }
    }

    multiplicity = __builtin_popcountl( multiplicity_set[0] );

    if ( free_set_size == 3 )
    {
      multiplicity += __builtin_popcountl( multiplicity_set[1] );    
      multiplicity += __builtin_popcountl( multiplicity_set[2] );
      multiplicity += __builtin_popcountl( multiplicity_set[3] );
    }

    return multiplicity;
  }

  template<typename Fn>
  std::pair<TT, std::vector<uint32_t>> enumerate_iset_combinations( uint32_t k, Fn&& fn, bool verbose = false )
  {
    TT tt = tt_start;

    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* special case */
    if ( num_vars <= k || k == 0 )
    {
      std::vector<uint32_t> res_perm( num_vars );
      std::iota( res_perm.begin(), res_perm.end(), 0u );
      return { tt, res_perm };
    }

    /* select k */
    k = std::min( k, num_vars - k );

    /* init permutation array */
    std::array<uint32_t, 16> perm, best_perm;
    std::iota( perm.begin(), perm.begin() + num_vars, 0u );
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

      if ( verbose )
      {
        kitty::print_hex( tt );
        std::cout << " " << cost << " ";
        print_perm( perm );
      }

      for ( uint32_t i = 1; i < num_vars; ++i )
      {
        std::swap( perm[0], perm[i] );
        kitty::swap_inplace( tt, 0, i );

        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm );
        }
      }
    }
    else if ( k == 2 )
    {
      for ( uint32_t i = 0; i < num_vars - 1; ++i)
      {
        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm );
        }

        for ( uint32_t j = 2; j < num_vars - i; ++j )
        {
          std::swap( perm[1], perm[j] );
          kitty::swap_inplace( tt, 1, j );

          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best_tt = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm );
          }
        }

        std::swap( perm[0], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, 0, num_vars - i - 1 );
      }
    }
    else if ( k == 3 )
    {
      for ( uint32_t i = 0; i < num_vars - 2; ++i )
      {
        for ( uint32_t j = i; j < num_vars - 2; ++j )
        {
          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best_tt = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm );
          }
  
          for ( uint32_t k = 3; k < num_vars - j; ++k )
          {
            std::swap( perm[2], perm[k] );
            kitty::swap_inplace( tt, 2, k );

            uint32_t cost = fn( tt );
            if ( cost < best_cost )
            {
              best_tt = tt;
              best_cost = cost;
              best_perm = perm;
            }

            if ( verbose )
            {
              kitty::print_hex( tt );
              std::cout << " " << cost << " ";
              print_perm( perm );
            }
          }

          std::swap( perm[1], perm[num_vars - j - 1] );
          kitty::swap_inplace( tt, 1, num_vars - j - 1 );
        }

        std::swap( perm[0], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, 0, num_vars - i - 1 );
      }
    }

    std::vector<uint32_t> res_perm( num_vars );
    std::copy( best_perm.begin(), best_perm.begin() + num_vars, res_perm.begin() );

    return { best_tt, res_perm };
  }

  template<typename Fn>
  std::pair<TT, std::vector<uint32_t>> enumerate_iset_combinations_offset( uint32_t k, uint32_t offset, Fn&& fn, bool verbose = false )
  {
    TT tt = tt_start;

    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* select k */
    k = std::min( k, num_vars - k );
    k -= offset;

    /* special case */
    if ( num_vars <= k || k == 0 )
    {
      std::vector<uint32_t> res_perm( num_vars );
      std::iota( res_perm.begin(), res_perm.end(), 0u );
      return { tt, res_perm };
    }

    /* init permutation array */
    std::array<uint32_t, 16> perm, best_perm;
    std::iota( perm.begin(), perm.begin() + num_vars, 0u );
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

      if ( verbose )
      {
        kitty::print_hex( tt );
        std::cout << " " << cost << " ";
        print_perm( perm );
      }

      for ( uint32_t i = offset + 1; i < num_vars; ++i )
      {
        std::swap( perm[offset], perm[i] );
        kitty::swap_inplace( tt, offset, i );

        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm );
        }
      }
    }
    else if ( k == 2 )
    {
      for ( uint32_t i = 0; i < num_vars - 1 - offset; ++i)
      {
        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best_tt = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm );
        }

        for ( uint32_t j = offset + 2; j < num_vars - i; ++j )
        {
          std::swap( perm[offset + 1], perm[j] );
          kitty::swap_inplace( tt, offset + 1, j );

          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best_tt = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm );
          }
        }

        std::swap( perm[offset], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, offset, num_vars - i - 1 );
      }
    }
    else if ( k == 3 )
    {
      for ( uint32_t i = 0; i < num_vars - 2 - offset; ++i )
      {
        for ( uint32_t j = i; j < num_vars - 2 - offset; ++j )
        {
          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best_tt = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm );
          }
  
          for ( uint32_t k = offset + 3; k < num_vars - j; ++k )
          {
            std::swap( perm[offset + 2], perm[k] );
            kitty::swap_inplace( tt, offset + 2, k );

            uint32_t cost = fn( tt );
            if ( cost < best_cost )
            {
              best_tt = tt;
              best_cost = cost;
              best_perm = perm;
            }

            if ( verbose )
            {
              kitty::print_hex( tt );
              std::cout << " " << cost << " ";
              print_perm( perm );
            }
          }

          std::swap( perm[offset + 1], perm[num_vars - j - 1] );
          kitty::swap_inplace( tt, offset + 1, num_vars - j - 1 );
        }

        std::swap( perm[offset], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, offset, num_vars - i - 1 );
      }
    }

    std::vector<uint32_t> res_perm( num_vars );
    std::copy( best_perm.begin(), best_perm.begin() + num_vars, res_perm.begin() );

    return { best_tt, res_perm };
  }

  void print_perm( std::array<uint32_t, 16>& perm )
  {
    std::cout << "[";
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      std::cout << perm.at( i ) << " ";
    }
    std::cout << "]\n";
  }

private:
  TT const& tt_start;
  uint32_t num_vars;
  ac_decomposition_params const& ps;
  std::vector<uint32_t> permutations;
};

} // namespace detail

} // namespace mockturtle