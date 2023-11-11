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

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/detail/constants.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>
#include <kitty/print.hpp>
#include <kitty/traits.hpp>

#include "../networks/klut.hpp"
#include "simulation.hpp"

namespace mockturtle
{

/*! \brief Parameters for ac_decomposition */
struct ac_decomposition_params
{
  /*! \brief LUT size for decomposition. */
  uint32_t lut_size{ 6 };
};

/*! \brief Statistics for ac_decomposition */
struct ac_decomposition_stats
{
  uint32_t num_luts{ 0 };
  uint32_t num_edges{ 0 };
  uint32_t num_levels{ 0 };

  // void report()
  // {
  //   std::cout << "";
  // }
};

struct ac_decomposition_result
{
  kitty::dynamic_truth_table tt;
  std::vector<uint32_t> support;
};

namespace detail
{

template<typename TT, typename = std::enable_if_t<kitty::is_complete_truth_table<TT>::value>>
class ac_decomposition_impl
{
private:
  struct encoding_matrix
  {
    uint64_t column{ 0 };
    uint32_t cost{ 0 };
    uint32_t index{ 0 };
    uint32_t sort_cost{ 0 };
  };

private:
  static constexpr uint32_t max_num_vars = 8;
  using STT = kitty::static_truth_table<max_num_vars>;

public:
  ac_decomposition_impl( TT const& tt, uint32_t num_vars, ac_decomposition_params const& ps, ac_decomposition_stats* pst = nullptr )
      : num_vars( num_vars ), ps( ps ), pst( pst ), permutations( num_vars )
  {
    tt_start = tt;
    std::iota( permutations.begin(), permutations.end(), 0 );
  }

  /*! \brief Runs ACD using late arriving variables */
  uint32_t run( std::vector<uint32_t> late_arriving )
  {
    /* truth table is too large for the settings */
    if ( num_vars > max_num_vars )
    {
      return UINT32_MAX;
    }

    /* return a high cost if too many late arriving variables */
    if ( late_arriving.size() > ps.lut_size / 2 || late_arriving.size() > 3 )
    {
      return UINT32_MAX;
    }

    /* convert to static TT */
    best_tt = kitty::extend_to<max_num_vars>( tt_start );
    best_multiplicity = UINT32_MAX;

    /* permute late arriving variables to be the least significant */
    reposition_late_arriving_variables( late_arriving );

    /* run ACD trying different bound sets and free sets */
    uint32_t free_set_size = late_arriving.size();
    uint32_t offset = std::max( static_cast<uint32_t>( late_arriving.size() ), 1u );
    for ( uint32_t i = offset; i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&]( STT const& tt ) { return column_multiplicity( tt, i ); };
      auto [tt_p, perm, cost] = enumerate_iset_combinations_offset( i, offset, evaluate_fn, false );

      /* check for feasible solution that improves the cost */
      if ( cost <= ( 1 << i ) && cost < best_multiplicity )
      {
        best_tt = tt_p;
        permutations = perm;
        best_multiplicity = cost;
        free_set_size = i;
      }
    }

    if ( best_multiplicity == UINT32_MAX )
      return UINT32_MAX;

    /* compute isets */
    std::vector<STT> isets = compute_isets( free_set_size );

    /* test for column multiplicity 4 */
    // if ( best_multiplicity == 4 )
    // {
    //   test_support_minimization_isets( isets, true );
    // }
    // else

    generate_support_minimization_encodings();
    solve_min_support_exact( isets, free_set_size );

    /* unfeasible decomposition */
    if ( best_bound_sets.empty() )
    {
      return UINT32_MAX;
    }

    /* TODO generate decomposition only when returning the result */
    dec_result = generate_decomposition( free_set_size );

    /* TODO: change return value */
    return best_multiplicity;
  }

  /*! \brief Runs ACD using late arriving variables and guaranteeing support minimization */
  uint32_t run_dsd( std::vector<uint32_t> late_arriving )
  {
    /* truth table is too large for the settings */
    if ( num_vars > max_num_vars )
    {
      return UINT32_MAX;
    }

    /* compute minimum number of variables in the free set */
    uint32_t dsd_vars = num_vars - ps.lut_size;

    /* return a high cost if too many late arriving variables */
    if ( late_arriving.size() > ps.lut_size / 2 || late_arriving.size() > 3 || dsd_vars > 3 )
    {
      return UINT32_MAX;
    }

    /* convert to static TT */
    best_tt = kitty::extend_to<max_num_vars>( tt_start );
    best_multiplicity = UINT32_MAX;

    /* permute late arriving variables to be the least significant */
    reposition_late_arriving_variables( late_arriving );

    /* run ACD trying different bound sets and free sets */
    uint32_t free_set_size = late_arriving.size();
    uint32_t offset = std::max( static_cast<uint32_t>( late_arriving.size() ), 1u );
    for ( uint32_t i = std::max( dsd_vars, offset ); i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&]( STT const& tt ) { return column_multiplicity( tt, i ); };
      auto [tt_p, perm, cost] = enumerate_iset_combinations_offset( i, offset, evaluate_fn, false );

      /* check for feasible solution that improves the cost */
      if ( cost <= ( 1 << i ) && cost < best_multiplicity )
      {
        best_tt = tt_p;
        permutations = perm;
        best_multiplicity = cost;
        free_set_size = i;
      }
    }

    /* TODO: change return value */
    // return ps.lut_size - free_set_size + 1;
    return best_multiplicity;
  }

  /*! \brief Runs ACD trying different bound sets and free sets */
  uint32_t run()
  {
    /* truth table is too large for the settings */
    if ( num_vars > max_num_vars )
    {
      return UINT32_MAX;
    }

    /* convert to static TT */
    best_tt = kitty::extend_to<max_num_vars>( tt_start );
    best_multiplicity = UINT32_MAX;

    /* run ACD trying different bound sets and free sets */
    uint32_t free_set_size = 1;
    for ( uint32_t i = 1; i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&]( STT const& tt ) { return column_multiplicity( tt, i ); };
      auto [tt_p, perm, cost] = enumerate_iset_combinations( i, evaluate_fn, false );

      /* check for feasible solution that improves the cost */
      if ( cost <= ( 1 << i ) && cost < best_multiplicity )
      {
        best_tt = tt_p;
        permutations = perm;
        best_multiplicity = cost;
        free_set_size = i;
      }
    }

    if ( best_multiplicity == UINT32_MAX )
      return UINT32_MAX;

    /* compute isets */
    std::vector<STT> isets = compute_isets( free_set_size );

    /* test for column multiplicity 4 */
    // if ( best_multiplicity == 4 )
    // {
    //   test_support_minimization_isets( isets, true );
    // }
    // else

    generate_support_minimization_encodings();
    solve_min_support_exact( isets, free_set_size );

    /* unfeasible decomposition */
    if ( best_bound_sets.empty() )
    {
      return UINT32_MAX;
    }

    /* TODO generate decomposition only when returning the result */
    dec_result = generate_decomposition( free_set_size );

    /* TODO: change return value */
    return best_multiplicity;
  }

  /*! \brief Runs ACD trying different bound sets */
  uint32_t run( uint32_t free_set_size )
  {
    /* truth table is too large for the settings */
    if ( num_vars > max_num_vars || free_set_size > ps.lut_size / 2 )
    {
      return UINT32_MAX;
    }

    /* convert to static TT */
    best_tt = kitty::extend_to<max_num_vars>( tt_start );
    best_multiplicity = UINT32_MAX;

    auto evaluate_fn = [&]( STT const& tt ) { return column_multiplicity( tt, free_set_size ); };
    auto [tt_p, perm, cost] = enumerate_iset_combinations( free_set_size, evaluate_fn, false );
    best_tt = tt_p;
    permutations = perm;
    best_multiplicity = cost;

    if ( best_multiplicity == UINT32_MAX )
      return UINT32_MAX;

    /* compute isets */
    std::vector<STT> isets = compute_isets( free_set_size );

    generate_support_minimization_encodings();
    solve_min_support_exact( isets, free_set_size );

    /* unfeasible decomposition */
    if ( best_bound_sets.empty() )
    {
      return UINT32_MAX;
    }

    /* TODO generate decomposition only when returning the result */
    dec_result = generate_decomposition( free_set_size );

    /* TODO: change return value */
    return best_multiplicity;
  }

  std::vector<ac_decomposition_result> get_result()
  {
    return dec_result;
  }

  std::optional<klut_network> get_result_ntk()
  {
    if ( dec_result.empty() )
      return std::nullopt;

    return get_result_ntk_impl();
  }

  bool verify_equivalence()
  {
    if ( dec_result.empty() )
      return false;

    return verify_equivalence_impl();
  }

private:
  uint32_t column_multiplicity( STT tt, uint32_t free_set_size )
  {
    uint64_t multiplicity_set[4] = { 0u, 0u, 0u, 0u };
    uint32_t multiplicity = 0;
    uint32_t num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;

    /* supports up to 64 values of free set (256 for |FS| == 3)*/
    assert( free_set_size <= 3 );

    /* extract iset functions */
    if ( free_set_size == 1 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 32; ++j )
        {
          multiplicity_set[0] |= UINT64_C( 1 ) << ( *it & 0x3 );
          *it >>= 2;
        }
        ++it;
      }
    }
    else if ( free_set_size == 2 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 16; ++j )
        {
          multiplicity_set[0] |= UINT64_C( 1 ) << ( *it & 0xF );
          *it >>= 4;
        }
        ++it;
      }
    }
    else /* free set size 3 */
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 8; ++j )
        {
          multiplicity_set[( *it >> 6 ) & 0x3] |= UINT64_C( 1 ) << ( *it & 0x3F );
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
  std::tuple<STT, std::vector<uint32_t>, uint32_t> enumerate_iset_combinations( uint32_t free_set_size, Fn&& fn, bool verbose = false )
  {
    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* special case */
    STT tt = best_tt;
    if ( num_vars <= free_set_size || free_set_size == 0 )
    {
      return { tt, permutations, UINT32_MAX };
    }

    /* select k */
    // free_set_size = std::min( free_set_size, num_vars - free_set_size );

    /* init permutation array */
    std::array<uint32_t, 16> perm, best_perm;
    std::copy( permutations.begin(), permutations.begin() + num_vars, perm.begin() );
    best_perm = perm;

    /* TT with best cost */
    STT best = tt;
    uint32_t best_cost = UINT32_MAX;

    /* enumerate combinations */
    if ( free_set_size == 1 )
    {
      uint32_t cost = fn( tt );
      if ( cost < best_cost )
      {
        best = tt;
        best_cost = cost;
        best_perm = perm;
      }

      if ( verbose )
      {
        kitty::print_hex( tt );
        std::cout << " " << cost << " ";
        print_perm( perm.begin(), free_set_size );
      }

      for ( uint32_t i = 1; i < num_vars; ++i )
      {
        std::swap( perm[0], perm[i] );
        kitty::swap_inplace( tt, 0, i );

        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm.begin(), free_set_size );
        }
      }
    }
    else if ( free_set_size == 2 )
    {
      for ( uint32_t i = 0; i < num_vars - 1; ++i )
      {
        uint32_t cost = fn( tt );
        if ( cost < best_cost )
        {
          best = tt;
          best_cost = cost;
          best_perm = perm;
        }

        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << cost << " ";
          print_perm( perm.begin(), free_set_size );
        }

        for ( uint32_t j = 2; j < num_vars - i; ++j )
        {
          std::swap( perm[1], perm[j] );
          kitty::swap_inplace( tt, 1, j );

          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm.begin(), free_set_size );
          }
        }

        std::swap( perm[0], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, 0, num_vars - i - 1 );
      }
    }
    else if ( free_set_size == 3 )
    {
      for ( uint32_t i = 0; i < num_vars - 2; ++i )
      {
        for ( uint32_t j = i; j < num_vars - 2; ++j )
        {
          uint32_t cost = fn( tt );
          if ( cost < best_cost )
          {
            best = tt;
            best_cost = cost;
            best_perm = perm;
          }

          if ( verbose )
          {
            kitty::print_hex( tt );
            std::cout << " " << cost << " ";
            print_perm( perm.begin(), free_set_size );
          }

          for ( uint32_t k = 3; k < num_vars - j; ++k )
          {
            std::swap( perm[2], perm[k] );
            kitty::swap_inplace( tt, 2, k );

            uint32_t cost = fn( tt );
            if ( cost < best_cost )
            {
              best = tt;
              best_cost = cost;
              best_perm = perm;
            }

            if ( verbose )
            {
              kitty::print_hex( tt );
              std::cout << " " << cost << " ";
              print_perm( perm.begin(), free_set_size );
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

    return std::make_tuple( best, res_perm, best_cost );
  }

  template<typename Fn>
  std::tuple<STT, std::vector<uint32_t>, uint32_t> enumerate_iset_combinations_offset( uint32_t free_set_size, uint32_t offset, Fn&& fn, bool verbose = false )
  {
    STT tt = best_tt;

    /* TT with best cost */
    STT best_tt = tt;
    uint32_t best_cost = UINT32_MAX;

    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* select k */
    free_set_size = std::min( free_set_size, num_vars - free_set_size );

    /* special case */
    if ( num_vars <= free_set_size || free_set_size <= offset )
    {
      if ( offset == free_set_size )
      {
        best_cost = fn( tt );
        if ( verbose )
        {
          kitty::print_hex( tt );
          std::cout << " " << best_cost << " ";
          print_perm( permutations.begin(), free_set_size );
        }

        return { tt, permutations, best_cost };
      }
      else
      {
        return { tt, permutations, UINT32_MAX };
      }
    }

    /* decrease combinations */
    free_set_size -= offset;

    /* init permutation array */
    std::array<uint32_t, 16> perm, best_perm;
    std::copy( permutations.begin(), permutations.begin() + num_vars, perm.begin() );
    best_perm = perm;

    /* enumerate combinations */
    if ( free_set_size == 1 )
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
        print_perm( perm.begin(), free_set_size + offset );
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
          print_perm( perm.begin(), free_set_size + offset );
        }
      }
    }
    else if ( free_set_size == 2 )
    {
      for ( uint32_t i = 0; i < num_vars - 1 - offset; ++i )
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
          print_perm( perm.begin(), free_set_size + offset );
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
            print_perm( perm.begin(), free_set_size + offset );
          }
        }

        std::swap( perm[offset], perm[num_vars - i - 1] );
        kitty::swap_inplace( tt, offset, num_vars - i - 1 );
      }
    }
    else if ( free_set_size == 3 )
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
            print_perm( perm.begin(), free_set_size + offset );
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
              print_perm( perm.begin(), free_set_size + offset );
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

    return std::make_tuple( best_tt, res_perm, best_cost );
  }

  std::vector<STT> compute_isets( uint32_t free_set_size, bool verbose = false )
  {
    /* construct isets involved in multiplicity */
    uint32_t isets_support = num_vars - free_set_size;
    std::vector<STT> isets( best_multiplicity );

    /* construct isets */
    std::unordered_map<uint64_t, uint32_t> column_to_iset;
    STT tt = best_tt;
    uint32_t offset = 0, block = 0;
    uint32_t num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;

    if ( free_set_size == 1 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 32; ++j )
        {
          uint64_t val = *it & 0x3;

          if ( auto el = column_to_iset.find( val ); el != column_to_iset.end() )
          {
            isets[el->second]._bits[i / 2] |= UINT64_C( 1 ) << ( j + offset );
          }
          else
          {
            isets[column_to_iset.size()]._bits[i / 2] |= UINT64_C( 1 ) << ( j + offset );
            column_to_iset[val] = column_to_iset.size();
          }

          *it >>= 2;
        }

        offset ^= 32;
        ++it;
      }
    }
    else if ( free_set_size == 2 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 16; ++j )
        {
          uint64_t val = *it & 0xF;

          if ( auto el = column_to_iset.find( val ); el != column_to_iset.end() )
          {
            isets[el->second]._bits[i / 4] |= UINT64_C( 1 ) << ( j + offset );
          }
          else
          {
            isets[column_to_iset.size()]._bits[i / 4] |= UINT64_C( 1 ) << ( j + offset );
            column_to_iset[val] = column_to_iset.size();
          }

          *it >>= 4;
        }

        offset = ( offset + 16 ) % 64;
        ++it;
      }
    }
    else /* free set size 3 */
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < num_blocks; ++i )
      {
        for ( auto j = 0; j < 8; ++j )
        {
          uint64_t val = *it & 0xFF;

          if ( auto el = column_to_iset.find( val ); el != column_to_iset.end() )
          {
            isets[el->second]._bits[i / 8] |= UINT64_C( 1 ) << ( j + offset );
          }
          else
          {
            isets[column_to_iset.size()]._bits[i / 8] |= UINT64_C( 1 ) << ( j + offset );
            column_to_iset[val] = column_to_iset.size();
          }

          *it >>= 8;
        }

        offset = ( offset + 8 ) % 64;
        ++it;
      }
    }

    /* extend isets to cover the whole truth table */
    for ( STT& iset : isets )
    {
      local_extend_to( iset, isets_support );
    }

    /* save free_set functions */
    std::vector<STT> free_set_tts( best_multiplicity );

    /* TODO: possible conflict */
    for ( auto const& pair : column_to_iset )
    {
      free_set_tts[pair.second]._bits[0] = pair.first;
      local_extend_to( free_set_tts[pair.second], free_set_size );
    }

    best_free_set_tts = std::move( free_set_tts );

    /* print isets  and free set*/
    if ( verbose )
    {
      std::cout << "iSets\n";
      uint32_t i = 0;
      for ( auto iset : isets )
      {
        kitty::print_hex( iset );
        std::cout << " of func ";
        kitty::print_hex( free_set_tts[i++] );
        std::cout << "\n";
      }
    }

    return isets;
  }

  void test_support_minimization_isets( std::vector<kitty::dynamic_truth_table> const& isets, bool verbose = false )
  {
    assert( best_multiplicity == 4 );

    std::vector<kitty::dynamic_truth_table> bound_sets( 2 );
    uint32_t best_cost_luts = UINT32_MAX;
    uint32_t best_cost_edges = UINT32_MAX;

    /* reset bound set values */
    best_bound_sets.clear();

    /* enumerate combinations */
    for ( int32_t i = 0; i < 6; ++i )
    {
      /* compute bound set */
      bound_sets[0] = isets[iset4_combinations[i][0][0]] | isets[iset4_combinations[i][0][1]];
      bound_sets[1] = isets[iset4_combinations[i][1][0]] | isets[iset4_combinations[i][1][1]];

      /* check support minimization */
      uint32_t vars0 = 0;
      uint32_t vars1 = 0;
      for ( uint32_t j = 0; j < isets[0].num_vars(); ++j )
      {
        vars0 += kitty::has_var( bound_sets[0], j ) ? 1 : 0;
        vars1 += kitty::has_var( bound_sets[1], j ) ? 1 : 0;
      }

      /* check cost */
      if ( vars0 > ps.lut_size || vars1 > ps.lut_size )
      {
        /* not feasible */
        continue;
      }

      uint32_t cost_luts = vars0 == 1 ? 0 : 1;
      cost_luts += vars1 == 1 ? 0 : 1;
      uint32_t cost_edges = vars0 == 1 ? 0 : vars0;
      cost_edges += vars1 == 1 ? 0 : vars1;

      if ( cost_luts < best_cost_luts || cost_luts == best_cost_luts && cost_edges < best_cost_edges )
      {
        best_bound_sets = bound_sets;
        best_cost_luts = cost_luts;
        best_cost_edges = cost_edges;

        /* load ONSET and OFFSET */
        best_care_sets.clear();
        best_iset_onset.clear();
        best_iset_offset.clear();
        for ( uint32_t k = 0; k < bound_sets.size(); ++k )
        {
          /* TODO: fix (currently considers only 2 members) */
          uint64_t onset = 0;
          for ( uint32_t t = 0; t < 2; ++t )
          {
            onset |= UINT64_C( 1 ) << iset4_combinations[i][k][t];
          }
          best_iset_onset.push_back( onset );

          uint64_t offset = 0;
          for ( uint32_t t = 0; t < 2; ++t )
          {
            offset |= UINT64_C( 1 ) << iset4_off_set[i][k][t];
          }
          best_iset_offset.push_back( offset );
          best_care_sets.push_back( ~( kitty::dynamic_truth_table( isets[0].num_vars() ) ) );
        }
      }
    }

    if ( verbose && best_bound_sets.size() )
    {
      std::cout << "Best bound sets:\n";
      kitty::print_hex( best_bound_sets[0] );
      std::cout << " with ONSET " << best_iset_onset[0];
      std::cout << ", OFFSET " << best_iset_offset[0] << "\n";
      kitty::print_hex( best_bound_sets[1] );
      std::cout << " with ONSET " << best_iset_onset[1];
      std::cout << ", OFFSET " << best_iset_offset[1] << "\n";
      std::cout << fmt::format( "Using {} LUTs and {} leaves\n", best_cost_luts, best_cost_edges );
    }
  }

  std::vector<ac_decomposition_result> generate_decomposition( uint32_t free_set_size )
  {
    std::vector<ac_decomposition_result> res;

    for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
    {
      ac_decomposition_result dec;
      auto tt = best_bound_sets[i];
      auto care = best_care_sets[i];

      /* compute and minimize support for bound set variables */
      uint32_t k = 0;
      for ( uint32_t j = 0; j < num_vars - free_set_size; ++j )
      {
        if ( !kitty::has_var( tt, j ) )
          continue;

        if ( !kitty::has_var( tt, care, j ) )
        {
          /* fix truth table */
          adjust_truth_table_on_dc( tt, care, j );
          continue;
        }

        if ( k < j )
        {
          kitty::swap_inplace( tt, k, j );
          kitty::swap_inplace( care, k, j );
        }
        dec.support.push_back( permutations[free_set_size + j] );
        ++k;
      }

      dec.tt = kitty::shrink_to( tt, dec.support.size() );
      res.push_back( dec );
    }

    /* compute the decomposition for the top-level LUT */
    compute_top_lut_decomposition( res, free_set_size );

    return res;
  }

  void compute_top_lut_decomposition( std::vector<ac_decomposition_result>& res, uint32_t free_set_size )
  {
    uint32_t top_vars = best_bound_sets.size() + free_set_size;
    assert( top_vars <= ps.lut_size );

    /* extend bound set functions with free_set_size LSB vars */
    kitty::dynamic_truth_table tt( top_vars );

    /* compute support */
    res.emplace_back();
    for ( uint32_t i = 0; i < free_set_size; ++i )
    {
      res.back().support.push_back( permutations[i] );
    }

    /* create functions for bound set */
    std::vector<kitty::dynamic_truth_table> bound_set_vars;
    for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
    {
      bound_set_vars.emplace_back( top_vars );
      kitty::create_nth_var( bound_set_vars[i], free_set_size + i );

      /* add bound-set variables to the support */
      res.back().support.push_back( num_vars + i );
    }

    /* create final function */
    for ( uint32_t i = 0; i < best_free_set_tts.size(); ++i )
    {
      kitty::dynamic_truth_table free_set_tt = kitty::shrink_to( best_free_set_tts[i], top_vars );

      /* find MUX assignments */
      for ( uint32_t j = 0; j < bound_set_vars.size(); ++j )
      {
        /* AND with ONSET or OFFSET */
        if ( ( ( best_iset_onset[j] >> i ) & 1 ) )
        {
          free_set_tt &= bound_set_vars[j];
        }
        else if ( ( ( best_iset_offset[j] >> i ) & 1 ) )
        {
          free_set_tt &= ~bound_set_vars[j];
        }
      }

      tt |= free_set_tt;
    }

    /* add top-level LUT to result */
    res.back().tt = tt;
  }

  inline void reposition_late_arriving_variables( std::vector<uint32_t> const& late_arriving )
  {
    if ( late_arriving.empty() )
      return;

    for ( uint32_t i = 0; i < late_arriving.size(); ++i )
    {
      if ( permutations[i] == late_arriving[i] )
        continue;

      uint32_t j = i + 1;
      while ( permutations[j] != late_arriving[i] )
      {
        ++j;
      }

      std::swap( permutations[i], permutations[j] );
      kitty::swap_inplace( best_tt, i, j );
    }
  }

  inline klut_network get_result_ntk_impl()
  {
    klut_network ntk;

    /* starting from index 2 */
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      ntk.create_pi();
    }

    /* starting from index 2 + num_vars */
    klut_network::signal f;
    for ( auto const& lut : dec_result )
    {
      std::vector<klut_network::signal> children;

      for ( auto index : lut.support )
      {
        children.push_back( index + 2 );
      }

      f = ntk.create_node( children, lut.tt );
    }

    ntk.create_po( f );

    return ntk;
  }

  template<class Iterator>
  void print_perm( Iterator begin, uint32_t free_set )
  {
    std::cout << "[";
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      if ( i == free_set )
      {
        std::cout << ", ";
      }
      std::cout << *begin << " ";
      ++begin;
    }
    std::cout << "]\n";
  }

  void generate_support_minimization_encodings()
  {
    uint32_t count = 0;

    /* enable don't cares only if not a power of 2 */
    uint32_t num_combs = 3;
    if ( __builtin_popcount( best_multiplicity ) == 1 )
    {
      num_combs = 1 << best_multiplicity;
      support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
      generate_support_minimization_encodings_rec<false>( 0, 0, 0, count );
    }
    else
    {
      for ( uint32_t i = 1; i < best_multiplicity; ++i )
      {
        num_combs = ( num_combs << 1 ) + num_combs;
      }
      support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
      generate_support_minimization_encodings_rec<true>( 0, 0, 0, count );
    }

    assert( count = num_combs );

    /* print combinations */
    // std::cout << "{ ";
    // for ( auto const& entry : support_minimization_encodings )
    // {
    //   std::cout << "{ " << entry[0] << ", " << entry[1] << " }, ";
    // }
    // std::cout << "}\n";
  }

  template<bool enable_dcset>
  bool generate_support_minimization_encodings_rec( uint64_t onset, uint64_t offset, uint32_t var, uint32_t& count )
  {
    if ( var == best_multiplicity )
    {
      support_minimization_encodings[count][0] = onset;
      support_minimization_encodings[count][1] = offset;
      ++count;
      return;
    }

    /* move var in DCSET */
    if constexpr ( enable_dcset )
    {
      generate_support_minimization_encodings_rec<enable_dcset>( onset, offset, var + 1, count );
    }

    /* move var in ONSET */
    onset |= 1 << var;
    generate_support_minimization_encodings_rec<enable_dcset>( onset, offset, var + 1, count );
    onset &= ~( 1 << var );

    /* move var in OFFSET */
    offset |= 1 << var;
    generate_support_minimization_encodings_rec<enable_dcset>( onset, offset, var + 1, count );
    offset &= ~( 1 << var );
  }

  void solve_min_support_exact( std::vector<STT> const& isets, uint32_t free_set_size )
  {
    std::vector<encoding_matrix> matrix;
    matrix.reserve( support_minimization_encodings.size() );
    best_bound_sets.clear();

    /* create covering matrix */
    if ( !create_covering_matrix( isets, matrix, free_set_size, best_multiplicity > 4 ) )
    {
      return;
    }

    /* solve the covering problem */
    std::array<uint32_t, 5> solution = covering_solve_exact<true>( matrix, 100 );

    /* check for failed decomposition */
    if ( solution[0] == UINT32_MAX )
    {
      return;
    }

    /* compute best bound sets */
    uint32_t num_luts = 1 + solution[4];
    uint32_t num_levels = 2;
    uint32_t num_edges = free_set_size + solution[4];
    uint32_t isets_support = num_vars - free_set_size;
    best_care_sets.clear();
    best_iset_onset.clear();
    best_iset_offset.clear();
    for ( uint32_t i = 0; i < solution[4]; ++i )
    {
      STT tt;
      STT care;

      const uint32_t onset = support_minimization_encodings[matrix[solution[i]].index][0];
      const uint32_t offset = support_minimization_encodings[matrix[solution[i]].index][1];
      for ( uint32_t j = 0; j < best_multiplicity; ++j )
      {
        if ( ( ( onset >> j ) & 1 ) )
        {
          tt |= isets[j];
        }
        if ( ( ( offset >> j ) & 1 ) )
        {
          care |= isets[j];
        }
      }

      care |= tt;
      num_edges += matrix[solution[i]].cost & ( ( 1 << isets_support ) - 1 );

      best_bound_sets.push_back( tt );
      best_care_sets.push_back( care );
      best_iset_onset.push_back( onset );
      best_iset_offset.push_back( offset );
    }

    if ( pst != nullptr )
    {
      pst->num_luts = num_luts;
      pst->num_levels = num_levels;
      pst->num_edges = num_edges;
    }
  }

  bool create_covering_matrix( std::vector<STT> const& isets, std::vector<encoding_matrix>& matrix, uint32_t free_set_size, bool sort )
  {
    assert( best_multiplicity < 12 );
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    uint64_t sol_existance = 0;
    uint32_t iset_support = num_vars - free_set_size;

    /* insert dichotomies */
    for ( uint32_t i = 0; i < support_minimization_encodings.size(); ++i )
    {
      uint32_t const onset = support_minimization_encodings[i][0];
      uint32_t const offset = support_minimization_encodings[i][1];

      uint32_t ones_onset = __builtin_popcount( onset );
      uint32_t ones_offset = __builtin_popcount( offset );

      /* filter columns that do not distinguish pairs */
      if ( ones_onset == 0 || ones_offset == 0 || ones_onset == best_multiplicity || ones_offset == best_multiplicity )
      {
        continue;
      }

      /* compute function and distinguishable seed dichotomies */
      uint64_t column = 0;
      STT tt;
      STT care;
      uint32_t pair_pointer = 0;
      for ( uint32_t j = 0; j < best_multiplicity; ++j )
      {
        if ( ( ( onset >> j ) & 1 ) )
        {
          tt |= isets[j];
        }

        if ( ( ( offset >> j ) & 1 ) )
        {
          care |= isets[j];
        }

        /* compute included seed dichotomies */
        for ( uint32_t k = j + 1; k < best_multiplicity; ++k )
        {
          /* if is not in ONSET and is in OFFSET */
          uint32_t test_pair = ( onset >> j ) & ( ( ~onset & offset ) >> k );

          if ( ( test_pair & 1 ) )
          {
            column |= UINT64_C( 1 ) << ( pair_pointer );
          }

          ++pair_pointer;
        }
      }

      care |= tt;

      /* compute cost */
      uint32_t cost = 0;
      for ( uint32_t j = 0; j < iset_support; ++j )
      {
        cost += kitty::has_var( tt, care, j ) ? 1 : 0;
      }

      /* discard solutions with support over LUT size */
      if ( cost > ps.lut_size )
        continue;

      if ( cost > 1 )
      {
        cost |= 1 << iset_support;
      }

      uint32_t sort_cost = cost + ( ( combinations - __builtin_popcountl( column ) ) << num_vars );

      /* insert */
      matrix.emplace_back( encoding_matrix{ column, cost, i, sort_cost } );

      sol_existance |= column;
    }

    /* necessary condition for the existance of a solution */
    if ( __builtin_popcountl( sol_existance ) != combinations )
    {
      return false;
    }

    if ( !sort )
    {
      return true;
    }

    std::sort( matrix.begin(), matrix.end(), [&]( auto const& a, auto const& b ) {
      return a.sort_cost < b.sort_cost;
    } );

    /* print */
    // if ( best_multiplicity < 6 )
    // {
    //   for ( uint32_t i = 0; i < columns.size(); ++i )
    //   {
    //     std::cout << indexes[i] << " " << costs[i] << " \t" << columns[i] << "\n";
    //   }
    // }

    return true;
  }

  template<bool limit_iter = false>
  std::array<uint32_t, 5> covering_solve_exact( std::vector<encoding_matrix>& matrix, uint32_t max_iter = 100 )
  {
    /* last value of res contains the size of the bound set */
    std::array<uint32_t, 5> res = { UINT32_MAX };
    uint32_t best_cost = UINT32_MAX;
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    bool looping = true;

    assert( best_multiplicity <= 16 );

    /* determine the number of needed loops*/
    if ( best_multiplicity <= 4 )
    {
      res[4] = 2;
      for ( uint32_t i = 0; i < matrix.size() - 1; ++i )
      {
        for ( uint32_t j = 1; j < matrix.size(); ++j )
        {
          /* filter by cost */
          if ( matrix[i].cost + matrix[j].cost >= best_cost )
            continue;

          /* check validity */
          if ( __builtin_popcountl( matrix[i].column | matrix[j].column ) == combinations )
          {
            res[0] = i;
            res[1] = j;
            best_cost = matrix[i].cost + matrix[j].cost;
          }
        }
      }
    }
    else if ( best_multiplicity <= 8 )
    {
      res[4] = 3;
      for ( uint32_t i = 0; i < matrix.size() - 2 && looping; ++i )
      {
        /* limit */
        if constexpr ( limit_iter )
        {
          if ( best_cost < UINT32_MAX && --max_iter == 0 )
          {
            looping = false;
          }
        }

        for ( uint32_t j = 1; j < matrix.size() - 1 && looping; ++j )
        {
          uint64_t current_columns = matrix[i].column | matrix[j].column;
          uint32_t current_cost = matrix[i].cost + matrix[j].cost;

          /* limit */
          if constexpr ( limit_iter )
          {
            if ( best_cost < UINT32_MAX && --max_iter == 0 )
            {
              looping = false;
            }
          }

          /* bound */
          if ( current_cost >= best_cost )
          {
            continue;
          }

          for ( uint32_t k = 2; k < matrix.size() && looping; ++k )
          {
            /* limit */
            if constexpr ( limit_iter )
            {
              if ( best_cost < UINT32_MAX && --max_iter == 0 )
              {
                looping = false;
              }
            }

            /* filter by cost */
            if ( current_cost + matrix[k].cost >= best_cost )
              continue;

            /* check validity */
            if ( __builtin_popcountl( current_columns | matrix[k].column ) == combinations )
            {
              res[0] = i;
              res[1] = j;
              res[2] = k;
              best_cost = current_cost + matrix[k].cost;
            }
          }
        }
      }
    }
    else
    {
      res[4] = 4;
      for ( uint32_t i = 0; i < matrix.size() - 3 && looping; ++i )
      {
        /* limit */
        if constexpr ( limit_iter )
        {
          if ( best_cost < UINT32_MAX && --max_iter == 0 )
          {
            looping = false;
          }
        }

        for ( uint32_t j = 1; j < matrix.size() - 2 && looping; ++j )
        {
          uint64_t current_columns0 = matrix[i].column | matrix[j].column;
          uint32_t current_cost0 = matrix[i].cost + matrix[j].cost;

          /* limit */
          if constexpr ( limit_iter )
          {
            if ( best_cost < UINT32_MAX && --max_iter == 0 )
            {
              looping = false;
            }
          }

          /* bound */
          if ( current_cost0 >= best_cost )
          {
            continue;
          }

          for ( uint32_t k = 2; k < matrix.size() - 1 && looping; ++k )
          {
            uint64_t current_columns1 = current_columns0 | matrix[k].column;
            uint32_t current_cost1 = current_cost0 + matrix[k].cost;

            /* limit */
            if constexpr ( limit_iter )
            {
              if ( best_cost < UINT32_MAX && --max_iter == 0 )
              {
                looping = false;
              }
            }

            /* bound */
            if ( current_cost1 >= best_cost )
            {
              continue;
            }

            for ( uint32_t t = 3; t < matrix.size() && looping; ++t )
            {
              /* limit */
              if constexpr ( limit_iter )
              {
                if ( best_cost < UINT32_MAX && --max_iter == 0 )
                {
                  looping = false;
                }
              }

              /* filter by cost */
              if ( current_cost1 + matrix[t].cost >= best_cost )
                continue;

              /* check validity */
              if ( __builtin_popcountl( current_columns1 | matrix[t].column ) == combinations )
              {
                res[0] = i;
                res[1] = j;
                res[2] = k;
                res[3] = t;
                best_cost = current_cost1 + matrix[t].cost;
              }
            }
          }
        }
      }
    }

    return res;
  }

  void adjust_truth_table_on_dc( STT& tt, STT& care, uint32_t var_index )
  {
    assert( var_index < tt.num_vars() );
    assert( tt.num_vars() == care.num_vars() );

    if ( tt.num_vars() <= 6 || var_index < 6 )
    {
      auto it_tt = std::begin( tt._bits );
      auto it_care = std::begin( care._bits );
      while ( it_tt != std::end( tt._bits ) )
      {
        uint64_t new_bits = *it_tt & *it_care;
        *it_tt = ( ( new_bits | ( new_bits >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index] ) |
                 ( ( new_bits | ( new_bits << ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections[var_index] );
        *it_care = *it_care | ( *it_care >> ( uint64_t( 1 ) << var_index ) );

        ++it_tt;
        ++it_care;
      }
      return;
    }

    const auto step = 1 << ( var_index - 6 );
    for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); i += 2 * step )
    {
      for ( auto j = 0; j < step; ++j )
      {
        tt._bits[i + j] = ( tt._bits[i + j] & care._bits[i + j] ) | ( tt._bits[i + j + step] & care._bits[i + j + step] );
        tt._bits[i + j + step] = tt._bits[i + j];
        care._bits[i + j] = care._bits[i + j] | care._bits[i + j + step];
        care._bits[i + j + step] = care._bits[i + j];
      }
    }
  }

  void local_extend_to( STT& tt, uint32_t real_num_vars )
  {
    if ( real_num_vars < 6 )
    {
      auto mask = *tt.begin();

      for ( auto i = real_num_vars; i < num_vars; ++i )
      {
        mask |= ( mask << ( 1 << i ) );
      }

      std::fill( tt.begin(), tt.end(), mask );
    }
    else
    {
      uint32_t num_blocks = ( 1u << ( real_num_vars - 6 ) );
      auto it = tt.begin();
      while ( it != tt.end() )
      {
        it = std::copy( tt.cbegin(), tt.cbegin() + num_blocks, it );
      }
    }
  }

  bool verify_equivalence_impl()
  {
    klut_network klut = get_result_ntk_impl();

    default_simulator<kitty::dynamic_truth_table> sim( num_vars );
    const auto tt_res = simulate<kitty::dynamic_truth_table>( klut, sim )[0];
    kitty::dynamic_truth_table cec = kitty::shrink_to( tt_start, tt_start.num_vars() );

    return tt_res == cec;
  }

private:
  uint32_t best_multiplicity{ UINT32_MAX };
  STT best_tt;
  std::vector<STT> best_bound_sets;
  std::vector<STT> best_care_sets;
  std::vector<STT> best_free_set_tts;
  std::vector<uint64_t> best_iset_onset;
  std::vector<uint64_t> best_iset_offset;
  std::vector<ac_decomposition_result> dec_result;

  std::vector<std::array<uint32_t, 2>> support_minimization_encodings;

  TT tt_start;
  uint32_t num_vars;
  ac_decomposition_params const& ps;
  ac_decomposition_stats* pst;
  std::vector<uint32_t> permutations;

  // static constexpr uint32_t iset3_combinations[][2] = { { 1, 6 }, { 2, 5 }, { 4, 3 } };

  static constexpr uint32_t iset3_combinations[][2][2] =
      { { { 0 }, { 1 } },
        { { 1 }, { 0 } },
        { { 2 }, { 0 } } };

  static constexpr uint32_t iset3_off_set[][2][2] =
      { { { 1, 2 }, { 2 } },
        { { 0, 2 }, { 2 } },
        { { 0, 1 }, { 1 } } };

  static constexpr uint32_t iset4_combinations[][2][2] =
      { { { 1, 3 }, { 2, 3 } },
        { { 1, 2 }, { 2, 3 } },
        { { 0, 2 }, { 2, 3 } },
        { { 0, 3 }, { 0, 1 } },
        { { 0, 3 }, { 0, 2 } },
        { { 0, 3 }, { 2, 3 } } };

  static constexpr uint32_t iset4_off_set[][2][2] =
      { { { 0, 2 }, { 0, 1 } },
        { { 0, 3 }, { 0, 1 } },
        { { 1, 3 }, { 0, 1 } },
        { { 1, 2 }, { 2, 3 } },
        { { 1, 2 }, { 1, 3 } },
        { { 1, 2 }, { 0, 1 } } };
};

} // namespace detail

} // namespace mockturtle