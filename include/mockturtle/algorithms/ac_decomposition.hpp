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
#include <optional>
#include <cassert>
#include <cstdint>
#include <type_traits>
#include <vector>
#include <unordered_map>

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>
#include <kitty/print.hpp>
#include <kitty/traits.hpp>

#include "../networks/klut.hpp"

namespace mockturtle
{

/*! \brief Parameters for ac_decomposition */
struct ac_decomposition_params
{
  /*! \brief LUT size for decomposition. */
  uint32_t lut_size{ 6 };
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
public:
  ac_decomposition_impl( TT const& tt, uint32_t num_vars, ac_decomposition_params const& ps )
      : tt_start( tt )
      , num_vars( num_vars )
      , ps( ps )
      , permutations( num_vars )
  {
    std::iota( permutations.begin(), permutations.end(), 0 );
  }

  /*! \brief Runs ACD using late arriving variables */
  uint32_t run( std::vector<uint32_t> late_arriving )
  {
    best_tt = tt_start;
    best_multiplicity = UINT32_MAX;

    /* return a high cost if too many late arriving variables */
    if ( late_arriving.size() > ps.lut_size / 2 || late_arriving.size() > 3 )
    {
      return UINT32_MAX;
    }

    /* permute late arriving variables to be the least significant */
    reposition_late_arriving_variables( late_arriving );

    /* run ACD trying different bound sets and free sets */
    uint32_t free_set_size = late_arriving.size();
    uint32_t offset = late_arriving.size();
    for ( uint32_t i = offset; i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, i ); };
      auto [tt_p, perm, cost] = enumerate_iset_combinations_offset( i, offset, evaluate_fn, true );

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
    std::vector<kitty::dynamic_truth_table> isets = compute_isets( free_set_size );

    /* test for column multiplicity 4*/
    std::vector<kitty::dynamic_truth_table> bound_sets;
    if ( best_multiplicity == 4 )
    {
      test_support_minimization_isets( isets, true );
    }

    /* if feasible decomposition */
    if ( best_bound_sets.empty() )
    {
      /* TODO: change in return empty */
      return best_multiplicity;
    }

    auto decomposition = generate_decomposition( free_set_size );

    dec_result = decomposition;

    /* TODO: change return value */
    return best_multiplicity;
  }

  /*! \brief Runs ACD using late arriving variables and guaranteeing support minimization */
  uint32_t run_dsd( std::vector<uint32_t> late_arriving )
  {
    best_tt = tt_start;
    best_multiplicity = UINT32_MAX;

    /* compute minimum number of variables in the free set */
    uint32_t dsd_vars = num_vars - ps.lut_size;

    /* return a high cost if too many late arriving variables */
    if ( late_arriving.size() > ps.lut_size / 2 || late_arriving.size() > 3 || dsd_vars > 3 )
    {
      return UINT32_MAX;
    }

    /* permute late arriving variables to be the least significant */
    reposition_late_arriving_variables( late_arriving );

    /* run ACD trying different bound sets and free sets */
    uint32_t free_set_size = late_arriving.size();
    uint32_t offset = late_arriving.size();
    for ( uint32_t i = std::max( dsd_vars, offset ); i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, i ); };
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
    assert( ps.lut_size > 3 && ps.lut_size < 13 );
    best_tt = tt_start;
    best_multiplicity = UINT32_MAX;

    uint32_t free_set_size = 1;
    for ( uint32_t i = 1; i <= ps.lut_size / 2 && i <= 3; ++i )
    {
      auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, i ); };
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

    return best_multiplicity;
  }

  /*! \brief Runs ACD trying different bound sets */
  uint32_t run_offset( uint32_t free_set_size, uint32_t offset )
  {
    best_tt = tt_start;
    auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, free_set_size ); };

    auto [tt_p, perm, cost] = enumerate_iset_combinations_offset( free_set_size, offset, evaluate_fn, false );
    best_tt = tt_p;
    permutations = perm;
    best_multiplicity = cost;

    return best_multiplicity;
  }

  /*! \brief Runs ACD trying different bound sets */
  uint32_t run( uint32_t free_set_size )
  {
    best_tt = tt_start;
    auto evaluate_fn = [&] ( TT const& tt ) { return column_multiplicity( tt, free_set_size ); };

    auto [tt_p, perm, cost] = enumerate_iset_combinations( free_set_size, evaluate_fn, false );
    best_tt = tt_p;
    permutations = perm;
    best_multiplicity = cost;

    return best_multiplicity;
  }

  /*! \brief Runs ACD without trying different bound sets */
  uint32_t run_no_permutations( uint32_t free_set_size )
  {
    best_tt = tt_start;
    best_multiplicity = column_multiplicity( tt_start, free_set_size );

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
  
private:
  uint32_t column_multiplicity( TT tt, uint32_t free_set_size )
  {
    uint64_t multiplicity_set[4] = { 0u, 0u, 0u, 0u };
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
          multiplicity_set[0] |= UINT64_C( 1 ) << ( *it & 0x3 );
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
          multiplicity_set[0] |= UINT64_C( 1 ) << ( *it & 0xF );
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
  std::tuple<TT, std::vector<uint32_t>, uint32_t> enumerate_iset_combinations( uint32_t free_set_size, Fn&& fn, bool verbose = false )
  {
    TT tt = best_tt;

    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* special case */
    if ( num_vars <= free_set_size || free_set_size == 0 )
    {
      return { tt, permutations, UINT32_MAX };
    }

    /* select k */
    free_set_size = std::min( free_set_size, num_vars - free_set_size );

    /* init permutation array */
    std::array<uint32_t, 16> perm, best_perm;
    std::copy( permutations.begin(), permutations.begin() + num_vars, perm.begin() );
    best_perm = perm;

    /* TT with best cost */
    TT best_tt = tt;
    uint32_t best_cost = UINT32_MAX;

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
        print_perm( perm.begin(), free_set_size );
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
          print_perm( perm.begin(), free_set_size );
        }
      }
    }
    else if ( free_set_size == 2 )
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
          print_perm( perm.begin(), free_set_size );
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
            best_tt = tt;
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
              best_tt = tt;
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

    return std::make_tuple( best_tt, res_perm, best_cost );
  }

  template<typename Fn>
  std::tuple<TT, std::vector<uint32_t>, uint32_t> enumerate_iset_combinations_offset( uint32_t free_set_size, uint32_t offset, Fn&& fn, bool verbose = false )
  {
    TT tt = best_tt;

    /* TT with best cost */
    TT best_tt = tt;
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

  std::vector<kitty::dynamic_truth_table> compute_isets( uint32_t free_set_size, bool verbose = false )
  {
    /* construct isets involved in multiplicity */
    std::vector<kitty::dynamic_truth_table> isets;
    isets.reserve( best_multiplicity );
    
    for ( uint32_t i = 0; i < best_multiplicity; ++i )
    {
      isets.push_back( kitty::create<kitty::dynamic_truth_table>( num_vars - free_set_size ) );
    }

    /* construct isets */
    std::unordered_map<uint64_t, uint32_t> column_to_iset;
    TT tt = best_tt;
    uint32_t offset = 0, block = 0;

    if ( free_set_size == 1 )
    {
      auto it = std::begin( tt );
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
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
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
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
      for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); ++i )
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

    /* save free_set functions */
    std::vector<kitty::dynamic_truth_table> free_set_tts;
    free_set_tts.reserve( best_multiplicity );

    for ( uint32_t i = 0; i < best_multiplicity; ++i )
    {
      free_set_tts.emplace_back( free_set_size );
    }
    for ( auto const& pair : column_to_iset )
    {
      free_set_tts[pair.second]._bits[0] = pair.first;
    }

    best_free_set_tts = free_set_tts;

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

      /* compute and minimize support for bound set variables */
      uint32_t k = 0;
      for ( uint32_t j = 0; j < num_vars - free_set_size; ++j )
      {
        if ( !kitty::has_var( best_bound_sets[i], j ) )
          continue;

        if ( k < j )
        {
          kitty::swap_inplace( tt, k, j );
        }
        dec.support.push_back( permutations[free_set_size + j] );
        ++k;
      }

      if ( dec.support.size() < tt.num_vars() )
      {
        dec.tt = kitty::shrink_to( tt, dec.support.size() );
      }
      else
      {
        dec.tt = tt;
      }

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
      auto free_set_tt = kitty::extend_to( best_free_set_tts[i], top_vars );

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

private:
  uint32_t best_multiplicity{ UINT32_MAX };
  TT best_tt;
  std::vector<kitty::dynamic_truth_table> best_bound_sets;
  std::vector<kitty::dynamic_truth_table> best_free_set_tts;
  std::vector<uint64_t> best_iset_onset;
  std::vector<uint64_t> best_iset_offset;
  std::vector<ac_decomposition_result> dec_result;

  TT const& tt_start;
  uint32_t num_vars;
  ac_decomposition_params const& ps;
  std::vector<uint32_t> permutations;

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