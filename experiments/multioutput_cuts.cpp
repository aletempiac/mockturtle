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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <kitty/kitty.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <parallel_hashmap/phmap.h>

#include <experiments.hpp>

#define CUT_K_MAX 3
#define CUT_L_MAX 2

using namespace mockturtle;

struct cut_hash
{
  uint64_t operator()( const std::array<uint32_t, CUT_K_MAX>& p ) const
  {
    uint64_t seed = hash_block( p[0] );

    for ( auto i = 1; i < CUT_K_MAX; ++i )
      hash_combine( seed, hash_block( p[i] ) );

    return seed;
  }
};

struct vtt_hash
{
  uint64_t operator()( const std::array<kitty::static_truth_table<CUT_K_MAX>, CUT_L_MAX>& vtt ) const
  {
    uint64_t seed = hash_block( vtt[0]._bits );

    for ( auto i = 1; i < CUT_L_MAX; ++i )
      hash_combine( seed, hash_block( vtt[i]._bits ) );

    return seed;
  }
};

using cuts_counter_t = phmap::flat_hash_map<std::array<kitty::static_truth_table<CUT_K_MAX>, CUT_L_MAX>, uint32_t, vtt_hash>;
using leaves_hash_t = phmap::flat_hash_map<std::array<uint32_t, CUT_K_MAX>, std::vector<uint64_t>, cut_hash>;
using multi_cuts_t = std::vector<std::vector<uint64_t>>;

template<class Ntk>
bool check_tfi_valid_rec( Ntk const& ntk, node<Ntk> const& n, node<Ntk> const& root, node<Ntk> const& target, bool& valid )
{
  /* leaf */
  if ( ntk.value( n ) )
    return false;

  /* already visited */
  if ( ntk.visited( n ) == ntk.trav_id() )
    return false;

  ntk.set_visited( n, ntk.trav_id() );

  if ( n == target )
  {
    valid = ntk.fanout_size( n ) != 1;
    return true;
  }

  bool found = false;
  ntk.foreach_fanin( n, [&]( auto const& f ) {
    found |= check_tfi_valid_rec( ntk, ntk.get_node( f ), root, target, valid );
    return valid;
  } );

  if ( found && n != root && ntk.fanout_size( n ) > 1 )
    valid = false;

  return found;
}

template<class Ntk, typename Cut_t>
bool check_compatibility( Ntk const& ntk, uint32_t index1, uint32_t index2, Cut_t const& cut )
{
  bool valid = true;

  /* check containment of cut1 in cut2 and viceversa */
  if ( index1 > index2 )
  {
    std::swap( index1, index2 );
  }

  ntk.foreach_fanin( ntk.index_to_node( index2 ), [&]( auto const& f ) {
    auto g = ntk.get_node( f );
    if ( ntk.node_to_index( g ) == index1 && ntk.fanout_size( g ) == 1 )
    {
      valid = false;
    }
    return valid;
  } );

  if ( !valid )
    return false;
  
  /* check containment when node is reachable from middle nodes */

  /* reference cut leaves */
  for ( auto leaf : cut )
  {
    ntk.incr_value( ntk.index_to_node( leaf ) );
  }

  ntk.incr_trav_id();
  check_tfi_valid_rec( ntk, ntk.index_to_node( index2 ), ntk.index_to_node( index2 ), ntk.index_to_node( index1 ), valid );

  /* dereference leaves */
  for ( auto leaf : cut )
  {
    ntk.decr_value( ntk.index_to_node( leaf ) );
  }

  return valid;
}

template<class Ntk, typename NetCuts>
void create_classes( Ntk const& ntk, NetCuts const& cuts, leaves_hash_t& cuts_classes, uint32_t k )
{
  std::array<uint32_t, CUT_K_MAX> leaves;

  ntk.foreach_gate( [&]( auto const& n ) {
    uint32_t cut_index = 0;
    for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
    {
      kitty::dynamic_truth_table tt = cuts.truth_table( *cut );

      if ( cut->size() != k )
      {
        ++cut_index;
        continue;
      }
      
      uint64_t data = ( static_cast<uint64_t>( ntk.node_to_index( n ) ) << 16 ) | cut_index;
      leaves.fill( 0 );

      uint32_t i = 0;
      for ( auto l : *cut )
        leaves[i++] = l;

      /* add to hash table */
      auto &v = cuts_classes[leaves];
      v.push_back( data );

      ++cut_index;
    }
  } );
}

template<class Ntk, typename NetCuts>
void combine_cuts( Ntk const& ntk, NetCuts const& cuts, leaves_hash_t const& cuts_classes, multi_cuts_t &multi_cuts, uint32_t k, uint32_t l )
{
  ( void )k;
  ( void )l;

  ntk.clear_values();

  for ( auto& it : cuts_classes )
  {
    /* not matched */
    if ( it.second.size() < 2 )
      continue;
    
    for ( uint32_t i = 0; i < it.second.size() - 1; ++i )
    {
      uint64_t data_i = it.second[i];
      uint32_t index_i = data_i >> 16;
      uint32_t cut_index_i = data_i & UINT16_MAX;
      auto const& cut_i = cuts.cuts( index_i )[cut_index_i];

      for ( uint32_t j = i + 1; j < it.second.size(); ++j )
      {
        uint64_t data_j = it.second[j];
        uint32_t index_j = data_j >> 16;
        
        /* check compatibility: no internal multi-fanout */
        assert( cut_i.size() == cut_j.size() );
        if ( !check_compatibility( ntk, index_i, index_j, cut_i ) )
          continue;
        
        /* add to solutions */
        multi_cuts.emplace_back();
        multi_cuts.back().push_back( data_i );
        multi_cuts.back().push_back( data_j );
      }
    }
  }

  std::cout << "[i] Valid " << multi_cuts.size() << "\n";
}

template<typename NetCuts>
void process_and_add_cuts( NetCuts const& cuts, cuts_counter_t &cuts_counter, multi_cuts_t const& multi_cuts )
{
  for ( auto it : multi_cuts )
  {
    std::array<kitty::static_truth_table<CUT_K_MAX>, CUT_L_MAX> vtt = {};

    for ( auto i = 0u; i < it.size(); ++i )
    {
      uint64_t data = it[i];
      uint32_t index = data >> 16;
      uint32_t cut_index = data & UINT16_MAX;
      auto const& cut = cuts.cuts( index )[cut_index];

      auto tt = cuts.truth_table( cut );
      auto tt_canon = std::get<0>( kitty::exact_npn_canonization( tt ) );
      
      vtt[i] = kitty::extend_to<CUT_K_MAX>( tt_canon );
    }

    std::sort( vtt.begin(), vtt.end(), [&]( auto const& a, auto const& b ) { return a._bits > b._bits; } );

    cuts_counter[vtt]++;
  }
}

int main( int argc, char **argv)
{
  using namespace experiments;

  if ( argc < 3 )
  {
    std::cout << "[e] two arguments required: K, L\n";
    return -1;
  }

  /* set params */
  const uint32_t k = static_cast<uint32_t>( atoi( argv[1] ) );
  const uint32_t l = static_cast<uint32_t>( atoi( argv[2] ) );

  if ( k > CUT_K_MAX )
  {
    std::cout << "[e] K is maximum " << CUT_K_MAX << " at compilation time\n";
    return -1;
  }  
  if ( l != CUT_L_MAX )
  {
    std::cout << "[e] L is different from " << CUT_L_MAX << " at compilation time\n";
    return -1;
  }

  /* map to count the multi-output cuts */
  cuts_counter_t cuts_counter;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( "optimized/" + benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    cut_enumeration_params ps;
    ps.cut_size = k;
    ps.minimize_truth_table = true;
    auto cuts = cut_enumeration<aig_network, true>( aig, ps );

    /* compute multi-output cuts */
    leaves_hash_t cuts_classes;
    multi_cuts_t multi_cuts;

    create_classes( aig, cuts, cuts_classes, k );
    combine_cuts( aig, cuts, cuts_classes, multi_cuts, k, l );

    /* add cuts after canonization */
    process_and_add_cuts( cuts, cuts_counter, multi_cuts );
  }

  /* copy and sort */
  std::vector<std::pair<std::array<kitty::static_truth_table<CUT_K_MAX>, CUT_L_MAX>, uint32_t>> cuts_instances;

  for ( auto it : cuts_counter )
    cuts_instances.push_back( it );
  
  std::sort( cuts_instances.begin(), cuts_instances.end(), [&]( auto const& a, auto const& b ) { return a.second > b.second; } );

  std::cout << "[i] Detected " << cuts_instances.size() << " unique multi-output gates\n";

  /* report top 10 */
  std::cout << "[i] Report of the detected 10-most occurrent multi-output functions\n";
  for ( auto i = 0u; i < cuts_instances.size() && i < 10; ++i )
  {
    auto const& data = cuts_instances[i];
    std::cout << data.second << "\t : ";
    for ( auto j = 0u; j < CUT_L_MAX; ++j )
    {
      std::cout << "(";
      // kitty::print_hex( data.first[j] );
      kitty::print_expression( data.first[j] );
      std::cout << ")\t ";
    }
    std::cout << "\n";
  }

  return 0;
}
