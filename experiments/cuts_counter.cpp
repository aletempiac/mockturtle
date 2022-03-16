/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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
#include <unordered_set>
#include <iostream>
#include <fstream>

#include <parallel_hashmap/phmap.h>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <kitty/npn.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/hash.hpp>
#include <kitty/print.hpp>
#include <kitty/operations.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/cut_enumeration/exact_map_cut.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <experiments.hpp>

template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 6 )
{
  mockturtle::write_verilog( ntk, "/tmp/network.v" );
  system( fmt::format( "abc -q \"/tmp/network.v; &get; &if -a -K {}; &put; write_blif /tmp/output.blif\"", k ).c_str() );
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR LUT" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}

void abc_map_functions( uint32_t k )
{
  using namespace experiments;
  using namespace mockturtle;
  using tt_hash = kitty::hash<kitty::static_truth_table<6>>;
  using lib_t = phmap::flat_hash_map<kitty::static_truth_table<6>, uint32_t, tt_hash>;

  lib_t functions;

  auto i = 0;
  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( ++i > 10 )
      break;

    fmt::print( "[i] processing {}\n", benchmark );
    mig_network ntk;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( ntk ) ) != lorina::return_code::success )
    {
      continue;
    }

    auto klut = lut_map( ntk, k );

    /* foreach cut obtain the NPN representative */
    klut.foreach_gate( [&]( auto const& n ) {
      if ( klut.fanin_size( n ) > 4 )
      {
        auto tt_d = klut.node_function( n );
        const auto support = kitty::min_base_inplace( tt_d );
        if ( support.size() > 4 )
        {
          const auto tt = kitty::extend_to<6>( tt_d );
          const auto config = exact_npn_canonization_minimized( tt, support.size() );
          const auto& npn_tt = std::get<0>( config );
          functions[npn_tt]++;
        }
      }
    } );
  }

  std::ofstream out( "abc_functions_" + std::to_string( k ) + ".txt" );

  std::vector<std::pair<kitty::static_truth_table<6>, uint32_t>> functions_sorted;
  functions_sorted.reserve( functions.size() );

  for ( auto const& p : functions )
  {
    functions_sorted.push_back( p );
  }

  std::sort( functions_sorted.begin(), functions_sorted.end(), []( auto const& a, auto const& b ) {
    return a.second > b.second;
  } );

  for ( auto const& p : functions_sorted )
  {
    out << kitty::to_hex( p.first ) << " " << p.second << "\n";
  }

  out.close();
}

void generate_functions()
{
  using namespace experiments;
  using namespace mockturtle;
  using tt_hash = kitty::hash<kitty::static_truth_table<6>>;
  using lib_t = phmap::flat_hash_map<kitty::static_truth_table<6>, uint32_t, tt_hash>;

  lib_t functions;

  auto i = 0;
  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( ++i > 10 )
      break;

    fmt::print( "[i] processing {}\n", benchmark );
    mig_network ntk;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( ntk ) ) != lorina::return_code::success )
    {
      continue;
    }

    cut_enumeration_params ps;
    ps.cut_limit = 24;
    ps.minimize_truth_table = true;

    auto cuts = fast_cut_enumeration<mig_network, 6, true, cut_enumeration_exact_map_cut>( ntk, ps );

    /* foreach cut obtain the NPN representative */
    ntk.foreach_gate( [&]( auto const& n ) {
      const auto index = ntk.node_to_index( n );
      
      for ( auto const& cut : cuts.cuts( index ) )
      {
        if ( cut->size() > 4 )
        {
          const auto tt = cuts.truth_table( *cut );
          const auto config = exact_npn_canonization_minimized( tt, cut->size() );
          const auto& npn_tt = std::get<0>( config );
          functions[npn_tt]++;
        }
      }
    } );
  }

  std::ofstream out( "functions_study_24.txt" );

  std::vector<std::pair<kitty::static_truth_table<6>, uint32_t>> functions_sorted;
  functions_sorted.reserve( functions.size() );

  for ( auto const& p : functions )
  {
    functions_sorted.push_back( p );
  }

  std::sort( functions_sorted.begin(), functions_sorted.end(), []( auto const& a, auto const& b ) {
    return a.second > b.second;
  } );

  for ( auto const& p : functions_sorted )
  {
    out << kitty::to_hex( p.first ) << " " << p.second << "\n";
  }

  out.close();
}

void read_functions( std::ifstream& in, std::unordered_set<kitty::static_truth_table<6>, kitty::hash<kitty::static_truth_table<6>>>& functions, uint32_t limit )
{
  uint32_t i = 0;
  std::string tt_s;
  uint32_t instances;

  while( ( in >> tt_s >> instances ) && i++ < limit )
  {
    kitty::static_truth_table<6> tt;
    kitty::create_from_hex_string( tt, tt_s );

    if ( functions.find( tt ) == functions.end() )
    {
      functions.insert( tt );
    }
  }
}

void merge_databases()
{
  using namespace experiments;
  using namespace mockturtle;
  using tt_hash = kitty::hash<kitty::static_truth_table<6>>;

  std::ifstream in_study( "functions_study_24.txt" );
  std::ifstream in_abc_5( "abc_functions_5.txt" );
  std::ifstream in_abc_6( "abc_functions_6.txt" );

  std::unordered_set<kitty::static_truth_table<6>, tt_hash> functions;
  functions.reserve( 500 );

  read_functions( in_study, functions, 180 );
  read_functions( in_abc_5, functions, 67 );
  read_functions( in_abc_6, functions, 97 );

  in_study.close();
  in_abc_5.close();
  in_abc_6.close();

  std::ofstream out( "functions_merge.txt" );

  for ( auto const& tt : functions )
  {
    out << kitty::to_hex( tt ) << "\n";
  }

  out.close();
}

int main()
{
  generate_functions();
  abc_map_functions( 5 );
  abc_map_functions( 6 );
  merge_databases();
  return 0;
}
