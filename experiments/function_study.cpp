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
#include <iostream>
#include <fstream>
#include <unordered_map>

#include <lorina/aiger.hpp>
#include <lorina/blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>

#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;

static const uint32_t K = 5u;
static const bool opt = true;

template<class Ntk>
klut_network abc_lut_map( Ntk const& ntk, uint32_t k = 4u )
{
  write_aiger( ntk, "/tmp/dtm.aig" );
  std::string command = fmt::format( "abc -q \"r /tmp/dtm.aig; &get; &if -a -K {}; &put; write_blif /tmp/res.blif\"", k );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  klut_network res;
  std::string string_path = ( "/tmp/res.blif" );
  if ( lorina::read_blif( string_path, blif_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_blif failed" << std::endl;

  return res;
}

template<class Ntk>
aig_network abc_opt( Ntk const& ntk )
{
  write_aiger( ntk, "/tmp/dto.aig" );
  std::string command = "abc -q \"r /tmp/dto.aig; ifraig; resyn2; resyn2rs; write_aiger /tmp/res.aig\"";

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  aig_network res;
  std::string string_path = ( "/tmp/res.aig" );
  if ( lorina::read_aiger( string_path, aiger_reader( res ) ) != lorina::return_code::success )
    std::cerr << "read_aiger failed" << std::endl;

  return res;
}

std::string generate_sop( kitty::static_truth_table<K> const& stt )
{
  /* create a dynamic truth table with minimized support */
  auto stt_min = stt;
  auto support = kitty::min_base_inplace( stt_min );

  kitty::dynamic_truth_table tt = kitty::shrink_to( stt_min, support.size() );

  /* translate into ISOP */
  auto cubes_p = isop( tt );
  auto cubes_n = isop( ~tt );

  /* minimize the #cubes up to output negation */
  if ( cubes_n.size() < cubes_p.size() )
  {
    cubes_p = cubes_n;
  }

  auto num_vars = support.size();

  /* minimize the negations NPN */
  std::vector<uint32_t> lit_p( num_vars );
  std::vector<uint32_t> lit_n( num_vars );
  std::vector<bool> flip( num_vars );
  for ( auto i = 0u; i < cubes_p.size(); ++i )
  {
    auto& cube = cubes_p[i];
    for ( auto j = 0u; j < num_vars; ++j )
    {
      if ( cube.get_mask( num_vars - 1 - j ) )
      {
        if ( cube.get_bit( num_vars - 1 - j ) )
          ++lit_p[j];
        else
          ++lit_n[j];
      }
    }
  }

  for ( auto j = 0u; j < num_vars; ++j )
  {
    if ( lit_n[j] > lit_p[j] )
      flip[j] = true;
  }

  /* generate SOP as a string */
  std::stringstream sop;

  for ( auto i = 0u; i < cubes_p.size(); ++i )
  {
    auto const& cube = cubes_p[i];
    for ( auto j = 0u; j < num_vars; ++j )
    {
      if ( cube.get_mask( num_vars - 1 - j ) )
      {
        char lit = 'a' + j;
        sop << lit;
        if ( !cube.get_bit( num_vars - 1 - j ) ^ flip[j] )
          sop << "\'";
      }
    }

    if ( i < cubes_p.size() - 1 )
      sop << "+";
  }

  return sop.str();
}

int main()
{
  using namespace kitty;

  std::unordered_map<static_truth_table<K>, uint32_t, kitty::hash<static_truth_table<K>>> functions;

  for ( auto const& benchmark : all_benchmarks() )
  {
    if ( benchmark == "leon2" || benchmark == "leon3" || benchmark == "leon3_opt" || benchmark == "leon3mp" || benchmark == "netcard" )
      continue;

    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if ( opt )
      aig = abc_opt( aig );

    /* map into k-LUTs */
    auto klut = abc_lut_map( aig, K );

    /* save functions (up to permutation) */
    klut.foreach_gate( [&]( auto const& n ) {
      dynamic_truth_table const tt = klut.node_function( n );
      static_truth_table<K> tt_s = kitty::extend_to<K>( tt );
      static_truth_table<K> const tt_p = std::get<0>( exact_npn_canonization( tt_s ) ); /* try with NPN as well */
      
      functions[tt_p]++;
    } );
  }

  /* sort functions */
  std::vector<std::pair<static_truth_table<K>, uint32_t>> functions_sorted;
  functions_sorted.reserve( functions.size() );

  for ( auto const& p : functions )
    functions_sorted.push_back( p );
  
  std::sort( functions_sorted.begin(), functions_sorted.end(), [&]( auto const& a, auto const& b ) {
    return a.second > b.second;
  } );

  /* report to file */
  std::ofstream out( "functions_" + std::to_string( K ) + ".txt" );

  for ( auto const& p : functions_sorted )
  {
    auto tt = p.first;
    out << fmt::format( "{} {} {}\n", kitty::min_base_inplace( tt ).size(), p.second, generate_sop( p.first ) );
  }

  // for ( auto const& p : functions_sorted )
  // {
  //   std::cout << kitty::to_hex( p.first ) << ": " << generate_sop( p.first ) << "\n";
  // }

  out.close();

  return 0;
}