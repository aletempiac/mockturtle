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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <mockturtle/algorithms/acd66.hpp>
#include <mockturtle/algorithms/acd666.hpp>
#include <mockturtle/algorithms/s66.h>

#include <kitty/dynamic_truth_table.hpp>

#include <experiments.hpp>

using namespace mockturtle;

float compute_delay( klut_network const& klut, bool skip_buffers = false )
{
  std::vector<float> delays( klut.size() );

  for ( uint32_t i = 0; i < klut.num_pis() + 2; ++i )
  {
    delays[i] = 0;
  }

  float max_delay = 0;
  klut.foreach_gate( [&]( auto const& g ) {
    float pin_delay = 0;
    klut.foreach_fanin( g, [&]( auto const& f ) {
      pin_delay = std::max( pin_delay, delays[f] );
    } );

    if ( skip_buffers && klut.fanin_size( g ) == 1 )
    {
      delays[g] = pin_delay;
      return;
    }

    // delays[g] = klut.fanin_size( g ) > 6 ? pin_delay + 1.2 : pin_delay + 1;
    delays[g] = pin_delay + 1;
    max_delay = std::max( max_delay, delays[g] );
  } );

  return max_delay;
}

uint32_t abc_acd( kitty::dynamic_truth_table const& tt )
{
  word Truth[CLU_WRD_MAX] = { 0 };

  for ( auto i = 0u; i < tt.num_blocks(); ++i )
  {
    Truth[i] = tt._bits[i];
  }
  word Func0, Func1, Func2;
  If_Grp_t G1 = { 0 }, G2 = { 0 }, R = { 0 };
  int nVarsNew = tt.num_vars();              // the number of variables afer support minimization
  int pVarPerm[CLU_VAR_MAX] = { 0 }; // the remaining variables after support minimization
  G1 = If_CluCheckTest( 2, 6, Truth, nVarsNew, &R, &G2, &Func0, &Func1, &Func2, &nVarsNew, pVarPerm );
  return tt.num_vars() + 1 + ( ( G1.nMyu > 2 ) ? 1 : 0 );
}

uint32_t mockturtle_acd66( kitty::dynamic_truth_table const& tt )
{
  uint64_t words[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];

  {
    acd66_impl acd( tt.num_vars(), false, false );

    bool res = acd.run( words );

    if ( res )
    {
      acd.compute_decomposition();
      return acd.get_num_edges();
    }
  }

  {
    acd66_impl acd( tt.num_vars(), true, false );

    acd.run( words );
    acd.compute_decomposition();
    return acd.get_num_edges();
  }

  assert( false );
  return 0;
}

std::tuple<uint32_t, uint32_t, float> compute_edges_for_s66( bool skip_buffers = false )
{
  klut_network klut;
  read_blif( "/tmp/tmp.blif", blif_reader( klut ) );

  /* compute edges */
  uint32_t num_luts = 0;
  uint32_t num_edges = 0;
  klut.foreach_gate( [&]( auto const& n ) {
    if ( skip_buffers && klut.fanin_size( n ) == 1 )
    {
      return;
    }

    ++num_luts;
    if ( klut.fanin_size( n ) > 6 )
    {
       num_edges += mockturtle_acd66( klut.node_function( n ) );
      //  num_edges += abc_acd( klut.node_function( n ) );
       ++num_luts;
    }
    else
    {
      num_edges += klut.fanin_size( n );
    }
  } );

  return std::make_tuple( num_luts, num_edges, compute_delay( klut, skip_buffers ) );
}

std::tuple<uint32_t, uint32_t, float> abc_map( aig_network const& aig, bool use_choices )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  // std::string command = fmt::format( "abc -q \"read_lut lut1.lib; read /tmp/tmp.aig; if -J 666; ps\"" );
  std::string command;

  if ( use_choices )
  {
    command = fmt::format( "abc -q \"read /tmp/tmp.aig; dch; if -Z 6 -K 8; ps\"" );
  }
  else
  {
    command = fmt::format( "abc -q \"read /tmp/tmp.aig; if -Z 6 -K 8; ps\"" );
  }

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "ABC: popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std::cout << result << std::endl;

  /* parse the result */
  uint32_t area = 0, edges = 0;
  float delay = 0;

  std::size_t pos = result.find( "nd" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( "e", pos + 1 ) - pos - 2 );
    lorina::detail::trim( area_res );

    // std::cout << area_res << std::endl;

    area = std::stoll( area_res );

    pos = result.find( "=", pos + 1 );
    std::string edges_res = result.substr( pos + 1, result.find( "l", pos + 1 ) - pos - 2 );

    edges = std::stoll( edges_res );

    pos = result.find( "l", pos + 1 );
    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 2 );

    delay = static_cast<float>( std::stoll( delay_res ) );
  }
  else
  {
    std::cout << "[e] failed to read the result\n";
  }

  return std::make_tuple( area, edges, delay );
}

int main()
{
  using namespace experiments;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, float, float> exp(
      "ABC_if", "benchmark", "size", "depth", "LUTs", "Edges", "Delay", "Time(s)" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    // if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    // {
    //   continue;
    // }
    if ( lorina::read_aiger( "lms/" + benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* balancing */
    // aig_balance( aig, { false } );

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    /* METHOD 1: map using ABC */
    stopwatch<>::duration time_abc{ 0 };
    auto [area_abc, edges_abc, delay_abc] = call_with_stopwatch( time_abc, [&]() {
      return abc_map( aig, benchmark != "hyp" );
    } );

    // auto correction = compute_edges_for_s66( false );
    // area_abc = std::get<0>( correction );
    // edges_abc = std::get<1>( correction );
    // delay_abc = std::get<2>( correction );

    exp( benchmark, size_before, depth_before, area_abc, edges_abc, delay_abc, to_seconds( time_abc ) );
  }

  exp.save();
  exp.table();

  return 0;
}