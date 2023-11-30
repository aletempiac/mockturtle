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
#include <mockturtle/algorithms/ac_decomposition.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <experiments.hpp>

uint32_t compute_num_edges( mockturtle::klut_network const& klut )
{
  uint32_t edges = 0;
  klut.foreach_gate( [&]( auto const& n ) {
    edges += klut.fanin_size( n );
  } );

  return edges;
}

void run_lut10()
{
  using namespace experiments;
  using namespace mockturtle;

  for ( auto const& benchmark : epfl_benchmarks( experiments::cavlc ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = 10u;
    ps.verbose = true;
    lut_map_stats st;
    const auto klut = lut_map<aig_network>( aig, ps, &st );

    klut.foreach_gate( [&]( auto const& g ) {
      kitty::print_hex( klut.node_function( g ) );
      std::cout << "\n";
    } );
  }
}

void run_mapper()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double, bool> exp( 
    "ACD", "benchmark", "luts", "luts_acd", "lut_depth", "lut_depth_acd", "edges", "edges_acd", "runtime", "runtime_acd", "equivalent"
  );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( "lms/" + benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    
    if ( benchmark == "hyp" )
      continue;
    
    // aig_balance( aig, { false } );

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = 6u;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = false;
    ps.area_share_rounds = 0;
    ps.edge_optimization = true;
    ps.cut_expansion = true;
    ps.verbose = false;
    lut_map_stats st;

    const klut_network klut = lut_map<aig_network>( aig, ps, &st );

    ps.delay_oriented_acd = true;
    ps.relax_required = 0;
    ps.acd_cut_size = 8;
    ps.verbose = true;
    lut_map_stats st_acd;
    const auto klut_acd = lut_map<aig_network, true>( aig, ps, &st_acd );

    uint32_t const luts = klut.num_gates();
    uint32_t const lut_depth = depth_view( klut ).depth();
    uint32_t const edges = compute_num_edges( klut );
    uint32_t const luts_acd = klut_acd.num_gates();
    uint32_t const lut_depth_acd = depth_view( klut_acd ).depth();
    uint32_t const edges_acd = compute_num_edges( klut_acd );
    // auto const cec = benchmark == "hyp" ? true : abc_cec( klut_acd, benchmark );

    exp( benchmark, luts, luts_acd, lut_depth, lut_depth_acd, edges, edges_acd, to_seconds( st.time_total ), to_seconds( st_acd.time_total ), true );

    // write_blif( klut_acd, benchmark + ".blif" );
  }

  exp.save();
  exp.table();
}

void run_lut8()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, double, bool> exp( "lut_mapper", "benchmark", "luts", "lut_depth", "edges", "runtime", "equivalent" );

  std::string benchmark = "test_lut8.aig";
  fmt::print( "[i] processing {}\n", benchmark );
  aig_network aig;
  if ( lorina::read_aiger( benchmark, aiger_reader( aig ) ) != lorina::return_code::success )
  {
    return;
  }

  /* test decomposing */
  kitty::dynamic_truth_table tt( 8 );
  kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

  ac_decomposition_params ac_ps;
  ac_ps.lut_size = 6;

  std::vector<uint32_t> late_arriving = {};
  detail::ac_decomposition_impl acd( tt, 8, ac_ps );
  acd.run( late_arriving );
  auto res = acd.get_result_ntk();

  if ( res.has_value() )
  {
    klut_network klut = *res;
    uint32_t const luts = klut.num_gates();
    uint32_t const lut_depth = depth_view( klut ).depth();
    uint32_t const edges = compute_num_edges( klut );

    bool const cec = abc_cec_impl( klut, benchmark );
    exp( benchmark, luts, lut_depth, edges, 0, cec );
  }

  /* test mapping */
  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = 6u;
  ps.cut_enumeration_ps.cut_limit = 8u;
  ps.recompute_cuts = true;
  ps.area_oriented_mapping = false;
  ps.edge_optimization = true;
  ps.cut_expansion = true;
  // ps.delay_oriented_acd = true;
  // ps.acd_cut_size = 8;
  ps.verbose = true;
  lut_map_stats st;
  const auto klut = lut_map<aig_network, true>( aig, ps, &st );

  depth_view<klut_network> klut_d{ klut };

  auto const cec = true;

  uint32_t const luts = klut.num_gates();
  uint32_t const lut_depth = depth_view( klut ).depth();
  uint32_t const edges = compute_num_edges( klut );

  exp( benchmark, luts, lut_depth, edges, to_seconds( st.time_total ), cec );

  exp.save();
  exp.table();
}

void test_new_enumeration()
{
  using namespace mockturtle;

  kitty::dynamic_truth_table tt( 6 );
  kitty::create_from_hex_string( tt, "1234123412341234" );

  ac_decomposition_params ac_ps;
  ac_ps.lut_size = 6;

  detail::ac_decomposition_impl acd( tt, 6, ac_ps );
  acd.test_enumeration( 3, 2 );
}

int main()
{
  // run_mapper();
  // run_lut10();
  // run_lut8();
  test_new_enumeration();
  return 0;
}
