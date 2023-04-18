/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/aig_collapse.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/multi_aig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, uint32_t, uint32_t, double, bool> exp(
      "lut_mapper_d", "benchmark", "size", "depth", "size_c", "depth_c", "luts", "lut_depth", "time", "luts_d", "luts_depth_d", "time_d", "equivalent_d" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( benchmark == "hyp" )
      continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    aig_balance( aig, { false } );

    uint32_t const initial_size = aig.num_gates();
    uint32_t const initial_depth = depth_view( aig ).depth();

    /* collapse logic */
    aig_collapse_params cps;
    cps.collapse_limit = 4u;
    multi_aig_network multi_aig = aig_collapse( aig, cps );

    uint32_t const collapsed_size = multi_aig.num_gates();
    uint32_t const collapsed_depth = depth_view( multi_aig ).depth();

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = 6u;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.recompute_cuts = true;
    ps.area_oriented_mapping = false;
    ps.cut_expansion = true;

    /* FLOW1: map AIG */
    lut_map_stats st1;
    mapping_view<aig_network, false> mapped_aig{ aig };
    lut_map<decltype( mapped_aig ), false>( mapped_aig, ps, &st1 );
    const auto klut1 = *collapse_mapped_network<klut_network>( mapped_aig );

    /* FLOW2: map collapsed AIG */
    lut_map_stats st2;
    mapping_view<multi_aig_network, false> mapped_multi_aig{ multi_aig };
    lut_map<decltype( mapped_multi_aig ), false>( mapped_multi_aig, ps, &st2 );
    const auto klut2 = *collapse_mapped_network<klut_network>( mapped_multi_aig );

    uint32_t const flow1_luts = klut1.num_gates();
    uint32_t const flow1_depth = depth_view( klut1 ).depth();
    uint32_t const flow2_luts = klut2.num_gates();
    uint32_t const flow2_depth = depth_view( klut2 ).depth();

    auto const cec = benchmark == "hyp" ? true : abc_cec( klut2, benchmark );

    exp( benchmark, initial_size, initial_depth, collapsed_size, collapsed_depth, flow1_luts, flow1_depth, to_seconds( st1.time_total ), flow2_luts, flow2_depth, to_seconds( st2.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
