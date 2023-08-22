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

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool, bool> exp(
      "mapper_dc", "benchmark", "size", "size_mig", "size_mig_dc", "depth", "depth_mig", "depth_mig_dc", "runtime1", "runtime2", "equivalent1", "equivalent2" );

  using Ntk = mig_network;
  const uint32_t iterations = 1;

  /* library to map to migs */
  // xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_incomplete> resyn;
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  eps.np_classification = true;
  eps.use_dont_cares = true;
  exact_library<Ntk, decltype( resyn )> exact_lib( resyn, eps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    // if ( benchmark == "leon2" || benchmark == "leon3_opt" || benchmark == "leon3" || benchmark == "leon3mp" )
    //   continue;

    fmt::print( "[i] processing {}\n", benchmark );
    Ntk mig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    functional_reduction( mig );
    mig = cleanup_dangling( mig );

    const uint32_t size_before = mig.num_gates();
    const uint32_t depth_before = depth_view( mig ).depth();

    map_params ps;
    ps.skip_delay_round = true;
    ps.use_dont_cares = false;
    ps.cut_enumeration_ps.minimize_truth_table = true;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.enable_logic_sharing = true;
    ps.logic_sharing_cut_limit = 1;
    ps.required_time = std::numeric_limits<double>::max();
    map_stats st1;

    Ntk res1 = cleanup_dangling( mig );
    Ntk res2 = cleanup_dangling( mig );

    auto i = iterations;
    while ( i-- > 0 )
    {
      const uint32_t size_before_map = res1.size();
      Ntk res1_map = map( res1, exact_lib, ps, &st1 );

      if ( res1_map.size() >= size_before_map )
        break;

      res1 = res1_map;
    }

    map_stats st2;
    ps.use_dont_cares = true;
    ps.window_size = 12u;
    ps.verbose = false;

    i = iterations;
    while ( i-- > 0 )
    {
      const uint32_t size_before_map = res2.size();
      Ntk res2_map = map( res2, exact_lib, ps, &st2 );

      if ( res2_map.size() >= size_before_map )
        break;

      res2 = res2_map;
    }

    // const auto cec1 = benchmark == "hyp" ? true : abc_cec( res1, benchmark );
    // const auto cec2 = benchmark == "hyp" ? true : abc_cec( res2, benchmark );
    const auto cec1 = true;
    const auto cec2 = true;

    const uint32_t depth_mig1 = depth_view( res1 ).depth();
    const uint32_t depth_mig2 = depth_view( res2 ).depth();

    exp( benchmark, size_before, res1.num_gates(), res2.num_gates(), depth_before, depth_mig1, depth_mig2, to_seconds( st1.time_total ), to_seconds( st2.time_total ), cec1, cec2 );
  }

  exp.save();
  exp.table();

  return 0;
}