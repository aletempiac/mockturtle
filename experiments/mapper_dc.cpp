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
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool, bool> exp(
      "mapper_dc", "benchmark", "size", "size_mig", "size_mig_dc", "depth", "depth_mig", "depth_mig_dc", "runtime1", "runtime2", "equivalent1", "equivalent2" );

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{false};
  exact_library_params eps;
  eps.use_dont_cares = true;
  exact_library<mig_network, mig_npn_resynthesis> exact_lib( resyn, eps );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    if ( benchmark != "c880" )
      continue;

    fmt::print( "[i] processing {}\n", benchmark );
    mig_network mig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = mig.num_gates();
    const uint32_t depth_before = depth_view( mig ).depth();

    map_params ps;
    ps.skip_delay_round = true;
    ps.use_dont_cares = true;
    ps.cut_enumeration_ps.minimize_truth_table = true;
    ps.required_time = std::numeric_limits<double>::max();
    map_stats st1;

    mig_network res1 = map( mig, exact_lib, ps, &st1 );

    map_stats st2;
    ps.use_dont_cares = true;
    ps.verbose = true;
    mig_network res2 = map( mig, exact_lib, ps, &st2 );

    const auto cec1 = benchmark == "hyp" ? true : abc_cec( res1, benchmark );
    const auto cec2 = benchmark == "hyp" ? true : abc_cec( res2, benchmark );

    const uint32_t depth_mig1 = depth_view( res1 ).depth();
    const uint32_t depth_mig2 = depth_view( res2 ).depth();

    exp( benchmark, size_before, res1.num_gates(), res2.num_gates(), depth_before, depth_mig1, depth_mig2, to_seconds( st1.time_total ), to_seconds( st2.time_total ), cec1, cec2 );
  }

  exp.save();
  exp.table();

  return 0;
}