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
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, float, bool, bool> exp(
      "mapper_dc", "benchmark", "size", "size_xag", "size_xag_dc", "depth", "depth_xag", "depth_xag_dc", "runtime1", "runtime2", "equivalent1", "equivalent2" );

  /* library to map to XAGs */
  xag_npn_resynthesis<xag_network> resyn;
  exact_library_params eps;
  exact_library<xag_network, xag_npn_resynthesis<xag_network>> exact_lib( resyn, eps );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    xag_network xag;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( xag ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = xag.num_gates();
    const uint32_t depth_before = depth_view( xag ).depth();

    map_params ps;
    ps.skip_delay_round = true;
    ps.required_time = std::numeric_limits<double>::max();
    map_stats st1;

    xag_network res1 = map( xag, exact_lib, ps, &st1 );

    ps.enable_logic_sharing = true;
    ps.logic_sharing_cut_limit = 1;
    map_stats st2;

    xag_network res2 = map( xag, exact_lib, ps, &st2 );

    const auto cec1 = benchmark == "hyp" ? true : abc_cec( res1, benchmark );
    const auto cec2 = benchmark == "hyp" ? true : abc_cec( res2, benchmark );

    const uint32_t depth_xag1 = depth_view( res1 ).depth();
    const uint32_t depth_xag2 = depth_view( res2 ).depth();

    exp( benchmark, size_before, res1.num_gates(), res2.num_gates(), depth_before, depth_xag1, depth_xag2, to_seconds( st1.time_total ), to_seconds( st2.time_total ), cec1, cec2 );
  }

  exp.save();
  exp.table();

  return 0;
}