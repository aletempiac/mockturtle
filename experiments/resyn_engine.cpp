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
#include <mockturtle/algorithms/experimental/resyn_engine.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float, bool> exp(
      "emap", "benchmark", "size_before", "size_after", "depth_before", "depth_after", "success", "fail", "runtime", "cec" );

  fmt::print( "[i] processing technology library\n" );

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    /* remap */
    resyn_engine_params rps;
    resyn_engine_stats rst;
    resyn_engine( aig, rps, &rst );

    const uint32_t size_after = aig.num_gates();
    const uint32_t depth_after = depth_view( aig ).depth();

    /* decompose multi-output cells for verification purposes */
    // klut_network klut = decompose_multioutput<block_network, klut_network>( res );
    const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );

    exp( benchmark, size_before, size_after, depth_before, depth_after, rst.num_success, rst.num_fail, to_seconds( rst.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}