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
#include <lorina/verilog.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/binate_covering.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/views/mapping_view.hpp>

#include <experiments.hpp>

std::vector<std::string> mcnc_benchmarks = {"5xp1", "c432", "c880", "count", "in5", "k2", "max512", "mlp4" "sqr6", "c1908", "c5315", "chkn", "dist", "in6", "m3", "misex3", "prom2", "x1dn"};
std::vector<std::string> small_benchmarks = {"b_1", "b_2", "b_3", "b_4", "b_5", "b_6", "b_7", "b_8", "b_9", "b_10", "b_11", "b_12" };

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, double, bool> exp( "binate_covering", "benchmark", "luts", "runtime", "equivalent" );

  for ( auto const& benchmark : small_benchmarks )
  {
    // if ( benchmark == "c5315" )
    //   break;

    aig_network aig;
    if ( lorina::read_aiger( benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }
    // if ( lorina::read_verilog( benchmark_path_verilog( benchmark ), verilog_reader( aig ) ) != lorina::return_code::success )
    // {
    //   continue;
    // }
    fmt::print( "[i] processing {}\t num gates {}\n", benchmark, aig.num_gates() );

    binate_covering_params ps;
    binate_covering_stats st;
    ps.cut_enumeration_ps.cut_size = 6;
    ps.cut_enumeration_ps.cut_limit = 4;
    // ps.bound = 60;
    ps.timeout = 5;
    ps.verbose = true;
    ps.debug = true;

    mapping_view<aig_network, false> mapped_aig{ aig };
    binate_covering<decltype( mapped_aig ), false>( mapped_aig, ps, &st );
    const auto klut = *collapse_mapped_network<klut_network>( mapped_aig );

    bool const cec = abc_cec_impl( klut, benchmark + ".aig" );
    // bool const cec = abc_cec_verilog( klut, benchmark );
    // bool const cec = true;

    exp( benchmark, klut.num_gates(), to_seconds( st.time_total ), cec );

    // break;
  }

  exp.save();
  exp.table();

  return 0;
}