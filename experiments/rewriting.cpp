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

#include <mockturtle/algorithms/rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <lorina/aiger.hpp>

#include <experiments.hpp>
#include <fmt/format.h>
#include <string>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, float, bool>
    exp( "refactoring", "benchmark", "size_before", "size_after", "runtime", "equivalent" );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  eps.use_dont_cares = false;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    // if ( benchmark != "bar" )
    //   continue;

    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    rewriting_params ps;
    rewriting_stats st;
    ps.use_dont_cares = false;
    ps.use_mffc = true;
    ps.allow_multiple_structures = true;
    ps.progress = false;
    ps.verbose = false;

    fanout_view fanout_aig{aig};
    uint32_t const size_before = fanout_aig.num_gates();
    rewrite( fanout_aig, exact_lib, ps, &st );
    aig = cleanup_dangling( aig );

    bool const cec = benchmark == "hyp" ? true : abc_cec( fanout_aig, benchmark );
    exp( benchmark, size_before, aig.num_gates(), to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();

  return 0;
}
