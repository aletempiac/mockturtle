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
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/xag_resyn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  /* XAG exact */
  xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
  exact_library<xag_network, decltype( resyn )> exact_lib( resyn );

  /* XAG heuristic */
  using engine_abc_t = xag_resyn_decompose<kitty::static_truth_table<12>>;
  engine_abc_t::stats st;

  std::vector<uint32_t> divisors;
  std::vector<kitty::static_truth_table<12>> divisor_functions;
  kitty::static_truth_table<12> x;
  for ( uint32_t i = 0; i < 4; ++i )
  {
    kitty::create_nth_var( x, i );
    divisor_functions.emplace_back( x );
    divisors.emplace_back( i );
  }

  uint32_t cost_exact = 0;
  uint32_t cost_heuristics = 0;
  uint32_t failures = 0;
  double percentage = 0;
  double count_classes = 0;
  uint32_t neq = 0;
  kitty::static_truth_table<4> tt;
  stopwatch<>::duration time_total{ 0 };
  do
  {
    /* exact method */
    const auto res = kitty::exact_npn_canonization( tt );
    auto const gates = exact_lib.get_supergates( std::get<0>( res ) );
    assert( gates != nullptr );
    cost_exact += gates->at( 0 ).area;

    /* heuristic method */
    kitty::static_truth_table<12> target = kitty::extend_to<12>( tt ); 
    kitty::static_truth_table<12> care;
    care = ~care;
    engine_abc_t engine{ st };
    auto const index =  call_with_stopwatch( time_total, [&]() {
      return engine( target, care, divisors.begin(), divisors.end(), divisor_functions );
    } );
    if ( !index.has_value() )
    {
      cost_heuristics += gates->at( 0 ).area;
      ++failures;
    }
    else
    {
      cost_heuristics += index->num_gates();
      if ( index->num_gates() == gates->at( 0 ).area )
        ++percentage;
      
      /* check correctness */
      xag_network xag_res;
      decode( xag_res, *index );

      default_simulator<kitty::static_truth_table<4>> sim;
      const auto tt_out = simulate<kitty::static_truth_table<4>>( xag_res, sim );
      if ( tt_out.front() != tt )
        ++neq;
    }

    ++count_classes;
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  std::cout << fmt::format( "[i] Cost exact     = {}\n[i] Cost heuristic = {}\n", cost_exact, cost_heuristics );
  std::cout << fmt::format( "[i] Percentage     = {:>5.2f}%\n[i] Failures       = {}\n[i] NEQ            = {}\n", percentage / count_classes * 100, failures, neq );
  std::cout << fmt::format( "[i] Time total     = {:>5.3f}\n[i] Average time   = {:>7.7f}\n", to_seconds( time_total ), to_seconds( time_total ) / count_classes );

  return 0;
}
