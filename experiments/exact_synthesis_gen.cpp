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
#include <iostream>
#include <fstream>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <kitty/kitty.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/io/write_verilog.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  exact_mig_resynthesis_params ps;
  ps.num_candidates = 1;

  exact_mig_resynthesis resyn( ps );
  resyn.set_bounds( 1, 6 );

  mig_network ntk;

  std::vector<mig_network::signal> leaves;

  for ( auto i = 0u; i < 6; ++i )
  {
    leaves.push_back( ntk.create_pi() );
  }

  std::ifstream in( "functions_merge.txt" );

  std::string tt_string;
  auto i = 0u;
  while ( in >> tt_string )
  {
    kitty::dynamic_truth_table tt_6( 6 );
    kitty::create_from_hex_string( tt_6, tt_string );

    /* shrink the truth table to the support */
    const auto support = kitty::min_base_inplace( tt_6 );
    auto tt = kitty::shrink_to( tt_6, support.size() );

    if ( support.size() == 6 )
      continue;

    /* generate the exact synthesis problem */
    resyn( ntk, tt, leaves.begin(), leaves.begin() + support.size(), [&]( auto const& f ) {
      ntk.create_po( f );
      return true;
    } );

    /* progress report */
    std::cout << fmt::format( "{} functions synthesized\r", i++ );
  }

  in.close();

  std::cout << fmt::format( "{} functions synthesized\n", i );

  write_verilog( ntk, "exact_synthesis_mig.v" );

  return 0;
}