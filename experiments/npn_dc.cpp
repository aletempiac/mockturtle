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
#include <map>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <kitty/kitty.hpp>
#include <experiments.hpp>

int main()
{
  using namespace mockturtle;
  using namespace kitty;

  constexpr uint32_t NInputs = 4;

  /* library to map to XAGs */
  xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
  exact_library_params ps;
  ps.np_classification = false;
  ps.use_dont_cares = true;
  exact_library<xag_network, decltype( resyn )> exact_lib( resyn, ps );

  /* create a function */
  static_truth_table<NInputs> tt, dc;
  create_from_hex_string( tt, "8000" );
  create_from_hex_string( dc, "8000" );

  /* get the NPN representative of tt */
  auto [tt_npn, phase, perm] = exact_npn_canonization( tt );

  /* permute the DC */
  auto dc_npn = create_from_npn_config( std::make_tuple( dc, phase, perm ) );

  auto const standard_match = exact_lib.get_supergates( tt_npn );
  auto const dc_match = exact_lib.get_supergates( tt_npn, dc_npn, phase, perm );

  std::cout << fmt::format( "Standard match size {}\n", (unsigned) standard_match->front().area );
  std::cout << fmt::format( "DC       match size {}\n", (unsigned) dc_match->front().area );

  return 0;
}