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
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <kitty/kitty.hpp>
#include <experiments.hpp>

int main()
{
  using namespace mockturtle;
  using namespace kitty;

  constexpr uint32_t NInputs = 4;

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{false};
  exact_library_params eps;
  eps.use_dont_cares = true;
  exact_library<mig_network, mig_npn_resynthesis> exact_lib( resyn, eps );

  /* create a function */
  static_truth_table<NInputs> tt, dc;
  create_from_hex_string( tt, "FCA8" );
  create_from_hex_string( dc, "0A0A" );

  /* get the NPN representative of tt */
  auto [tt_npn, neg, perm] = exact_npn_canonization( tt );

  /* permute the DC */
  auto dc_npn = create_from_npn_config2( std::make_tuple( dc, neg & ~( 1 << NInputs ), perm ) );

  neg ^= 1 << NInputs;

  /* report phase and permutation */
  std::cout << "NPN: "; print_hex( tt_npn ); std::cout << "\n";
  std::cout << "DC NPN: "; print_hex( dc_npn ); std::cout << "\n";
  std::cout << fmt::format( "Phase: {0:x}\n", neg );
  std::cout << "Perm : ";
  for ( auto const& p : perm )
    std::cout << fmt::format( "{0:d} ", p );
  std::cout << "\n";

  auto const standard_match = exact_lib.get_supergates( ~tt_npn );
  auto const dc_match = exact_lib.get_supergates( ~tt_npn, dc_npn, neg, perm );

  std::cout << fmt::format( "Standard match size {}\n", (unsigned) standard_match->front().area );
  std::cout << fmt::format( "DC       match size {}\n", (unsigned) dc_match->front().area );

  /* report updated phase and permutation */
  std::cout << fmt::format( "Phase: {0:x}\n", neg );
  std::cout << "Perm : ";
  for ( auto const& p : perm )
    std::cout << fmt::format( "{0:d} ", p );
  std::cout << "\n";

  /* simulate the permutation of the mapper */
  std::array<uint8_t, NInputs> permutation;
  uint32_t phase = 0;

  for ( auto j = 0u; j < perm.size() && j < NInputs; ++j )
  {
    permutation[perm[j]] = j;
    phase |= ( ( neg >> perm[j] ) & 1 ) << j;
  }

  /* dump the network and check the functionality */
  mig_network mig;
  auto const a = mig.create_pi();
  auto const b = mig.create_pi();
  auto const c = mig.create_pi();
  auto const d = mig.create_pi();

  std::vector<uint32_t> best_cut = {1, 2, 3, 4};
  std::vector<typename mig_network::signal> children( NInputs, mig.get_constant( false ) );

  auto ctr = 0u;
  for ( auto l : best_cut )
  {
    children[permutation[ctr++]] = mig.make_signal( mig.index_to_node( l ) );
  }
  for ( auto i = 0u; i < NInputs; ++i )
  {
    if ( ( phase >> i ) & 1 )
    {
      children[i] = !children[i];
    }
  }

  topo_view topo{ exact_lib.get_database(), dc_match->front().root };
  auto f = cleanup_dangling( topo, mig, children.begin(), children.end() ).front();
  mig.create_po( !f );

  /* simulate to verify the correctness of the result */
  default_simulator<kitty::dynamic_truth_table> sim( 4 );
  auto sim_res = simulate_nodes<kitty::dynamic_truth_table>( mig, sim );
  print_hex( sim_res[f] );

  return 0;
}