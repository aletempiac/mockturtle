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
#include <lorina/verilog.hpp>
#include <kitty/kitty.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/utils/index_list.hpp>

#include <experiments.hpp>

template<class Ntk = mockturtle::aig_network>
void generate( bool allow_xor = false, uint32_t num_candidates = 10u )
{
  using namespace mockturtle;

  /* XAG synthesis */
  exact_aig_resynthesis<Ntk> resyn( allow_xor );
  resyn.set_num_candidates( num_candidates );

  Ntk ntk;

  std::vector<typename Ntk::signal> leaves;

  for ( auto i = 0u; i < 4; ++i )
  {
    leaves.push_back( ntk.create_pi() );
  }

  std::unordered_set<kitty::static_truth_table<4>, kitty::hash<kitty::static_truth_table<4>>> classes;
  kitty::static_truth_table<4> tt;

  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  /* Constuct supergates */
  auto i = 0u;
  for ( auto const& entry : classes )
  {
    auto tt = kitty::shrink_to( entry, 4 );
    resyn( ntk, tt, leaves.begin(), leaves.begin() + 4, [&]( auto const& f ) {
      ntk.create_po( f );
    } );

    /* progress report */
    std::cout << fmt::format( "{} functions synthesized\r", i++ ) << std::flush;
  }
  std::cout << fmt::format( "{} functions synthesized\n", i );

  write_verilog( ntk, "exact_synthesis_aig.v" );
}

void create_index_list()
{
  using namespace mockturtle;

  xag_network ntk;

  if ( lorina::read_verilog( "exact_synthesis_xag.v", verilog_reader( ntk ) ) != lorina::return_code::success )
  {
    return;
  }

  /* create xag index list */
  xag_index_list list( 4 );
  encode( list, ntk );

  /* extract raw */
  auto const raw = list.raw();

  std::ofstream out( "xags_raw.txt" );

  for ( auto const& v : raw )
  {
    out << v << ", ";
  }

  out.close();
}

int main()
{
  generate<mockturtle::aig_network>( false, 50 );
  // create_index_list();
}