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

  /* NPN classes */
  std::unordered_set<static_truth_table<NInputs>, hash<static_truth_table<NInputs>>> classes;
  static_truth_table<NInputs> tt;
  do
  {
    const auto res = exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    next_inplace( tt );
  } while ( !is_const0( tt ) );

  /* library to map to XAGs */
  xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_complete> resyn;
  exact_library_params ps;
  ps.np_classification = false;
  ps.use_dont_cares = true;
  exact_library<xag_network, decltype( resyn )> exact_lib( resyn, ps );

  /* save the size for each NPN class */
  std::unordered_map<static_truth_table<NInputs>, uint32_t, hash<static_truth_table<NInputs>>> class_sizes;
  for ( auto const& entry : classes )
  {
    // print_binary( entry );
    const auto res = exact_lib.get_supergates( entry );
    const unsigned numgates = res ? static_cast<unsigned>( res->front().area ) : 0;
    // std::cout << fmt::format( "\t {}\n", numgates );
    class_sizes.insert( {entry, numgates} );
  }

  /* study don't cares */
  for ( auto entry1 = classes.begin(); entry1 != classes.end(); ++entry1 )
  {
    // print_binary( *entry1 );

    std::unordered_map<static_truth_table<NInputs>, uint32_t, hash<static_truth_table<NInputs>>> dc_sets;
    for ( auto entry2 = classes.begin(); entry2 != classes.end();  ++entry2 )
    {
      uint32_t size = class_sizes[*entry2];

      /* evaluate DC only for size improvement */
      if ( size >= class_sizes[*entry1] )
        continue;

      exact_npn_canonization( *entry2, [&]( auto const& tt ) {
        /* extract the DC set */
        const auto dc = *entry1 ^ tt;

        // if ( count_ones( dc ) > 3 )
        //   return;

        /* check existance */
        if ( auto const& p = dc_sets.find( dc ); p != dc_sets.end() )
        {
          if ( size < std::get<1>( *p ) )
            dc_sets[dc] = size;

          return;
        }

        /* check dominance */
        auto it = dc_sets.begin();
        while ( it != dc_sets.end() )
        {
          auto const& dc_set_tt = std::get<0>( *it );
          auto const& and_tt = dc_set_tt & dc;

          if ( dc_set_tt == and_tt && std::get<1>( *it ) <= size )
          {
            return;
          }
          else if ( dc == and_tt && size <= std::get<1>( *it ) )
          {
            it = dc_sets.erase( it );
          }
          else
          {
            ++it;
          }
        }

        /* insert in the dc_sets */
        dc_sets[dc] = size;
      } );
    }

    uint32_t max_gain = 0;
    float avg_gain = 0;
    for ( auto const& dc : dc_sets )
    {
      max_gain = std::max( max_gain, class_sizes[*entry1] - std::get<1>( dc ) );
      avg_gain += class_sizes[*entry1] - std::get<1>( dc );
    }
    // std::cout << fmt::format( "\t {:>5d}\t {:>5d}\t {:>5.3f}\n", dc_sets.size(), max_gain, avg_gain / dc_sets.size() );
  }

  return 0;
}