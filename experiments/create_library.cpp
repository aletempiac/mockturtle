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
#include <deque>
#include <map>

#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_genlib.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/utils/super_utils.hpp>

#include <experiments.hpp>

using namespace mockturtle;

std::vector<gate> filter_gates( std::vector<gate> const& gates, uint32_t max_inputs )
{
  std::vector<gate> filtered_gates;
  filtered_gates.reserve( gates.size() );

  /* remove multi-output gates */
  super_utils<6> super( gates, {}, { false, true } );
  std::deque<composed_gate<6>> single_output_lib = super.get_super_library();
  std::cout << fmt::format( "[i] Filtered multi-output gates; new size {}\n", single_output_lib.size() );

  /* select gates between 2 and max_inputs */
  for ( auto const& sg : single_output_lib )
  {
    if ( sg.num_vars < 2 || sg.num_vars > max_inputs )
      continue;
    filtered_gates.push_back( *sg.root );
  }
  std::cout << fmt::format( "[i] Filtered based on the number of inputs; new size {}\n", filtered_gates.size() );

  /* define functionality classes and select the smallest implementation */
  std::unordered_map<kitty::dynamic_truth_table, std::vector<gate>, kitty::hash<kitty::dynamic_truth_table>> classes;
  for ( gate const& g : filtered_gates )
    classes[g.function].push_back( g );

  std::cout << fmt::format( "[i] Found {} classes\n", classes.size() );
  filtered_gates.clear();

  for ( auto const& pair : classes )
  {
    std::vector<gate> const& func_class = pair.second;
    uint32_t smallest_i = 0;
    double smallest_area = func_class[0].area;
    for ( uint32_t i = 1; i < func_class.size(); ++i )
    {
      if ( func_class[i].area < smallest_area )
      {
        smallest_area = func_class[i].area;
        smallest_i = i;
      }
    }

    filtered_gates.push_back( func_class[smallest_i] );
  }

  std::cout << fmt::format( "[i] Filtered based on functionality classes; new size {}\n", filtered_gates.size() );

  return filtered_gates;
}

std::vector<gate> combine_gates( std::vector<gate> const& gates )
{
  std::vector<gate> combined_gates;
  uint32_t ids = 0;

  for ( uint32_t i = 0; i < gates.size() - 1; ++i )
  {
    /* combinations without repetitions */
    gate const& g_i = gates[i];
    for ( uint32_t j = i + 1; j < gates.size(); ++j )
    {
      gate const& g_j = gates[j];

      if ( g_i.num_vars != g_j.num_vars )
        continue;

      /* create a new gate g_new */
      std::string name = g_i.name + g_j.name;
      double area = ( g_i.area + g_j.area ) * 0.8;
      gate g_oi = g_i;
      gate g_oj = g_j;
      g_oi.name = name;
      g_oj.name = name;
      g_oi.id = ids++;
      g_oj.id = ids++;
      g_oi.area = area;
      g_oj.area = area;

      combined_gates.push_back( g_oi );
      combined_gates.push_back( g_oj );
    }
  }

  std::cout << fmt::format( "[i] Created {} 2-output gates\n", combined_gates.size() );

  return combined_gates;
}

int main()
{
  
  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( "asap7.genlib" );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  /* extract single-output gates only with 2 to 3 inputs */
  std::vector<gate> filtered_gates = filter_gates( gates, 3 );

  /* combine extracted gates */
  std::vector<gate> combined_gates = combine_gates( filtered_gates );

  write_genlib( combined_gates, "combined_lib.genlib" );

  return 0;
}
