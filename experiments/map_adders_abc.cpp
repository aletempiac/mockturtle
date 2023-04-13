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
#include <lorina/genlib.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/experimental/decompose_multioutput.hpp>
#include <mockturtle/algorithms/experimental/emap.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/map_adders.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/dont_touch_view.hpp>

#include <experiments.hpp>

using namespace mockturtle;

std::pair<double, double> abc_map( aig_network const& aig, std::string const& library )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; read {}; &get; &nf -p; &put; print_stats;\"", library );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "ABC: popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  /* parse the result */
  double area = -1, delay = -1;

  std::size_t pos = result.find( "area" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 1 );
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 1 );

    delay = std::stod( delay_res );
  }
  else
  {
    std::cout << "[e] failed to read the result\n";
  }

  return std::make_pair( area, delay );
}

int main()
{
  using namespace experiments;

  experiment<std::string, uint32_t, uint32_t, double, double, float, double, double, uint32_t, float> exp(
      "map_adders_ABC", "benchmark", "size", "depth", "area_abc", "delay_abc", "runtime_abc", "area_emap",  "delay_emap", "used_adders", "runtime_emap" );

  /* library to map to technology */
  fmt::print( "[i] processing technology library\n" );
  std::vector<gate> gates;
  std::string cell_library = "/Users/tempia/Documents/phd/libraries/aletempiac_merge/mockturtle/build/asap7.genlib";
  std::ifstream in( cell_library );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tps.verbose = true;
  tps.load_multioutput_gates = true;
  tps.load_multioutput_gates_single = true;
  tech_library<6, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* balancing */
    aig_balance( aig, { false } );

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    /* METHOD 1: map using ABC */
    stopwatch<>::duration time_abc{ 0 };
    auto [area_abc, delay_abc] = call_with_stopwatch( time_abc, [&]() {
      return abc_map( aig, cell_library );
    } );

    /* METHOD 2: map using emap */
    emap_params ps;
    ps.map_multioutput = true;
    ps.area_oriented_mapping = false;
    emap_stats st;
    binding_view<klut_network> res_emap = emap<aig_network, 6>( aig, tech_lib, ps, &st );

    exp( benchmark, size_before, depth_before, area_abc, delay_abc, to_seconds( time_abc ), st.area, st.delay, st.multioutput_gates, to_seconds( st.time_total ) );
  }

  exp.save();
  exp.table();

  return 0;
}
