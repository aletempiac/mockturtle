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
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <experiments.hpp>

using namespace mockturtle;

std::tuple<uint32_t, uint32_t, uint32_t> abc_map( aig_network const& aig )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read_lut lut1.lib; read /tmp/tmp.aig; if -J 66; ps\"" );

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

  std::cout << result << std::endl;

  /* parse the result */
  uint32_t area = 0, edges = 0, delay = 0;

  std::size_t pos = result.find( "nd" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( "e", pos + 1 ) - pos - 2 );
    lorina::detail::trim( area_res );

    // std::cout << area_res << std::endl;

    area = std::stoll( area_res );

    pos = result.find( "=", pos + 1 );
    std::string edges_res = result.substr( pos + 1, result.find( "l", pos + 1 ) - pos - 2 );

    edges = std::stoll( edges_res );

    pos = result.find( "l", pos + 1 );
    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( " ", pos + 1 ) - pos - 2 );

    delay = std::stoll( delay_res );
  }
  else
  {
    std::cout << "[e] failed to read the result\n";
  }

  return std::make_tuple( area, edges, delay );
}

int main()
{
  using namespace experiments;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, float> exp(
      "ABC_if", "benchmark", "size", "depth", "LUTs", "Edges", "Depth", "Time(s)" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* skip very large designs */
    // if ( benchmark == "hyp" )
    //   continue;

    /* balancing */
    // aig_balance( aig, { false } );

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    /* METHOD 1: map using ABC */
    stopwatch<>::duration time_abc{ 0 };
    auto [area_abc, edges_abc, delay_abc] = call_with_stopwatch( time_abc, [&]() {
      return abc_map( aig );
    } );

    exp( benchmark, size_before, depth_before, area_abc, edges_abc, delay_abc, to_seconds( time_abc ) );
  }

  exp.save();
  exp.table();

  return 0;
}