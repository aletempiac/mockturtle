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
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/stopwatch.hpp>

#include <experiments.hpp>

using namespace mockturtle;

std::tuple<uint32_t, uint32_t, uint32_t> abc_opt( klut_network const& klut, std::string const& script )
{
  write_blif( klut, "/tmp/tmp.blif" );
  std::string command = fmt::format( "abc -q \"read_blif /tmp/tmp.blif; {}; write_blif /tmp/tmp.blif; ps\"", script );

  uint32_t area, edges = 0, delay = 0;
  area = klut.num_gates();
  uint32_t area_before;
  do
  {
    area_before = area;

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
    area = 0;
    edges = 0;
    delay = 0;

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
  } while ( area < area_before );

  return std::make_tuple( area, edges, delay );
}

int main( int argc, char** argv )
{
  using namespace experiments;

  if ( argc != 2 )
    return 1;
  
  std::string benchmark{ argv[1] };

  fmt::print( "[i] processing {}\n", benchmark );

  klut_network klut;
  if ( lorina::read_blif( benchmark, blif_reader( klut ) ) != lorina::return_code::success )
  {
    return 2;
  }

  /* METHOD 1: map using ABC */
  stopwatch<>::duration time_abc{ 0 };
  auto [area_abc, edges_abc, delay_abc] = call_with_stopwatch( time_abc, [&]() {
    return abc_opt( klut, "mfs2; &get -nm; &satlut -d -N 64 -C 5000; &put; lutpack" );
  } );

  return 0;
}