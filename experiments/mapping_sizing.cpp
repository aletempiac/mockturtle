/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2024  EPFL
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
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/utils/name_utils.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/names_view.hpp>

#include <experiments.hpp>

using namespace mockturtle;

std::pair<double, double> abc_map( aig_network const& aig, std::string const& library )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; read {}; &get; &nf -p; &put; print_stats; write_verilog /tmp/tmp.v\"", library );

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
    std::string area_res = result.substr( pos + 1, result.find( " ", pos + 2 ) - pos - 1 );
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( " ", pos + 2 ) - pos - 1 );

    delay = std::stod( delay_res );
  }
  else
  {
    std::cout << "[e] failed to read the result\n";
  }

  return std::make_pair( area, delay );
}

std::pair<double, double> abc_size( std::string const& liberty )
{
  std::string command = fmt::format( "abc -q \"read_lib {}; read -m /tmp/tmp.v; buffer; upsize; dnsize; stime\"", liberty );

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

  std::size_t pos = result.find( "Area" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( "(", pos + 1 ) - pos - 1 );
    std::cout << area_res << std::endl;
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( "p", pos + 1 ) - pos - 1 );
    std::cout << delay_res << std::endl;
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

  experiment<std::string, uint32_t, uint32_t, double, double, double, double> exp(
      "mapping_sizing", "benchmark", "size", "depth", "area_abc", "delay_abc", "area_emap", "delay_emap" );

  /* library to map to technology */
  fmt::print( "[i] processing technology library\n" );
  std::string library = "tsmc28";
  // std::string liberty = "/Users/tempia/Documents/phd/libraries/asap7_lib/lib/asap7_merged.lib";
  std::string liberty = "/Users/tempia/Documents/phd/libraries/tsmc28/tcbn28hpcplusbwp30p140ffg0p88v0c.lib";
  // std::string liberty = "/Users/tempia/Documents/phd/libraries/arm28nm/sc9mcpp140z_cln28ht_base_svt_c30_ffg_cbestt_min_0p77v_0c.lib";
  std::string cell_library = cell_libraries_path( library );
  std::vector<gate> gates;
  std::ifstream in( cell_library );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tps.ignore_symmetries = false; // set to true to drastically speed-up mapping with minor delay increase
  tps.verbose = true;
  tech_library<9> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* remove structural redundancies */
    aig_balancing_params bps;
    bps.minimize_levels = false;
    bps.fast_mode = true;
    aig_balance( aig, bps );

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    if ( size_before > 100000 )
      continue;

    /* METHOD 1: map using ABC */
    abc_map( aig, cell_library );
    auto [area_abc, delay_abc] = abc_size( liberty );

    /* METHOD 2: map using emap */
    emap_params ps;
    ps.matching_mode = emap_params::hybrid;
    ps.area_oriented_mapping = false;
    ps.map_multioutput = false;
    ps.use_match_alternatives = true;
    ps.relax_required = 0;
    emap_stats st;
    cell_view<block_network> res = emap<9>( aig, tech_lib, ps, &st );
    write_verilog_with_cell( res, "/tmp/tmp.v" );
    auto [area_emap, delay_emap] = abc_size( liberty );
    // double area_emap = res.compute_area();
    // double delay_emap = res.compute_worst_delay();

    exp( benchmark, size_before, depth_before, area_abc, delay_abc, area_emap, delay_emap );
  }

  exp.save();
  exp.table();

  return 0;
}
