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
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>

#include <experiments.hpp>

using namespace mockturtle;

aig_network abc_opt( aig_network const& aig, std::string const& script )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; {}; write_aiger /tmp/tmp.aig\"", script );

  aig_network res = aig.clone();

  while ( true )
  {
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

    aig_network tmp;
    if( lorina::read_aiger( "/tmp/tmp.aig", aiger_reader( tmp ) ) != lorina::return_code::success )
    {
      std::cerr << "read_aiger failed" << std::endl;
      return res;
    }

    if ( tmp.num_gates() < res.num_gates() )
      res = tmp;
    else
      break;
    break;
  }

  return res;
}

int main()
{
  using namespace experiments;

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    if ( aig.num_gates() > 650000 )
      continue;

    // aig_network res = abc_opt( aig, "dfraig; resyn; resyn2; resyn2rs; resyn; resyn2;" );
    aig_network res = abc_opt( aig, "rec_start3 rec6Lib_final_filtered3_recanon.aig; dfraig; resyn; resyn2; resyn2rs; &get; &if -y -K 6; &put; resyn2rs" );
    write_aiger( res, "lms/" + benchmark + ".aig" );
  }

  return 0;
}