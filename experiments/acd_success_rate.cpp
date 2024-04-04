/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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

#include <chrono>
#include <random>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/acd.hpp>
#include <mockturtle/algorithms/acd66.hpp>
#include <mockturtle/algorithms/acd666.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/s66.h>
#include <mockturtle/algorithms/spfd_utils.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/truth_table_cache.hpp>

#include <experiments.hpp>

std::tuple<uint32_t, uint32_t, uint32_t> abc_map( std::string const& tt, std::string const& map_flag, uint32_t cut_size )
{
  std::string command = fmt::format( "abc -q \"read_truth {}; if -{} 66 -K {}; ps\"", tt, map_flag, cut_size );

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

  // std::cout << result << std::endl;

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

bool abc_acd( std::string const& tt_string )
{
  // std::string command = fmt::format( "./cascade/s66dec31112 {}", tt );

  // std::array<char, 128> buffer;
  // std::string result;
  // std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  // if ( !pipe )
  // {
  //   throw std::runtime_error( "ABC: popen() failed" );
  // }
  // while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  // {
  //   result += buffer.data();
  // }

  // std::cout << result << std::endl;

  // return result[0] == 'S';
  int nVars = std::log2( 4 * tt_string.size() );

  kitty::dynamic_truth_table tt( nVars );
  kitty::create_from_hex_string( tt, tt_string );
  word Truth[CLU_WRD_MAX] = { 0 };

  for ( auto i = 0u; i < tt.num_blocks(); ++i )
  {
    Truth[i] = tt._bits[i];
  }
  word Func0, Func1, Func2;
  If_Grp_t G1 = { 0 }, G2 = { 0 }, R = { 0 };
  int nVarsNew = nVars;              // the number of variables afer support minimization
  int pVarPerm[CLU_VAR_MAX] = { 0 }; // the remaining variables after support minimization
  G1 = If_CluCheckTest( 2, 6, Truth, nVars, &R, &G2, &Func0, &Func1, &Func2, &nVarsNew, pVarPerm );
  return G1.nVars > 0;
}

bool mockturtle_acd66( std::string const& tt_string, uint32_t delay_profile )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );

  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  uint64_t words[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];

  acd66_impl acd( tt.num_vars(), true, false );

  int res = acd.run( words, delay_profile );

  if ( res == 0 )
    return false;

  // int correct = acd.compute_decomposition();

  // if ( correct == 1 )
  // {
  //   std::cout << fmt::format( "[e] incorrect decomposition of {}\n", tt_string );
  // }

  return true;
}

bool mockturtle_acd666( std::string const& tt_string )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );

  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  uint64_t words[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];

  acd666_impl acd( tt.num_vars(), false );

  bool res = acd.run( words );

  if ( !res )
    return false;

  int correct = acd.compute_decomposition();

  if ( correct == 1 )
  {
    std::cout << fmt::format( "[e] incorrect decomposition of {}\n", tt_string );
  }

  return true;
}

int32_t mockturtle_acd_generic( std::string const& tt_string, uint32_t delay_profile )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );

  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  uint64_t words[1024];

  for ( uint32_t i = 0; i < tt.num_blocks(); ++i )
    words[i] = tt._bits[i];

  acd_params ps;
  ps.use_first = false;
  ps.max_multiplicity = 16;
  acd_stats st;
  acd_impl acd( num_vars, ps, &st );

  int res = acd.run( words, delay_profile );

  if ( res < 0 )
    return -1;

  // int correct = acd.compute_decomposition();

  // std::vector<acd_result> const& dec_res = acd.get_acd_result();

  // /* simulate and check decomposition */
  // klut_network ntk;
  // for ( uint32_t i = 0; i < num_vars; ++i )
  //   ntk.create_pi();
  // for ( acd_result const& lut : dec_res )
  // {
  //   std::vector<klut_network::signal> children;
  //   std::transform( lut.support.begin(), lut.support.end(), std::back_inserter( children ), [&]( uint32_t a ) { return a + 2; } );
  //   ntk.create_node( children, lut.tt );
  // }
  // ntk.create_po( ntk.get_node( dec_res.size() + 1 + num_vars ) );

  // default_simulator<kitty::dynamic_truth_table> sim( num_vars );
  // auto sim_tt = simulate<kitty::dynamic_truth_table>( ntk, sim )[0];
  // if ( sim_tt != tt )
  // {
  //   std::cout << fmt::format( "[e] incorrect decomposition of {}\n", tt_string );
  //   std::cout << fmt::format( "[e] obtained  decomposition    {}\n", kitty::to_hex( sim_tt ) );
  // }

  // if ( correct == -1 )
  // {
  //   std::cout << fmt::format( "[e] incorrect decomposition of {}\n", tt_string );
  // }

  return st.num_luts;
}

bool acd_andrea( std::string const& tt_string )
{
  using namespace mockturtle;

  uint32_t num_vars = std::log2( 4 * tt_string.size() );
  kitty::dynamic_truth_table tt( num_vars );
  kitty::create_from_hex_string( tt, tt_string );

  lut_resynthesis_t<6, 10> acd;
  auto ll_tt = acd.decompose( tt, 20 );

  if ( ll_tt && acd.num_luts() > 2 )
    return false;
  return true;
}

void compute_functions( uint32_t cut_size )
{
  using namespace experiments;
  using namespace mockturtle;

  truth_table_cache<kitty::dynamic_truth_table> cache( 200000 );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    lut_map_params ps;
    ps.cut_enumeration_ps.cut_size = cut_size;
    ps.cut_enumeration_ps.cut_limit = 8u;
    ps.area_share_rounds = 0;
    ps.recompute_cuts = true;
    ps.cut_expansion = false;
    lut_map_stats st;

    detail::lut_map_impl<aig_network, true, lut_unitary_cost> p( aig, ps, st );
    const auto klut = p.run();

    truth_table_cache<kitty::dynamic_truth_table> const& cache_map = p.get_truth_cache();

    /* load content into cache */
    for ( uint32_t i = 0; i < cache_map.size(); ++i )
    {
      kitty::dynamic_truth_table tt = cache_map[i << 1];

      if ( tt.num_vars() != cut_size )
        continue;

      // auto res = kitty::exact_n_canonization( tt );
      auto res = kitty::sifting_npn_canonization( tt );
      cache.insert( std::get<0>( res ) );
    }
  }

  /* print truth tables to file */
  std::string filename = "cuts_" + std::to_string( cut_size ) + ".txt";
  std::ofstream out( filename );

  for ( uint32_t i = 0; i < cache.size(); ++i )
  {
    kitty::dynamic_truth_table tt = cache[i << 1];

    kitty::print_hex( tt, out );
    out << "\n";
  }

  out.close();
}

void compute_success_rate( uint32_t cut_size )
{
  /* read file */
  std::ifstream in( "cuts_" + std::to_string( cut_size ) + ".txt" );

  /* count the number of lines */
  uint32_t num_lines = 0;
  std::string line;
  while ( std::getline( in, line ) )
  {
    ++num_lines;
  }

  in.close();
  in.open( "cuts_" + std::to_string( cut_size ) + ".txt" );

  std::ofstream out( "cuts_" + std::to_string( cut_size ) + "_fail.txt" );

  using clock = typename std::chrono::steady_clock;
  using time_point = typename clock::time_point;

  time_point time_begin = clock::now();

  /* compute */
  uint32_t successS = 0, successJ = 0, successJ2 = 0, successG = 0, successA = 0;
  uint32_t visit = 0;
  uint32_t num_luts_acd = 0;
  while ( in.good() )
  {
    // std::cout << fmt::format( "[i] Progress {:8d} / {}\r", visit, num_lines );
    std::string tt;
    in >> tt;

    ++visit;

    if ( tt.size() < 16 )
      continue;

    /* run evaluation */
    // bool resS = abc_acd( tt );
    bool resS = false;
    bool resJ = false;
    // bool resJ2 = mockturtle_acd66( tt, 0 );
    bool resJ2 = false;
    int resG = mockturtle_acd_generic( tt, 0 );
    // int resG = false;
    bool resA = false;

    if ( resG > 0 )
      num_luts_acd += resG;

    if ( resS )
      ++successS;
    if ( resJ )
      ++successJ;
    // else
    //   out << tt << "\n";
    if ( resJ2 )
      ++successJ2;
    if ( resG > 0 )
      ++successG;
    if ( resA )
      ++successA;
  }

  std::cout << "\n";

  /* print stats */
  std::cout << fmt::format( "[i] Run a total of {} truth tables on {} variables\n", num_lines, cut_size );
  std::cout << fmt::format( "[i] Success of -S 66  = {} \t {:>5.2f}%\n", successS, ( (double)successS ) / num_lines * 100 );
  std::cout << fmt::format( "[i] Success of -J 66  = {} \t {:>5.2f}%\n", successJ, ( (double)successJ ) / num_lines * 100 );
  std::cout << fmt::format( "[i] Success of -J 666 = {} \t {:>5.2f}%\n", successJ2, ( (double)successJ2 ) / num_lines * 100 );
  std::cout << fmt::format( "[i] Success of -Z 6   = {} \t {:>5.2f}% \t {} luts\n", successG, ( (double)successG ) / num_lines * 100, num_luts_acd );
  std::cout << fmt::format( "[i] Success of -A 6   = {} \t {:>5.2f}%\n", successA, ( (double)successA ) / num_lines * 100 );
  std::cout << fmt::format( "[i] Time = {:>5.2f} s\n", std::chrono::duration_cast<std::chrono::duration<double>>( clock::now() - time_begin ).count() );

  in.close();
  out.close();
}

void compute_success_rate_delay( uint32_t cut_size, uint32_t late_vars = 2, uint32_t const repeat = 10 )
{
  /* read file */
  std::ifstream in( "cuts_" + std::to_string( cut_size ) + ".txt" );

  /* count the number of lines */
  uint32_t num_lines = 0;
  std::string line;
  while ( std::getline( in, line ) )
  {
    ++num_lines;
  }

  in.close();
  in.open( "cuts_" + std::to_string( cut_size ) + ".txt" );

  std::ofstream out( "cuts_" + std::to_string( cut_size ) + "_fail.txt" );

  using clock = typename std::chrono::steady_clock;
  using time_point = typename clock::time_point;

  time_point time_begin = clock::now();
  std::default_random_engine::result_type seed{ 1 };
  std::uniform_int_distribution<uint32_t> dist( 0u, cut_size - 1 );

  /* compute */
  uint32_t successJ = 0, successG = 0;
  uint32_t visit = 0;
  while ( in.good() )
  {
    std::cout << fmt::format( "[i] Progress {:8d} / {}\r", visit, num_lines );
    std::string tt;
    in >> tt;

    ++visit;

    if ( tt.size() < 16 )
      continue;
    
    /* generate random delay profile with late variables */
    for ( uint32_t i = 0; i < repeat; ++i )
    {
      uint32_t delay_profile = 0;
      for ( uint32_t i = 0; i < late_vars; ++i )
      {
        std::default_random_engine gen( seed++ );
        uint32_t var = dist( gen );

        while ( ( delay_profile >> var ) & 1 )
        {
          std::default_random_engine gen2( seed++ );
          var = dist( gen2 );
        }

        delay_profile |= 1u << var;
      }

      /* run evaluation */
      bool resJ = mockturtle_acd66( tt, delay_profile );
      bool resG = mockturtle_acd_generic( tt, delay_profile ) > 0;

      if ( resJ )
        ++successJ;
      if ( resG )
        ++successG;
    }
  }

  std::cout << "\n";

  /* print stats */
  std::cout << fmt::format( "[i] Run a total of {} truth tables on {} variables\n", num_lines, cut_size );
  std::cout << fmt::format( "[i] Success of -J 66  = {} \t {:>5.2f}%\n", successJ, ( (double)successJ ) / num_lines * 100 / repeat );
  std::cout << fmt::format( "[i] Success of -Z 6   = {} \t {:>5.2f}%\n", successG, ( (double)successG ) / num_lines * 100 / repeat );
  std::cout << fmt::format( "[i] Time = {:>5.2f} s\n", std::chrono::duration_cast<std::chrono::duration<double>>( clock::now() - time_begin ).count() );

  in.close();
  out.close();
}

int main( int argc, char** argv )
{
  if ( argc != 2 )
    return -1;

  uint32_t cut_size = atoi( argv[1] );

  // compute_functions( cut_size );
  compute_success_rate( cut_size );
  // compute_success_rate_delay( cut_size, 4, 10 );

  // mockturtle_acd66( "003f4cd9c000264cffffb3b300004c4c", 0 );
  // mockturtle_acd_generic( "00000000000000000000000000000000000000000000000000000000000000000000000000000000001717ffffffffff001717ffffffffff0000000000000000", 81 );

  return 0;
}