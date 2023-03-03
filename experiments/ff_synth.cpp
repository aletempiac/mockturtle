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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/det_randomization.hpp>
#include <mockturtle/algorithms/factor_resub.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <experiments.hpp>

using namespace mockturtle;

template<class Ntk>
uint32_t count_literals( Ntk& ntk )
{
  ntk.clear_values();
  ntk.foreach_po( [&]( auto const& f ) {
    ntk.incr_value( ntk.get_node( f ) );
  } );

  uint32_t lits = 0;
  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_constant( n ) )
    {
      return;
    }
    else if ( ntk.is_pi( n ) || ntk.fanout_size( n ) > 1 )
    {
      lits += ntk.fanout_size( n ) - ntk.value( n );
      if ( ntk.fanout_size( n ) == ntk.value( n ) )
        ++lits;
    }
  } );

  return lits;
}

std::pair<double, double> abc_map( aig_network const& aig, std::string const& script )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; {}; stime;\"", script );

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
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( "ps", pos + 1 ) - pos - 1 );

    delay = std::stod( delay_res );
  }
  else
  {
    std::cout << "[e] failed to read the stime result\n";
  }

  return std::make_pair( area, delay );
}

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

#pragma region explore
template<class Lib>
void low_effort_optimization( aig_network& aig, int32_t opt_iterations, bool report_steps, Lib const& exact_lib, sop_factoring<aig_network>& sop_resyn )
{
  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 6u;
      sps.max_inserts = 1u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    {
      rewrite_params ps;
      rewrite_stats st;
      ps.use_mffc = false;
      ps.optimize_literal_cost = true;
      ps.allow_zero_gain = false;

      fanout_view fanout_aig{aig};
      rewrite( fanout_aig, exact_lib, ps, &st );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rw\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    // {
    //   fanout_view fanout_aig{aig};
    //   refactoring_params fps;
    //   fps.max_pis = 10;
    //   fps.allow_zero_gain = false;

    //   refactoring( aig, sop_resyn, fps );

    //   aig = cleanup_dangling( aig );
    // }

    // if ( report_steps )
    //   std::cout << fmt::format( "rf\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    if ( count_literals( aig ) >= lits_loop_before )
      break;

    aig = det_randomize( aig, opt_i );
  }
}

template<class Lib>
bool medium_effort_optimization( aig_network& aig, int32_t opt_i, bool report_steps, Lib const& exact_lib, sop_factoring<aig_network>& sop_resyn )
{
  const uint32_t lits_before = count_literals( aig );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 6u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    fanout_view fanout_aig{aig};
    refactoring_params fps;
    fps.max_pis = 10;
    fps.allow_zero_gain = false;

    refactoring( aig, sop_resyn, fps );

    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rf\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 8u;
    sps.max_inserts = 1u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  aig_balance( aig );

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 8u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    rewrite_params ps;
    rewrite_stats st;
    ps.use_mffc = false;
    ps.optimize_literal_cost = true;
    ps.allow_zero_gain = false;

    fanout_view fanout_aig{aig};
    rewrite( fanout_aig, exact_lib, ps, &st );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rw\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  aig = det_randomize( aig, opt_i );

  return count_literals( aig ) < lits_before;
}

template<class Lib>
void high_effort_optimization( aig_network& aig, int32_t opt_i, bool report_steps, Lib const& exact_lib, sop_factoring<aig_network>& sop_resyn, bool last = false )
{
  // aig_balance( aig );
  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 10u;
    sps.max_inserts = 1u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  if ( !last )
  {
    rewrite_params ps;
    rewrite_stats st;
    ps.use_mffc = false;
    ps.optimize_literal_cost = false;
    ps.allow_zero_gain = true;

    fanout_view fanout_aig{aig};
    rewrite( fanout_aig, exact_lib, ps, &st );
    aig = cleanup_dangling( aig );

    if ( report_steps )
      std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  }
  
  // {
  //   fanout_view fanout_aig{aig};
  //   refactoring_params fps;
  //   fps.max_pis = 8;
  //   fps.allow_zero_gain = false;

  //   refactoring( aig, sop_resyn, fps );

  //   aig = cleanup_dangling( aig );
  // }

  // if ( report_steps )
  //   std::cout << fmt::format( "rf\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 10u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  
  aig_balance( aig );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 12u;
    sps.max_inserts = 1u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  // {
  //   rewrite_params ps;
  //   rewrite_stats st;
  //   ps.use_mffc = false;
  //   ps.optimize_literal_cost = true;
  //   ps.allow_zero_gain = false;

  //   fanout_view fanout_aig{aig};
  //   rewrite( fanout_aig, exact_lib, ps, &st );
  //   aig = cleanup_dangling( aig );
  // }

  // if ( report_steps )
  //   std::cout << fmt::format( "rw\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    fanout_view fanout_aig{aig};
    refactoring_params fps;
    fps.max_pis = 10;
    fps.allow_zero_gain = true;

    refactoring( aig, sop_resyn, fps );

    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rfz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 12u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rs\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

  // if ( !last )
  // {
  //   rewrite_params ps;
  //   rewrite_stats st;
  //   ps.use_mffc = false;
  //   ps.optimize_literal_cost = false;
  //   ps.allow_zero_gain = true;

  //   fanout_view fanout_aig{aig};
  //   rewrite( fanout_aig, exact_lib, ps, &st );
  //   aig = cleanup_dangling( aig );

  //   if ( report_steps )
  //     std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  // }

  aig = det_randomize( aig, opt_i );
}

// void optimizer( aig_network& aig, uint32_t opt_iterations, bool report_steps )
// {
//   sop_factoring_params sop_ps;
//   sop_ps.use_boolean_division = false;
//   sop_ps.minimize_with_espresso = false;
//   sop_factoring<aig_network> sop_resyn( sop_ps );

//   xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
//   exact_library_params eps;
//   exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

//   int opt_i = opt_iterations;
//   while( opt_i-- > 0 )
//   {
//     const uint32_t lits_loop_before = count_literals( aig );

//     low_effort_optimization( aig, 10, report_steps, exact_lib, sop_resyn );

//     medium_effort_optimization( aig, opt_i, report_steps, exact_lib, sop_resyn );
    
//     high_effort_optimization( aig, opt_i, report_steps, exact_lib, sop_resyn, opt_i == 0 );

//     // if ( count_literals( aig ) >= lits_loop_before )
//     //   break;
//     // aig = det_randomize( aig, opt_i );
//   }

//   aig_balance( aig );

//   /* size recovery */
//   {
//     rewrite_params ps;
//     rewrite_stats st;
//     ps.use_mffc = false;
//     ps.optimize_literal_cost = true;
//     ps.allow_zero_gain = false;

//     fanout_view fanout_aig{aig};
//     rewrite( fanout_aig, exact_lib, ps, &st );
//     aig = cleanup_dangling( aig );
//   }

//   aig_balance( aig );

//   if ( report_steps )
//     std::cout << fmt::format( "rw; b\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
// }

void resub_opt( aig_network& aig, int32_t opt_iterations, bool report_steps )
{
  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    // aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 8u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rs -K 8 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    // aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    if ( count_literals( aig ) >= lits_loop_before )
      break;

    aig = det_randomize( aig, opt_i );
  }
}

template<class Lib>
void rewrite_opt( aig_network& aig, int32_t opt_iterations, bool report_steps, Lib const& exact_lib )
{
  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    // aig_balance( aig );

    {
      rewrite_params ps;
      rewrite_stats st;
      ps.use_mffc = false;
      ps.optimize_literal_cost = true;
      ps.allow_zero_gain = false;

      fanout_view fanout_aig{aig};
      rewrite( fanout_aig, exact_lib, ps, &st );
      aig = cleanup_dangling( aig );
    }

    if ( report_steps )
      std::cout << fmt::format( "rw\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    if ( count_literals( aig ) >= lits_loop_before )
      break;

    // aig = det_randomize( aig, opt_i );
  }
}

void refactor_opt( aig_network& aig, bool report_steps, sop_factoring<aig_network>& sop_resyn, bool zero_gain )
{
  {
    fanout_view fanout_aig{aig};
    refactoring_params fps;
    fps.max_pis = 10;
    fps.allow_zero_gain = zero_gain;

    refactoring( aig, sop_resyn, fps );

    aig = cleanup_dangling( aig );
  }

  if ( report_steps )
    std::cout << fmt::format( "rfz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
}

void optimizer2( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  low_effort_optimization( aig, 100, true, exact_lib, sop_resyn );

  aig_balance( aig );

  resub_opt( aig, 1, report_steps );

  rewrite_opt( aig, 3, report_steps, exact_lib );

  refactor_opt( aig, report_steps, sop_resyn, true );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 10u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );

    if ( report_steps )
      std::cout << fmt::format( "rs -K 2 -N 10\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  }

  // {
  //   rewrite_params ps;
  //   rewrite_stats st;
  //   ps.use_mffc = false;
  //   ps.optimize_literal_cost = false;
  //   ps.allow_zero_gain = true;

  //   fanout_view fanout_aig{aig};
  //   rewrite( fanout_aig, exact_lib, ps, &st );
  //   aig = cleanup_dangling( aig );

  //   if ( report_steps )
  //     std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  // }

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 3, report_steps, exact_lib );

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 2 -N 10\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    /* zero rewrite */
    // if ( opt_i > 1 )
    // {
    //   rewrite_params ps;
    //   rewrite_stats st;
    //   ps.use_mffc = false;
    //   ps.optimize_literal_cost = false;
    //   ps.allow_zero_gain = true;

    //   fanout_view fanout_aig{aig};
    //   rewrite( fanout_aig, exact_lib, ps, &st );
    //   aig = cleanup_dangling( aig );

    //   if ( report_steps )
    //     std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    // }

    // aig_balance( aig );
    // if ( count_literals( aig ) >= lits_loop_before )
    //   break;
    // aig = det_randomize( aig, opt_i );
  }

  /* size recovery */
  {
    rewrite_params ps;
    rewrite_stats st;
    ps.use_mffc = false;
    ps.optimize_literal_cost = true;
    ps.allow_zero_gain = false;

    fanout_view fanout_aig{aig};
    rewrite( fanout_aig, exact_lib, ps, &st );
    aig = cleanup_dangling( aig );
  }

  aig_balance( aig );

  if ( report_steps )
    std::cout << fmt::format( "rw; b\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
}

void optimizer3( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  aig_balance( aig );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 6u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );

    if ( report_steps )
      std::cout << fmt::format( "rs -K 6 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  }

  rewrite_opt( aig, 1, report_steps, exact_lib );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 10u;
    sps.max_inserts = 3u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );

    if ( report_steps )
      std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  }

  aig_balance( aig );

  {
    resubstitution_params sps;
    resubstitution_stats sst;

    sps.max_pis = 12u;
    sps.max_inserts = 2u;
    sps.progress = false;
    factor_resubstitution( aig, sps, &sst );
    aig = cleanup_dangling( aig );

    if ( report_steps )
      std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  }

  rewrite_opt( aig, 1, report_steps, exact_lib );

  refactor_opt( aig, report_steps, sop_resyn, true );

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 8u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 8 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 1, report_steps, exact_lib );
  }
}

void optimizer4( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  // aig_balance( aig );

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 6u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 6 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 1, report_steps, exact_lib );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    // rewrite_opt( aig, 1, report_steps, exact_lib );

    refactor_opt( aig, report_steps, sop_resyn, true );

    const uint32_t lits_loop_before = count_literals( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 8u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 8 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    // aig_balance( aig );

    // {
    //   resubstitution_params sps;
    //   resubstitution_stats sst;

    //   sps.max_pis = 10u;
    //   sps.max_inserts = 3u;
    //   sps.progress = false;
    //   factor_resubstitution( aig, sps, &sst );
    //   aig = cleanup_dangling( aig );

    //   if ( report_steps )
    //     std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    // }

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 1, report_steps, exact_lib );

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    {
      rewrite_params ps;
      rewrite_stats st;
      ps.use_mffc = false;
      ps.optimize_literal_cost = false;
      ps.allow_zero_gain = true;

      fanout_view fanout_aig{aig};
      rewrite( fanout_aig, exact_lib, ps, &st );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 1, report_steps, exact_lib );

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    aig_balance( aig );

    // aig = det_randomize( aig, opt_i );
  }
}

void optimizer5( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  // aig_balance( aig );

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    return;

    refactor_opt( aig, report_steps, sop_resyn, true );

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 10u;
      sps.max_inserts = 3u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 10 -N 3\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    rewrite_opt( aig, 1, report_steps, exact_lib );

    aig_balance( aig );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    {
      rewrite_params ps;
      rewrite_stats st;
      ps.use_mffc = false;
      ps.optimize_literal_cost = false;
      ps.allow_zero_gain = true;

      fanout_view fanout_aig{aig};
      rewrite( fanout_aig, exact_lib, ps, &st );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    refactor_opt( aig, report_steps, sop_resyn, true );

    if ( opt_i == 0 && opt_iterations != 1 )
      break;

    {
      rewrite_params ps;
      rewrite_stats st;
      ps.use_mffc = false;
      ps.optimize_literal_cost = false;
      ps.allow_zero_gain = true;

      fanout_view fanout_aig{aig};
      rewrite( fanout_aig, exact_lib, ps, &st );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rwz\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }

    refactor_opt( aig, report_steps, sop_resyn, true );

    {
      resubstitution_params sps;
      resubstitution_stats sst;

      sps.max_pis = 12u;
      sps.max_inserts = 2u;
      sps.progress = false;
      factor_resubstitution( aig, sps, &sst );
      aig = cleanup_dangling( aig );

      if ( report_steps )
        std::cout << fmt::format( "rs -K 12 -N 2\t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
    }
  }
}
#pragma endregion

#pragma region lastest code
void resub_opt( aig_network& aig, uint32_t k, uint32_t n, bool report_steps )
{
  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_pis = k;
  ps.max_inserts = n;
  ps.progress = false;
  factor_resubstitution( aig, ps, &st );
  aig = cleanup_dangling( aig );

  if ( report_steps )
    std::cout << fmt::format( "rs -K {} -N {}\t gates = {};\t lits = {}\n", k, n, aig.num_gates(), count_literals( aig ) );
}

template<class Lib>
void rewrite_opt( aig_network& aig, Lib const& exact_lib, bool lit_cost, bool zero_gain, bool report_steps )
{
  rewrite_params ps;
  rewrite_stats st;
  ps.use_mffc = false;
  ps.optimize_literal_cost = lit_cost;
  ps.allow_zero_gain = zero_gain;

  fanout_view fanout_aig{aig};
  rewrite( fanout_aig, exact_lib, ps, &st );
  aig = cleanup_dangling( aig );

  if ( report_steps )
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {}\n", zero_gain ? "rwz" : "rw", aig.num_gates(), count_literals( aig ) );
}

void refactor_opt( aig_network& aig, sop_factoring<aig_network>& sop_resyn, bool zero_gain, bool report_steps )
{
  fanout_view fanout_aig{aig};
  refactoring_params fps;
  fps.max_pis = 10;
  fps.allow_zero_gain = zero_gain;

  refactoring( aig, sop_resyn, fps );

  aig = cleanup_dangling( aig );

  if ( report_steps )
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {}\n", zero_gain ? "rfz" : "rf", aig.num_gates(), count_literals( aig ) );
}

void refactor_opt_new( aig_network& aig, sop_factoring<aig_network>& sop_resyn, bool zero_gain, bool report_steps )
{
  fanout_view fanout_aig{aig};
  refactoring_params fps;
  fps.max_pis = 10;
  fps.allow_zero_gain = zero_gain;

  refactoring( aig, sop_resyn, fps );

  aig = cleanup_dangling( aig );

  if ( report_steps )
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {}\n", zero_gain ? "rfz" : "rf", aig.num_gates(), count_literals( aig ) );
}

void optimizer_old( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );


  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );

    resub_opt( aig, 6, 1, report_steps);

    rewrite_opt( aig, exact_lib, true, false, report_steps );

    resub_opt( aig, 6, 2, report_steps);

    refactor_opt(aig, sop_resyn, false, report_steps );

    resub_opt( aig, 8, 1, report_steps);

    aig_balance( aig );

    resub_opt( aig, 8, 2, report_steps);

    rewrite_opt( aig, exact_lib, true, false, report_steps );

    resub_opt( aig, 10, 2, report_steps);

    rewrite_opt( aig, exact_lib, true, true, report_steps );

    resub_opt( aig, 10, 2, report_steps);

    aig_balance( aig );

    resub_opt( aig, 12, 1, report_steps);

    refactor_opt(aig, sop_resyn, true, report_steps );

    resub_opt( aig, 12, 2, report_steps);

    rewrite_opt( aig, exact_lib, true, true, report_steps );

    if ( count_literals( aig ) >= lits_loop_before )
      break;
  }

  /* size recovery */
  rewrite_opt( aig, exact_lib, true, false, report_steps );
}

void optimizer( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  uint32_t lits_loop_before = count_literals( aig );

  int opt_i = opt_iterations;
  while( opt_i-- > 0 )
  {
    resub_opt( aig, 12, 2, report_steps );

    refactor_opt_new( aig, sop_resyn, true, report_steps );

    refactor_opt_new( aig, sop_resyn, true, report_steps );

    resub_opt( aig, 10, 3, report_steps );

    rewrite_opt( aig, exact_lib, true, false, report_steps );

    aig_balance( aig );
    if ( report_steps )
        std::cout << fmt::format( "b          \t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );

    resub_opt( aig, 12, 2, report_steps );

    if ( opt_i == 0 || count_literals( aig ) >= lits_loop_before )
      break;
    
    lits_loop_before = count_literals( aig );

    rewrite_opt( aig, exact_lib, false, true, report_steps );

    resub_opt( aig, 12, 2, report_steps );

    refactor_opt_new( aig, sop_resyn, true, report_steps );

    resub_opt( aig, 12, 2, report_steps );

    refactor_opt_new( aig, sop_resyn, true, report_steps );

    rewrite_opt( aig, exact_lib, false, true, report_steps );

    refactor_opt_new( aig, sop_resyn, true, report_steps );
  }
}
#pragma endregion

int main( int argc, char **argv )
{
  using namespace experiments;

  if ( argc < 2 )
    return -1;

  std::string benchmark{ argv[1] };
  fmt::print( "[i] processing {}\n", benchmark );
  aig_network aig;
  if ( lorina::read_aiger( benchmark, aiger_reader( aig ) ) != lorina::return_code::success )
  {
    return -1;
  }

  const uint32_t size_before = aig.num_gates();
  const uint32_t depth_before = depth_view( aig ).depth();

  /* optimize */
  aig_network aig_abc = abc_opt( aig, "compress2rs" );
  // aig_network aig_abc = cleanup_dangling( aig );
  // aig = abc_opt( aig, "resyn2rs" );
  optimizer5( aig, 1, true );

  /* tech mapping */
  // auto map_abc = depth_before > 4000 ? abc_map( aig_abc, "read_lib -G 1 tcbn03e_bwph169l3p48cpd_base_lvttt_0p75v_25c_typical_ccs.lib; &get; &nf -R 1000; &put" ) : abc_map( aig_abc, "read_lib -G 1 tcbn03e_bwph169l3p48cpd_base_lvttt_0p75v_25c_typical_ccs.lib; amap" );
  // auto map_ff = depth_before > 4000 ? abc_map( aig, "read_lib -G 1 tcbn03e_bwph169l3p48cpd_base_lvttt_0p75v_25c_typical_ccs.lib; &get; &nf -R 1000; &put" ) : abc_map( aig, "read_lib -G 1 tcbn03e_bwph169l3p48cpd_base_lvttt_0p75v_25c_typical_ccs.lib; &get; amap" );
  std::pair<double, double> map_abc = {0, 0};
  std::pair<double, double> map_ff = {0, 0};

  /* report */
  fmt::print( "[i] ABC:\t gates = {}\t lits = {}\t area = {:>8.5f}\n[i] FFL:\t gates = {}\t lits = {}\t area = {:>8.5f}\n",
              aig_abc.num_gates(), count_literals( aig_abc ), map_abc.first, aig.num_gates(), count_literals( aig ), map_ff.first );

  write_aiger( aig, "res.aig" );


  return 0;
}
