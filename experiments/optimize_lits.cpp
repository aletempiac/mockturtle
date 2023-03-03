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
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/det_randomization.hpp>
#include <mockturtle/algorithms/factor_resub.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/write_aiger.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/properties/litcost.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <experiments.hpp>

using namespace mockturtle;

struct lut_literals_cost
{
  std::pair<uint32_t, uint32_t> operator()( uint32_t num_leaves ) const
  {
    if ( num_leaves < 2u )
      return {0u, 0u};
    return {num_leaves, 1u}; /* area, delay */
  }

  std::pair<uint32_t, uint32_t> operator()( kitty::dynamic_truth_table const& tt ) const
  {
    if ( tt.num_vars() < 2u )
      return {0u, 0u};
    return {factored_literal_cost( tt, false ), 1u}; /* area, delay */
  }
};

template<class Ntk>
uint32_t count_literals( Ntk& ntk )
{
  // ntk.clear_values();
  // ntk.foreach_po( [&]( auto const& f ) {
  //   ntk.incr_value( ntk.get_node( f ) );
  // } );

  // uint32_t lits = 0;
  // ntk.foreach_node( [&]( auto const& n ) {
  //   if ( ntk.is_constant( n ) )
  //   {
  //     return;
  //   }
  //   else if ( ntk.is_pi( n ) || ntk.fanout_size( n ) > 1 )
  //   {
  //     lits += ntk.fanout_size( n ) - ntk.value( n );
  //     if ( ntk.fanout_size( n ) == ntk.value( n ) )
  //       ++lits;
  //   }
  // } );
  uint32_t lits = 2 * ntk.num_gates() + ntk.num_pos();

  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.fanout_size( n ) == 1 )
      --lits;
  } );

  return lits;
}

std::pair<double, double> abc_map( aig_network const& aig, std::string const& script )
{
  write_aiger( aig, "/tmp/tmp.aig" );
  std::string command = fmt::format( "abc -q \"read /tmp/tmp.aig; {}; ps;\"", script );

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

  // std::size_t pos = result.find( "Area" );

  // if ( pos != std::string::npos )
  // {
  //   pos = result.find( "=", pos + 1 );
  //   std::string area_res = result.substr( pos + 1, result.find( "(", pos + 1 ) - pos - 1 );
  //   lorina::detail::trim( area_res );

  //   area = std::stod( area_res );

  //   pos = result.find( "=", pos + 1 );
  //   std::string delay_res = result.substr( pos + 1, result.find( "ps", pos + 1 ) - pos - 1 );

  //   delay = std::stod( delay_res );
  // }
  // else
  // {
  //   std::cout << "[e] failed to read the stime result\n";
  // }

  std::size_t pos = result.find( "area" );

  if ( pos != std::string::npos )
  {
    pos = result.find( "=", pos + 1 );
    std::string area_res = result.substr( pos + 1, result.find( "d", pos + 1 ) - pos - 1 );
    lorina::detail::trim( area_res );

    area = std::stod( area_res );

    pos = result.find( "=", pos + 1 );
    std::string delay_res = result.substr( pos + 1, result.find( "l", pos + 1 ) - pos - 1 );

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
void resub_opt( aig_network& aig, uint32_t k, uint32_t n, bool report_steps = false, bool depth_opt = false )
{
  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_pis = k;
  ps.max_inserts = n;
  ps.preserve_depth = depth_opt;
  ps.progress = false;
  factor_resubstitution( aig, ps, &st );
  aig = cleanup_dangling( aig );

  if ( report_steps )
    std::cout << fmt::format( "rs -K {} -N {}\t gates = {};\t lits = {};\t depth = {}\n", k, n, aig.num_gates(), count_literals( aig ), depth_view( aig ).depth() );
}

template<class Lib>
void rewrite_opt( aig_network& aig, Lib const& exact_lib, bool lit_cost, bool zero_gain, bool report_steps, bool depth_opt = false, bool aggressive = false )
{
  rewrite_params ps;
  rewrite_stats st;
  ps.use_mffc = false;
  ps.optimize_literal_cost = lit_cost;
  ps.allow_zero_gain = zero_gain;
  ps.preserve_depth = depth_opt;
  ps.aggressive_zero_gain = aggressive;

  aig_network aig_tmp = cleanup_dangling( aig );

  rewrite( aig_tmp, exact_lib, ps, &st );
  aig_tmp = cleanup_dangling( aig_tmp );

  // if ( depth_opt && depth_view( aig_tmp ).depth() <= depth_view( aig ).depth() ) /* TODO: remove workaround */
    aig = cleanup_dangling( aig_tmp );

  if ( report_steps )
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {};\t depth = {}\n", zero_gain ? "rwz" : "rw", aig.num_gates(), count_literals( aig ), depth_view( aig ).depth() );
}

void refactor_opt( aig_network& aig, sop_factoring<aig_network>& sop_resyn, bool zero_gain, bool report_steps, bool depth_opt = false )
{
  fanout_view fanout_aig{aig};
  refactoring_params fps;
  fps.max_pis = 10;
  fps.allow_zero_gain = zero_gain;
  fps.preserve_depth = depth_opt;

  aig_network aig_tmp = cleanup_dangling( aig );

  refactoring( aig_tmp, sop_resyn, fps );

  aig_tmp = cleanup_dangling( aig_tmp );

  if ( depth_opt && depth_view( aig_tmp ).depth() <= depth_view( aig ).depth() ) /* TODO: remove workaround */
  {
    aig = cleanup_dangling( aig_tmp );
  }
  else
  {
    aig = abc_opt( aig, "rfz" );
  }

  if ( report_steps )
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {};\t depth = {}\n", zero_gain ? "rfz" : "rf", aig.num_gates(), count_literals( aig ), depth_view( aig ).depth() );
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
    std::cout << fmt::format( "{}         \t gates = {};\t lits = {};\t depth = {}\n", zero_gain ? "rfz" : "rf", aig.num_gates(), count_literals( aig ), depth_view( aig ).depth() );
}

void optimizer_old_mid( aig_network& aig, uint32_t opt_iterations, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  int opt_i = opt_iterations;
  // refactor_opt(aig, sop_resyn, false, report_steps );

  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );
    const uint32_t size_loop_before = aig.num_gates();

    // aig_balance( aig );

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

    // if ( opt_i != 0 )
      rewrite_opt( aig, exact_lib, true, true, report_steps );

    aig_balance( aig );

    // if ( count_literals( aig ) >= lits_loop_before && aig.num_gates() >= size_loop_before )
    //   break;
  }

  /* size recovery */
  // rewrite_opt( aig, exact_lib, true, false, report_steps );
}

void optimizer_old( aig_network& aig, uint32_t opt_iterations, bool depth_opt, bool report_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete> resyn;
  exact_library_params eps;
  exact_library<aig_network, decltype( resyn )> exact_lib( resyn, eps );

  int opt_i = opt_iterations;
  // refactor_opt(aig, sop_resyn, false, report_steps );

  while( opt_i-- > 0 )
  {
    const uint32_t lits_loop_before = count_literals( aig );
    const uint32_t size_loop_before = aig.num_gates();

    // aig_balance( aig );

    resub_opt( aig, 6, 1, report_steps, depth_opt);

    rewrite_opt( aig, exact_lib, true, true, report_steps, depth_opt );

    resub_opt( aig, 6, 2, report_steps, depth_opt);

    refactor_opt(aig, sop_resyn, false, report_steps, depth_opt );

    resub_opt( aig, 8, 1, report_steps, depth_opt);

    aig_balance( aig );

    resub_opt( aig, 8, 2, report_steps, depth_opt);

    rewrite_opt( aig, exact_lib, true, false, report_steps, depth_opt );

    resub_opt( aig, 10, 1, report_steps, depth_opt);

    rewrite_opt( aig, exact_lib, true, true, report_steps, depth_opt, true );

    resub_opt( aig, 10, 2, report_steps, depth_opt);

    aig_balance( aig );

    resub_opt( aig, 12, 1, report_steps, depth_opt);

    refactor_opt(aig, sop_resyn, true, report_steps, depth_opt );

    resub_opt( aig, 12, 2, report_steps, depth_opt);

    // if ( opt_i != 0 )
    rewrite_opt( aig, exact_lib, true, true, report_steps, depth_opt );

    aig_balance( aig );

    // if ( count_literals( aig ) >= lits_loop_before && aig.num_gates() >= size_loop_before )
    //   break;
  }

  /* size recovery */
  // rewrite_opt( aig, exact_lib, true, false, report_steps );
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

  exact_library_params eps2;
  eps2.np_classification = true;
  exact_library<aig_network, decltype( resyn )> exact_lib2( resyn, eps2 );

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

aig_network remap( aig_network& aig )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  lut_map_params ps;
  ps.cut_enumeration_ps.cut_size = 5u;
  ps.cut_enumeration_ps.cut_limit = 25u;
  ps.recompute_cuts = false;
  ps.remove_dominated_cuts = false;
  ps.area_oriented_mapping = true;
  ps.cut_expansion = false;
  lut_map_stats st;
  mapping_view<aig_network, true> mapped_aig{aig};
  lut_map<decltype( mapped_aig ), true, lut_literals_cost>( mapped_aig, ps, &st );
  const auto klut = *collapse_mapped_network<klut_network>( mapped_aig );

  node_resynthesis_stats nst;
  aig_network res = node_resynthesis<aig_network>( klut, sop_resyn, {}, &nst );

  std::cout << fmt::format( "pre-map   \t gates = {};\t lits = {}\n", aig.num_gates(), count_literals( aig ) );
  std::cout << fmt::format( "remap     \t gates = {};\t lits = {}\n", res.num_gates(), count_literals( res ) );
  return res;
}

int main()
{
  using namespace experiments;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, double, double, double, double, double> exp( "ff_opt", "benchmark", "size_before", "size_abc", "size_ff","depth_before", "depth_abc", "depth_ff", "literals_before", "lits_abc", "lits_ff", "area_before", "area_abc", "area_ff", "delay_before", "delay_abc", "delay_ff" );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( ( benchmark == "leon2" || benchmark == "leon3" || benchmark == "leon3mp" || benchmark == "leon3_opt" || benchmark == "netcard" ) )
      continue;

    // if ( benchmark != "DMA" )
    //   continue;
    
    // if ( benchmark != "ac97_ctrl" )
    //   continue;

    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    // const uint32_t size_before = aig.num_gates();
    // const uint32_t depth_before = depth_view( aig ).depth();

    /* optimize */
    aig_network aig_abc = abc_opt( aig, "compress2rs; compress2rs" );
    // aig_network aig_abc = abc_opt( aig, "resyn2rs; resyn2rs" );
    // aig_network aig_abc = cleanup_dangling( aig );

    aig = cleanup_dangling( aig_abc );
    aig_abc = abc_opt( aig_abc, "compress2rs; compress2rs" );
    // aig_abc = abc_opt( aig_abc, "resyn2rs; resyn2rs" );
    // aig = abc_opt( aig, "resyn2rs" );

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();
    const uint32_t lits_before = count_literals( aig );
    // auto map_before_1 = depth_before > 4000 ? abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &dch -f; &nf -R 1000; &put;" ) : abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &dch -f; &nf -R 1000; &put;" );
    std::pair<double, double> map_before_1 = {2000000.0, 2000000.0};
    auto map_before_2 = depth_before > 4000 ? abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" ) : abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" );
    // std::pair<double, double> map_before_ff = {std::min( map_before_1.first, map_before_2.first), std::min( map_before_1.second, map_before_2.second)};
    double map_before_ff_delay, map_before_ff_area;

    if ( map_before_1.first < map_before_2.first || ( map_before_1.first == map_before_2.first && map_before_1.second < map_before_2.second ) )
    {
      map_before_ff_delay = map_before_1.second;
      map_before_ff_area = map_before_1.first;
    }
    else
    {
      map_before_ff_delay = map_before_2.second;
      map_before_ff_area = map_before_2.first;
    }

    // const uint32_t size_before = aig.num_gates();
    // const uint32_t depth_before = depth_view( aig ).depth();
    // const uint32_t lits_before = count_literals( aig );
    // // auto map_before_1 = abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &dch; &nf -C 25; &put" );
    // std::pair<double, double> map_before_1 = {2000000.0, 2000000.0};
    // auto map_before_2 = abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &nf -C 25; &put" );
    // std::pair<double, double> map_before_ff = {std::min( map_before_1.first, map_before_2.first), std::min( map_before_1.second, map_before_2.second)};
    // double map_before_ff_delay, map_before_ff_area;

    // if ( map_before_1.second < map_before_2.second || ( map_before_1.second == map_before_2.second && map_before_1.first < map_before_2.first ) )
    // {
    //   map_before_ff_delay = map_before_1.second;
    //   map_before_ff_area = map_before_1.first;
    // }
    // else
    // {
    //   map_before_ff_delay = map_before_2.second;
    //   map_before_ff_area = map_before_2.first;
    // }

    //  aig = remap( aig );
    // optimizer_old_mid( aig, 1, true );
    // aig = abc_opt( aig, "compress2rs" );
    optimizer_old( aig, 2, false, true );

    /* tech mapping */
    // auto map_abc_1 = depth_before > 4000 ? abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &dch -f; &nf -R 1000; &put;" ) : abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &dch -f; &nf -R 1000; &put;" );
    // auto map_ff_1 = depth_before > 4000 ? abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; dch -f; &get; &dch -f; &nf -R 1000; &put;" ) : abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &dch -f; &nf -R 1000; &put;" );
    std::pair<double, double> map_abc_1 = {20000000.0, 200000000.0};
    std::pair<double, double> map_ff_1 = {20000000.0, 200000000.0};
    auto map_abc_2 = depth_before > 4000 ? abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" ) : abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" );
    auto map_ff_2 = depth_before > 4000 ? abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" ) : abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b -l; &get; &nf -R 1000; &put;" );
    // std::pair<double, double> map_abc = {std::min( map_abc_1.first, map_abc_2.first), std::min( map_abc_1.second, map_abc_2.second)};
    // std::pair<double, double> map_ff = {std::min( map_ff_1.first, map_ff_2.first), std::min( map_ff_1.second, map_ff_2.second)};
    double map_abc_delay, map_abc_area;
    double map_ff_delay, map_ff_area;

    if ( map_abc_1.first < map_abc_2.first || ( map_abc_1.first == map_abc_2.first && map_abc_1.second < map_abc_2.second ) )
    {
      map_abc_delay = map_abc_1.second;
      map_abc_area = map_abc_1.first;
    }
    else
    {
      map_abc_delay = map_abc_2.second;
      map_abc_area = map_abc_2.first;
    }

    if ( map_ff_1.first < map_ff_2.first || ( map_ff_1.first == map_ff_2.first && map_ff_1.second < map_ff_2.second ) )
    {
      map_ff_delay = map_ff_1.second;
      map_ff_area = map_ff_1.first;
    }
    else
    {
      map_ff_delay = map_ff_2.second;
      map_ff_area = map_ff_2.first;
    }

    // auto map_abc_1 = abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &dch; &nf -C 25; &put" );
    // auto map_ff_1 = abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &dch; &nf -C 25; &put" );
    // std::pair<double, double> map_abc_1 = {20000000.0, 200000000.0};
    // std::pair<double, double> map_ff_1 = {20000000.0, 200000000.0};
    // auto map_abc_2 = abc_map( aig_abc, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &nf -C 25; &put" );
    // auto map_ff_2 = abc_map( aig, "read_lib -G 1 ../../../asap7_lib/asap7_clean.lib; b; &get; &nf -C 25; &put" );
    // double map_abc_delay, map_abc_area;
    // double map_ff_delay, map_ff_area;

    // if ( map_abc_1.second < map_abc_2.second || ( map_abc_1.second == map_abc_2.second && map_abc_1.first < map_abc_2.first ) )
    // {
    //   map_abc_delay = map_abc_1.second;
    //   map_abc_area = map_abc_1.first;
    // }
    // else
    // {
    //   map_abc_delay = map_abc_2.second;
    //   map_abc_area = map_abc_2.first;
    // }

    // if ( map_ff_1.second < map_ff_2.second || ( map_ff_1.second == map_ff_2.second && map_ff_1.first < map_ff_2.first ) )
    // {
    //   map_ff_delay = map_ff_1.second;
    //   map_ff_area = map_ff_1.first;
    // }
    // else
    // {
    //   map_ff_delay = map_ff_2.second;
    //   map_ff_area = map_ff_2.first;
    // }

    /* report */
    fmt::print( "[i] ABC:\t gates = {}\t lits = {}\t area = {:>8.5f}\t delay = {:>8.5f}\n[i] FFL:\t gates = {}\t lits = {}\t area = {:>8.5f}\t delay = {:>8.5f}\n",
                aig_abc.num_gates(), count_literals( aig_abc ), map_abc_area, map_abc_delay, aig.num_gates(), count_literals( aig ), map_ff_area, map_ff_delay );


    // const auto cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );

    exp( benchmark, size_before, aig_abc.num_gates(), aig.num_gates(), depth_before, depth_view( aig_abc ).depth(), depth_view( aig ).depth(), lits_before, count_literals( aig_abc ), count_literals( aig ), map_before_ff_area, map_abc_area, map_ff_area, map_before_ff_delay, map_abc_delay, map_ff_delay );
  }

  exp.save();
  exp.table();

  return 0;
}
