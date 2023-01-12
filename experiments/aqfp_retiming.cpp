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

#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <lorina/verilog.hpp>

#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_dot.hpp>

#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/retime.hpp>
#include <mockturtle/algorithms/aqfp/buffer_insertion.hpp>
#include <mockturtle/algorithms/aqfp/aqfp_retiming.hpp>
#include <mockturtle/algorithms/aqfp/buffer_verification.hpp>
#include <mockturtle/algorithms/aqfp/aqfp_network_convertion.hpp>
#include <mockturtle/algorithms/aqfp/aqfp_depth_optimization.hpp>

#include <mockturtle/networks/aqfp.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/buffered.hpp>

#include <mockturtle/utils/tech_library.hpp>

#include <mockturtle/views/depth_view.hpp>

#include <experiments.hpp>

using namespace experiments;
using namespace mockturtle;

// template<class Ntk>
// class dot_drawer_mig() : default_dot_drawer
// {
//   dot_drawer_mig()

// };

template<typename T>
auto count_majorities( T& ntk )
{
  std::unordered_map<uint32_t, uint32_t> counts;
  ntk.foreach_gate( [&]( auto n ) { counts[ntk.fanin_size( n )]++; } );
  return counts;
}

template<class Ntk>
struct aqfp_depth_cost
{
  uint32_t operator()( Ntk const& ntk, node<Ntk> const& node ) const
  {
    if ( ntk.is_buf( node ) && ntk.fanout_size( node ) == 1 )
      return 0u;
    else
      return 1u;
  }
};

template<class Ntk>
struct aqfp_depth_lower_bound_cost
{
  uint32_t operator()( Ntk const& ntk, node<Ntk> const& node ) const
  {
    if ( ntk.is_buf( node ) )
    {
      if ( ntk.fanout_size( node ) > 1 )
      {
        uint32_t level = 0;
        ntk.foreach_fanin( node, [&]( auto const& f ) {
          if ( !ntk.is_buf( ntk.get_node( f ) ) )
            level = 1;
        } );
        return level;
      }
      else
      {
        return 0u;
      }
    }
    else
    {
        return 1u;
    }
  }
};

struct opt_params_t
{
  uint32_t optimization_rounds{ 1 };
  uint32_t max_remapping_rounds{ 1 };
  uint32_t max_resynthesis_rounds{ 10 };
  std::unordered_map<uint32_t, double> gate_costs{ { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters{ { 1u, 2.0 }, { 4u, 2.0 } };
  bool balance_pis{ true };
  bool branch_pis{ true };
  bool balance_pos{ true };
};

mig_network remapping_round( mig_network const& ntk, exact_library<mig_network, mig_npn_resynthesis>& exact_lib, opt_params_t const& opt_params )
{
  map_params psm;
  psm.skip_delay_round = false;
  map_stats stm;

  mig_network mig = cleanup_dangling( ntk );

  /* initial mig mapping, depth-oriented */
  for ( auto i = 0u; i < opt_params.max_remapping_rounds; ++i )
  {
    uint32_t old_mig_depth = depth_view( ntk ).depth();
    uint32_t old_mig_size = ntk.num_gates();

    mig_network mig_map = map( mig, exact_lib, psm, &stm );

    if ( depth_view( mig_map ).depth() > old_mig_depth ||
         ( depth_view( mig_map ).depth() == old_mig_depth && mig_map.num_gates() >= old_mig_size ) )
    {
      break;
    }
    mig = cleanup_dangling( mig_map );
  }

  return mig;
}

int main()
{
  opt_params_t opt_params;

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  exact_library<mig_network, mig_npn_resynthesis> exact_lib( resyn, eps );

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, uint32_t, double, uint32_t, uint32_t, uint32_t, double, bool> exp(
      "aqfp_retiming", "bench", "size_init", "Depth_init", "B/S_sched", "JJs_sched", "Depth_sched", "Time_sched (s)", "B/S_ret", "JJs_ret", "Depth_ret", "Time (s)", "cec" );

  uint32_t total_jjs = 0;
  uint32_t total_bufs = 0;

  double retiming_opt_ratio = 0;
  double num_benchmarks = 0;
  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    // if ( benchmark == "hyp" || benchmark == "div" || benchmark == "sqrt" )
    //   continue;
    // if ( benchmark == "hyp" )
    //   continue;

    // if ( benchmark != "c5315" )
    //   continue;

    mig_network mig;
    // if ( lorina::read_verilog( benchmark_aqfp_iscas_path( benchmark ), verilog_reader( mig ) ) != lorina::return_code::success )
    // {
    //   continue;
    // }

    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    /* main optimization loop */

    auto mig_opt = remapping_round( mig, exact_lib, opt_params );
    // auto mig_opt = cleanup_dangling( mig );

    const uint32_t size_before = mig_opt.num_gates();
    const uint32_t depth_before = depth_view( mig_opt ).depth();

    double total_runtime = 0;

    /* convert network to AQFP */
    aqfp_network aqfp = cleanup_dangling<mig_network, aqfp_network>( mig_opt );

    /* Buffer insertion params */
    buffer_insertion_params buf_ps;
    buf_ps.scheduling = buffer_insertion_params::depth_optimal;
    buf_ps.optimization_effort = buffer_insertion_params::none;
    buf_ps.max_chunk_size = 100;
    buf_ps.assume.splitter_capacity = 4u;
    buf_ps.assume.branch_pis = true;
    buf_ps.assume.balance_pis = true;
    buf_ps.assume.balance_pos = true;

    /* buffer insertion */
    buffer_insertion buf_inst( aqfp, buf_ps );
    buffered_aqfp_network buffered_aqfp;
    // unordered_node_map<std::pair<uint32_t, uint32_t>, buffered_aqfp_network> mobility( buffered_aqfp );
    uint32_t num_bufs = buf_inst.run( buffered_aqfp );
    uint32_t num_jjs =  aqfp.num_gates() * 6 + num_bufs * 2;
    uint32_t jj_depth = buf_inst.depth();
    const double sched_runtime = buf_inst.get_runtime();
    total_runtime += sched_runtime;

    aqfp_assumptions aqfp_ps;
    aqfp_ps.splitter_capacity = buf_ps.assume.splitter_capacity;
    aqfp_ps.branch_pis = buf_ps.assume.branch_pis;
    aqfp_ps.balance_pis = buf_ps.assume.balance_pis;
    aqfp_ps.balance_pos = buf_ps.assume.balance_pos;

    // depth_view buffered_aqfp_d{buffered_aqfp};
    // buffered_aqfp.foreach_gate( [&]( auto const& n ) {
    //   std::cout << fmt::format( "Node {} at level {}: \t min {}\t max {}\n", n, buffered_aqfp_d.level( n ), mobility[n].first, mobility[n].second );
    // } );

    // write_dot<buffered_mig_network, gate_dot_drawer<buffered_mig_network>>( buffered_mig, benchmark + ".dot" );

    /* retiming params */
    aqfp_retiming_params aps;
    aps.aqfp_assumptions_ps = aqfp_ps;
    aps.backwards_first = !is_ALAP;
    aps.iterations = 250;
    aps.verbose = true;
    aps.retime_splitters = true;

    /* chunk movement params */
    buffer_insertion_params buf_ps2 = buf_ps;
    buf_ps2.scheduling = buffer_insertion_params::provided;
    buf_ps2.optimization_effort = buffer_insertion_params::one_pass;

    double retiming_saved = buffered_aqfp.size();

    /* first retiming */
    {
      aqfp_retiming_stats ast;
      auto buf_aqfp_ret = aqfp_retiming( buffered_aqfp, aps, &ast );
      total_runtime += to_seconds( ast.time_total );
      buffered_aqfp = buf_aqfp_ret;
    }

    retiming_saved -= buffered_aqfp.size();

    aps.det_randomization = true;

    /* repeat loop */
    uint32_t iterations = 1;
    while ( iterations-- > 0 )
    {
      uint32_t size_previous = buffered_aqfp.size();

      /* chunk movement */
      aqfp_reconstruct_splitter_trees_params reconstruct_ps;
      reconstruct_ps.buffer_insertion_ps = buf_ps2;
      double runtime_restruct = 0;
      auto buf_aqfp_chunk = aqfp_reconstruct_splitter_trees( buffered_aqfp, reconstruct_ps, nullptr, &runtime_restruct );
      total_runtime += runtime_restruct;

      retiming_saved += buffered_aqfp.size();

      /* retiming */
      aqfp_retiming_stats ast;
      auto buf_aqfp_ret = aqfp_retiming( buf_aqfp_chunk, aps, &ast );
      total_runtime += to_seconds( ast.time_total );

      retiming_saved -= buffered_aqfp.size();

      if ( buf_aqfp_ret.size() >= size_previous )
        break;

      buffered_aqfp = buf_aqfp_ret;
    }

    // write_dot<buffered_aqfp_network, gate_dot_drawer<buffered_aqfp_network>>( buf_aqfp_ret, benchmark + "_ret.dot" );

    // /* depth optimization */
    // aqfp_optimize_depth_params adps;
    // adps.aqfp_assumptions_ps = aqfp_ps;
    // adps.verbose = false;
    // adps.iterations = 1;
    // auto buf_mig_final2_temp = aqfp_optimize_depth( buf_aqfp_ret, adps );
    // auto buf_mig_final2 = buf_aqfp_final;

    // aps.backward_only = true;
    // aps.retime_splitters = false;

    // // aps.iterations = 100;

    /* cec */
    // auto cec = abc_cec_aqfp( buffered_aqfp, benchmark );
    // auto cec = abc_cec( buffered_aqfp, benchmark );

    auto cec = verify_aqfp_buffer( buffered_aqfp, aqfp_ps );

    uint32_t num_jjs_ret = 0;
    uint32_t num_bufs_ret = 0;
    uint32_t jj_depth_ret = depth_view<buffered_aqfp_network>( buffered_aqfp ).depth();

    buffered_aqfp.foreach_node( [&]( auto const& n ) {
      if ( buffered_aqfp.is_pi( n ) || buffered_aqfp.is_constant( n ) )
        return;
      if ( buffered_aqfp.is_buf( n ) )
      {
        ++num_bufs_ret;
        num_jjs_ret += 2;
      }
      else
      {
        num_jjs_ret += 6;
      }
    } );

    total_bufs += num_bufs_ret;
    total_jjs += num_jjs_ret;

    if ( ( num_bufs - num_bufs_ret ) > 0)
    {
      retiming_opt_ratio += retiming_saved / static_cast<double>( num_bufs - num_bufs_ret );
      ++num_benchmarks;
    }

    exp( benchmark, size_before, depth_before, num_bufs, num_jjs, jj_depth, sched_runtime, num_bufs_ret, num_jjs_ret, jj_depth_ret, total_runtime, cec );
  }

  exp.save();
  exp.table();

  std::cout << fmt::format( "[i] Total B/S = {} \tTotal JJs = {}\n", total_bufs, total_jjs );
  std::cout << "Ratio: " << retiming_opt_ratio * 100.0 / num_benchmarks << "\n";

  return 0;
}
