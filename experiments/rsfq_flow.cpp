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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
#include <mockturtle/algorithms/xag_balancing.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/lut_mapper.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/xag_optimization.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_minmc.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_minmc2.hpp>
#include <mockturtle/algorithms/xag_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/xmg_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/retiming.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_path_balancing.hpp>
// #include <mockturtle/algorithms/rsfq/rsfq_splitter.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_network_conversion.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/networks/generic.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>

#include <experiments.hpp>

static constexpr uint32_t splitter_jj = 3;

void aig_prepare()
{
  using namespace experiments;
  using namespace mockturtle;

  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    std::string command = fmt::format( "abc -q \"&read {}; &fraig -x; &put; compress2rs; compress2rs; if -g; resyn2rs; write_aiger {}\"", benchmark_path( benchmark ), "rsfq_opt/" + benchmark + ".aig" );

    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
    if ( !pipe )
    {
      throw std::runtime_error( "popen() failed" );
    }
    while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
    {
      result += buffer.data();
    }
  }
}

template<class Ntk>
struct rsfq_cost
{
  uint32_t operator()( Ntk const& ntk, mockturtle::node<Ntk> const& node ) const
  {
    if ( ntk.is_and( node ) )
    {
      return 11u;
    }
    else
    {
      return 9u;
    }
  }
};

mockturtle::xag_network depth_opt( mockturtle::xag_network const& xag_start, bool xor_opt = false )
{
  using namespace mockturtle;

  /* exact XAG database */
  using xag_resyn = xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_incomplete>;
  xag_resyn resyn;
  exact_library_params ps;
  ps.np_classification = true;
  exact_library<xag_network, xag_resyn> exact_lib( resyn, ps );
  ps.np_classification = false;
  exact_library<xag_network, xag_resyn> rw_lib( resyn, ps );

  /* MC database */
  future::xag_minmc_resynthesis mc_resyn;

  xag_network xag = cleanup_dangling( xag_start );

  /*  XAG algebraic rewriting */
  {
    auto xag_rw = cleanup_dangling( xag );
    fanout_view xag_fout{xag_rw};
    depth_view d_xag{xag_fout};
    fmt::print("Pre RW XAG:      size = {}\t depth = {}\n", xag_rw.num_gates(), d_xag.depth() );
    xag_algebraic_depth_rewriting_params ps;
    ps.allow_area_increase = true;
    xag_algebraic_depth_rewriting( d_xag, ps );
    xag_rw = cleanup_dangling( xag_rw );

    if ( d_xag.depth() < depth_view( xag ).depth() )
      xag = cleanup_dangling( xag_rw );

    fmt::print("Post RW XAG:     size = {}\t depth = {}\n", xag.num_gates(), depth_view( xag ).depth() );
  }

  /* delay-oriented remapping */
  for ( auto i = 0; i < 5; ++i )
  {
    uint32_t old_xag_depth = depth_view( xag ).depth();
    uint32_t old_xag_size = xag.num_gates();

    auto xag_map = cleanup_dangling( xag );
    xag_balance( xag_map, { false } );
    xag_network new_xag = map( xag_map, exact_lib );
    
    if ( depth_view( new_xag ).depth() > old_xag_depth ||
            ( depth_view( new_xag ).depth() == old_xag_depth && new_xag.num_gates() >= old_xag_size ) )
    {
      break;
    }
    xag = cleanup_dangling( new_xag );
  }
  fmt::print("Map XAG:     size = {}\t depth = {}\n", xag.num_gates(), depth_view( xag ).depth() );

  /* recover area */
  {
    rewrite_params cps;
    cps.preserve_depth = true;
    cps.allow_zero_gain = true;
    for ( auto i = 0; i < 2; ++i )
    {
      uint32_t xag_gates_before = xag.num_gates();
      rewrite( xag, rw_lib, cps );
      xag = cleanup_dangling( xag );

      if ( xag.num_gates() >= xag_gates_before )
        break;
    }
    fmt::print("ARec RW XAG:     size = {}\t depth = {}\n", xag.num_gates(), depth_view( xag ).depth() );
  }

  /* XOR optimization: multiplicative complexity reduction */
  if ( xor_opt )
  {
    /* Do multiplicative complexity optimization */
    {
      uint32_t num_and = 0u, num_xor = 0u;
      xag.foreach_gate( [&]( auto f ) {
        if ( xag.is_and( f ) )
        {
          num_and++;
        }
        else if ( xag.is_xor( f ) )
        {
          num_xor++;
        }
      } );
      std::cout << "Pre  MC opt: num ANDs = " << num_and << " \t" << "num XORs = " << num_xor << "\n";

      cut_rewriting_params cps;
      cps.cut_enumeration_ps.cut_size = 5;
      cps.preserve_depth = true;
      xag = cut_rewriting<xag_network, future::xag_minmc_resynthesis<xag_network>, rsfq_cost<xag_network>>( xag, mc_resyn, cps );

      num_and = 0u, num_xor = 0u;
      xag.foreach_gate( [&]( auto f ) {
        if ( xag.is_and( f ) )
        {
          num_and++;
        }
        else if ( xag.is_xor( f ) )
        {
          num_xor++;
        }
      } );
      std::cout << "Post MC opt: num ANDs = " << num_and << " \t" << "num XORs = " << num_xor << "\n";
    }

    /* Don't cares based multiplicative complexity reduction */
    // {
    //   xag = xag_dont_cares_optimization_windowed( xag );

    //   uint32_t num_and = 0u, num_xor = 0u;
    //   xag.foreach_gate( [&]( auto f ) {
    //     if ( xag.is_and( f ) )
    //     {
    //       num_and++;
    //     }
    //     else if ( xag.is_xor( f ) )
    //     {
    //       num_xor++;
    //     }
    //   } );
    //   std::cout << "Post DC opt: num ANDs = " << num_and << " \t" << "num XORs = " << num_xor << "\n";
    // }
  }

  return xag;
}

void rsfq_flow( int opt_iter )
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, double, float, bool> exp(
      "RFSQ flow", "benchmark", "size", "depth", "size_opt", "depth_opt", "area", "delay", "runtime", "equivalent" );

  fmt::print( "[i] processing RSFQ technology library\n" );

  /* library to map to technology */
  std::vector<gate> gates;
  std::ifstream in( "/Users/tempia/Documents/phd/libraries/aletempiac/mockturtle/experiments/cell_libraries/suny_rsfq_cell_library.genlib" );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return;
  }

  super_lib super_data;
  std::ifstream in_super( "/Users/tempia/Documents/phd/libraries/aletempiac/mockturtle/experiments/cell_libraries/suny_rsfq_cell_library.super" );

  if ( lorina::read_super( in_super, super_reader( super_data ) ) != lorina::return_code::success )
  {
    return;
  }

  tech_library_params tps;
  tech_library<5, classification_type::np_configurations> tech_lib( gates, super_data, tps );
  // tech_library<5, classification_type::np_configurations> tech_lib( gates, tps );


  /* exact XAG database */
  using xag_resyn = xag_npn_resynthesis<xag_network, xag_network, xag_npn_db_kind::xag_incomplete>;
  xag_resyn resyn;
  exact_library<xag_network, xag_resyn> exact_lib( resyn );

  /* SOP balancing */
  sop_rebalancing<aig_network> sop_balancing;

  map_params mps;
  mps.enable_logic_sharing = true;
  mps.logic_sharing_cut_limit = 1;
  mps.skip_delay_round = true;

  generic_network net;

  /* flow */
  for ( auto const& benchmark : iscas_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    if ( benchmark == "hyp" || benchmark == "sqrt" )
      continue;

    xag_network aig;
    if ( lorina::read_aiger( "rsfq_opt/" + benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    // if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    // {
    //   continue;
    // }

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    fmt::print("Initial AIG: size = {}\t depth = {}\n", size_before, depth_before );

    xag_network fraig = cleanup_dangling( aig );

    /* functional reduction */
    functional_reduction( fraig );
    fraig = cleanup_dangling( fraig );

    /* Map to XAG */
    aig_balance( fraig, { false } );
    xag_network xag = map( fraig, exact_lib );

    /* optimization steps */
    for ( auto i = 0; i < opt_iter; ++i )
    {
      xag_network xag_opt = depth_opt( xag, false );

      if ( depth_view( xag_opt ).depth() > depth_view( xag ).depth() ||
          ( depth_view( xag_opt ).depth() == depth_view( xag ).depth() && xag_opt.num_gates() >= xag.num_gates() ) )
      {
        break;
      }

      xag = cleanup_dangling( xag_opt );
    }

    /* Technology mapping */
    map_params ps;
    ps.cut_enumeration_ps.minimize_truth_table = true;
    ps.cut_enumeration_ps.cut_limit = 49;
    map_stats st;
    xag_balance( xag, { true } );
    binding_view<klut_network> res = map( xag, tech_lib, ps, &st );

    // write_verilog( res, benchmark + ".v" );

    /* path balancing with buffers */
    auto res_test = rsfq_path_balancing( res );

    fmt::print("TMap XAG:    area = {:>5.2f}\t delay = {:>5.2f}\n", res_test.compute_area(), res_test.compute_worst_delay() );

    /* retime registers */
    retime_params rps;
    retime_stats rst;
    rps.verbose = false;
    auto net = rsfq_generic_network_create_from_mapped( res_test );
    retime( net, rps, &rst );
    auto retime_res = rsfq_mapped_create_from_generic_network( net );

    /* splitter insertion */
    double area_final = retime_res.compute_area();
    retime_res.foreach_node( [&]( auto const& n ) {
      if ( !retime_res.is_constant( n ) )
        area_final += splitter_jj * ( retime_res.fanout_size( n ) - 1 );
    } );

    // const auto cec = benchmark == "hyp" ? true : abc_cec( retime_res, benchmark );
    const auto cec = true;

    std::cout << "Area after retime and splitters: " << area_final << " check: " << rsfq_check_buffering( retime_res );
    std::cout << " cec: " << cec << "\n";

    exp( benchmark, size_before, depth_before, xag.num_gates(), depth_view( xag ).depth(), area_final, retime_res.compute_worst_delay(), to_seconds( rst.time_total ), cec );
  }

  exp.save();
  exp.table();
}


int main( int argc, char** argv)
{
  int opt_iter = 1;

  if ( argc > 1 )
  {
    opt_iter = atoi( argv[1] );
  }

  // aig_prepare();
  rsfq_flow( opt_iter );

  return 0;
}