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

#include "experiments_eleonora.hpp"

#include <mockturtle/algorithms/aqfp/mig_algebraic_rewriting_splitters.hpp>
#include <mockturtle/algorithms/aqfp/mig_resub_splitters.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/akers.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/aqfp_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_limit_view.hpp>

#include <fmt/format.h>
#include <lorina/lorina.hpp>

#include <future>
#include <mutex>
#include <string>
#include <vector>

template<class Ntk>
struct jj_cost
{
  uint32_t operator()( Ntk const& ntk, mockturtle::node<Ntk> const& n ) const
  {
    uint32_t cost = 0;
    if ( ntk.fanout_size( n ) == 1 )
      cost = ntk.fanout_size( n );
    else if ( ntk.fanout_size( n ) <= 4 )
      cost = 3;
    else
      cost = 11;
    return cost;
  }
};

template<class Ntk>
struct fanout_cost_depth_local
{
  uint32_t operator()( Ntk const& ntk, mockturtle::node<Ntk> const& n ) const
  {
    uint32_t cost = 0;
    if ( ntk.is_pi( n ) )
      cost = 0;
    else if ( ntk.fanout_size( n ) == 1 )
      cost = 1;
    else if ( ( ntk.fanout_size( n ) > 1 ) && ( ( ntk.fanout_size( n ) < 5 ) ) )
      cost = 2;
    else if ( ( ( ntk.fanout_size( n ) > 4 ) ) )
      cost = 3;
    return cost;
  }
};

using namespace mockturtle;

using cost_fn_t = fanout_cost_depth_local<mig_network>;
using limit_view_t = fanout_limit_view<mig_network>;
using aqfp_view_t = aqfp_view<limit_view_t, true>;
using depth_view_t = depth_view<limit_view_t>;
using jj_depth_view_t = depth_view<limit_view_t, cost_fn_t>;

void get_statistics( mig_network const& mig, uint32_t& size, uint32_t& depth, uint32_t& jj_count, uint32_t& jj_depth )
{
  limit_view_t mig_limited = cleanup_dangling<mig_network, limit_view_t>( mig );
  aqfp_view_t mig_aqfp( mig_limited );
  depth_view_t mig_depth( mig_limited );
  jj_depth_view_t mig_jj_depth( mig_limited, cost_fn_t() );

  size = mig_limited.num_gates();
  depth = mig_depth.depth();
  jj_count = size * 6 + mig_aqfp.num_buffers() * 2;
  jj_depth = mig_jj_depth.depth();
}

template<typename Fn1, typename Fn2>
void do_experiment( std::string benchmark_path, Fn1& callback1, Fn2& callback2 )
{
  auto last_dir_delim_ind = benchmark_path.find_last_of( "/" );
  std::string benchmark_name = ( last_dir_delim_ind == std::string::npos ) ? benchmark_path : benchmark_path.substr( last_dir_delim_ind + 1 );

  auto last_ext_delim_ind = benchmark_name.find_last_of( "." );
  assert( last_ext_delim_ind != std::string::npos );

  bool is_verilog = ( benchmark_name.substr( last_ext_delim_ind + 1 ) == "v" );
  if ( !is_verilog )
  {
    assert( benchmark_name.substr( last_ext_delim_ind + 1 ) == "aig" );
  }
  benchmark_name = benchmark_name.substr( 0, last_ext_delim_ind );

  std::cerr << "reading the benchmark\n";

  mockturtle::mig_network mig;
  if ( is_verilog )
  {
    lorina::read_verilog( benchmark_path, mockturtle::verilog_reader( mig ) );
  }
  else
  {
    lorina::read_aiger( benchmark_path, mockturtle::aiger_reader( mig ) );
  }

  std::cerr << "reading the benchmark done\n";

  uint32_t size_before, depth_before, jj_before, jj_levels_before;
  {
    get_statistics( mig, size_before, depth_before, jj_before, jj_levels_before );
  }

  return;
  auto t1 = std::chrono::high_resolution_clock::now();

  auto iteration = 0u;
  std::cerr << fmt::format( "benchmark {} starting point: size = {}, depth = {}, JJ count = {}, JJ depth = {}\n", benchmark_name, size_before, depth_before, jj_before, jj_levels_before );

  auto iter = 0;
  while ( true )
  {
    uint32_t size, depth, jj_count, jj_depth;
    get_statistics( mig, size, depth, jj_count, jj_depth );
    auto const jj_depth_before_rewrite = jj_depth;

    ++iteration;
    std::cerr << fmt::format( "benchmark {} iteration {}: size = {}, JJ depth = {}\n", benchmark_name, iteration, size, jj_depth );

    /* Section 3.2: Depth optimization with algebraic rewriting -- limiting fanout size increase */
    {
      std::cerr << fmt::format( "benchmark {} starting algebraic rewriting\n", benchmark_name );

      mig_algebraic_depth_rewriting_params ps_alg_rewrite;
      ps_alg_rewrite.overhead = 1.5;
      ps_alg_rewrite.strategy = mig_algebraic_depth_rewriting_params::dfs;
      ps_alg_rewrite.allow_area_increase = true;

      limit_view_t mig_limited = cleanup_dangling<mig_network, limit_view_t>( mig );
      jj_depth_view_t mig_jj_depth{ mig_limited, cost_fn_t() };
      mig_algebraic_depth_rewriting_splitters( mig_jj_depth, ps_alg_rewrite );
      mig = cleanup_dangling<jj_depth_view_t, mig_network>( mig_jj_depth );

      // std::cerr << fmt::format( "benchmark {} done algebraic rewriting\n", benchmark_name );
    }

    get_statistics( mig, size, depth, jj_count, jj_depth );
    auto const jj_depth_after_rewrite = jj_depth;
    auto const size_before_resub = size;

    /* Section 3.3: Size optimization with Boolean resubstitution -- considering fanout size limitation */
    {
      std::cerr << fmt::format( "benchmark {} starting resubstitution\n", benchmark_name );

      resubstitution_params ps_resub;
      ps_resub.max_divisors = 250;
      ps_resub.max_inserts = 1;
      ps_resub.preserve_depth = true;

      limit_view_t mig_limited = cleanup_dangling<mig_network, limit_view_t>( mig );
      jj_depth_view_t mig_jj_depth{ mig_limited, cost_fn_t() };
      mig_resubstitution_splitters( mig_jj_depth, ps_resub );
      mig = cleanup_dangling<jj_depth_view_t, mig_network>( mig_jj_depth );

      //  std::cerr << fmt::format( "benchmark {} done resubstitution\n", benchmark_name );
    }

    get_statistics( mig, size, depth, jj_count, jj_depth );

    /* Section 3.4: Further size optimization with refactoring */

    mig_network mig_copy = mig;
    auto const size_before_refactor = size;
    auto const depth_before_refactor = depth;
    auto const jj_depth_before_refactor = jj_depth;

    {
      std::cerr << fmt::format( "benchmark {} starting akers synthesis\n", benchmark_name );

      limit_view_t mig_limited = cleanup_dangling<mig_network, limit_view_t>( mig );
      akers_resynthesis<mig_network> resyn;
      refactoring( mig_limited, resyn, {}, nullptr, jj_cost<mig_network>() );
      mig = cleanup_dangling<limit_view_t, mig_network>( mig_limited );

            std::cerr << fmt::format( "benchmark {} done akers synthesis\n", benchmark_name );
    }

    get_statistics( mig, size, depth, jj_count, jj_depth );

    /* Undo refactoring if (1) size increases; or (2) JJ depth increases; or (3) depth increases */
    if ( ( size > size_before_refactor ) || ( jj_depth > jj_depth_before_refactor ) || ( depth > depth_before_refactor ) )
    {
      mig = mig_copy;
    }

    {
      get_statistics( mig, size, depth, jj_count, jj_depth );
    }

    /* Terminate when (1) [resub + refactor] cannot decrease size anymore; or (2) rewriting cannot decrease JJ depth anymore */
    if ( ( size >= size_before_resub ) || ( jj_depth_after_rewrite >= jj_depth_before_rewrite ) )
    {
      break;
    }
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  bool const cec = true; // experiments::abc_cec_impl( mig, benchmark_path );

  auto t3 = std::chrono::high_resolution_clock::now();

  auto exp_time = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  auto ver_time = std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 ).count();

  write_verilog( mig, benchmark_name + "_after_eleonora.v" );

  uint32_t size_after, depth_after, jj_after, jj_levels_after;
  {
    get_statistics( mig, size_after, depth_after, jj_after, jj_levels_after );
  }


  std::cerr << fmt::format( "benchmark {} after AQFP flow: size = {}, depth = {}, JJ count = {}, JJ depth = {}\n", benchmark_name, size_after, depth_after, jj_after, jj_levels_after );

  float const impr_size = ( (float)size_before - (float)size_after ) / (float)size_before * 100;
  float const impr_depth = ( (float)depth_before - (float)depth_after ) / (float)depth_before * 100;
  float const impr_jj = ( (float)jj_before - (float)jj_after ) / (float)jj_before * 100;
  float const impr_levels = ( (float)jj_levels_before - (float)jj_levels_after ) / (float)jj_levels_before * 100;

  callback1( benchmark_name, size_before, size_after, impr_size, depth_before, depth_after, impr_depth, cec );
  callback2( benchmark_name, jj_before, jj_after, impr_jj, jj_levels_before, jj_levels_after, impr_levels, exp_time / 1000.0, ver_time / 1000.0, cec );
}

int main( int argc, char** argv )
{
  using namespace experiments;

  bool const verbose = true; // turn on/off verbose printing of intermediate results
  int iteration;

  experiment<std::string, uint32_t, uint32_t, float, uint32_t, uint32_t, float, bool>
      exp1( "mcnc_table1", "benchmark", "size MIG", "Size Opt MIG", "Impr. Size", "depth MIG", "depth Opt MIG", "Impr. depth", "eq cec" );
  experiment<std::string, uint32_t, uint32_t, float, uint32_t, uint32_t, float, double, double, bool>
      exp2( "mcnc_table3", "benchmark", "jj MIG", "jj Opt MIG", "Impr. jj", "jj levels MIG", "jj levels Opt MIG", "Impr. jj levels", "exp time", "verif time", "eq cec" );

  std::vector<std::string> benchmark_paths;

  if ( argc > 1 )
  {
    std::string b( argv[1] );
    bool is_verilog = std::string( argv[4] ) == "verilog";
    if ( is_verilog )
    {
      benchmark_paths.push_back( fmt::format( "./benchmarks/{}.v", b ) );
    }
    else
    {
      benchmark_paths.push_back( experiments::benchmark_path( b ) );
    }
  }
  else
  {
    for ( auto benchmark : epfl_benchmarks() )
    {
      benchmark_paths.push_back( experiments::benchmark_path( benchmark ) );
      if ( benchmark_paths.size() == 3 )
        break;
    }
  }

  std::mutex mu;

  std::vector<std::thread> res;
    auto cb1 = [&]( std::string a, uint32_t b, uint32_t c, float d, uint32_t e, uint32_t f, float g, float h ) {
      fmt::print( "Table 1" );
      exp1( a, b, c, d, e, f, g, h );
      exp1.save();
      exp1.table();
    };

    auto cb2 = [&]( std::string a, uint32_t b, uint32_t c, float d, uint32_t e, uint32_t f, float g, double h, double i, bool j ) {
      fmt::print( "Table 2" );
      exp2( a, b, c, d, e, f, g, h, i, j );
      exp2.save();
      exp2.table();
    };

  for ( auto const& benchmark_path : benchmark_paths )
  {
    do_experiment(benchmark_path, cb1, cb2);
  }

  fmt::print( "Table 1: Results for size and depth optimization over MIG\n" );
  exp1.save();
  exp1.table();

  fmt::print( "Table 3: Results for area, delay, and number of buffers & splitters for MIGs mapped into AQFP technology\n" );
  exp2.save();
  exp2.table();

  return 0;
}
