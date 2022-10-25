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
#include <thread>
#include <mutex>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/algorithms/aig_balancing.hpp>
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
using aig_resyn = xag_npn_resynthesis<aig_network, aig_network, xag_npn_db_kind::aig_complete>;

static constexpr uint32_t steps = 5;
static constexpr uint32_t num_benchmarks = 3;

#pragma region mutex
std::atomic<uint32_t> move_id{0};
#pragma endregion

#pragma shared memory
std::array<double, 1u << ( 3 * steps )> flow_reward;
#pragma endregion

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

void resub_opt( aig_network& aig, uint32_t k, uint32_t n )
{
  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_pis = k;
  ps.max_inserts = n;
  ps.progress = false;
  factor_resubstitution( aig, ps, &st );
  aig = cleanup_dangling( aig );
}

template<class Lib>
void rewrite_opt( aig_network& aig, Lib const& exact_lib, uint32_t zero_gain, bool optimize_literals )
{
  rewrite_params ps;
  rewrite_stats st;
  ps.use_mffc = false;
  ps.optimize_literal_cost = optimize_literals;
  ps.allow_zero_gain = zero_gain;

  fanout_view fanout_aig{aig};
  rewrite( fanout_aig, exact_lib, ps, &st );
  aig = cleanup_dangling( aig );
}

void refactor_opt( aig_network& aig, sop_factoring<aig_network>& sop_resyn, uint32_t zero_gain )
{
  fanout_view fanout_aig{aig};
  refactoring_params fps;
  fps.max_pis = 10;
  fps.allow_zero_gain = zero_gain;

  refactoring( aig, sop_resyn, fps );

  aig = cleanup_dangling( aig );
}

template<class Lib>
double execute( aig_network& aig, uint32_t move, Lib const& lib, sop_factoring<aig_network>& sop_resyn )
{
  double num_lits = count_literals( aig );
  double num_gates = aig.num_gates();

  if ( move == 0 )
  {
    aig_balance( aig );
  }
  else if ( move == 1 )
  {
    resub_opt( aig, 6, 2 );
  }
  else if ( move == 2 )
  {
    resub_opt( aig, 8, 2 );
  }
  else if ( move == 3 )
  {
    resub_opt( aig, 10, 3 );
  }
  else if ( move == 4 )
  {
    resub_opt( aig, 12, 2 );
  }
  else if ( move == 5 )
  {
    rewrite_opt( aig, lib, false, true );
  }
  else if ( move == 6 )
  {
    rewrite_opt( aig, lib, true, false );
  }
  else if ( move == 7 )
  {
    refactor_opt( aig, sop_resyn, true );
  }

  double reward = 0.9 * ( num_lits - static_cast<double>( count_literals( aig ) ) ) / num_lits + 0.1 * ( num_gates - static_cast<double>( aig.num_gates() ) ) / num_gates;

  return reward;
}

void goto_starting_point( std::vector<aig_network>& init_nets, uint64_t const starting_point, uint32_t starting_steps )
{
  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  aig_resyn resyn;
  exact_library_params eps;
  exact_library<aig_network, aig_resyn> exact_lib( resyn, eps );

  double reward = 0;

  for ( auto& aig : init_nets )
  {
    uint64_t moves = starting_point;
    uint32_t local_steps = starting_steps;

    while ( local_steps-- > 0 )
    {
      /* fetch instruction */
      uint32_t move = static_cast<uint32_t>( moves & 7 );
      moves = moves >> 3;

      /* execute */
      reward += execute( aig, move, exact_lib, sop_resyn );
    }
  }

  fmt::print( "[i] Starting point reached in {} steps with reward of {}\n", starting_steps, reward );
}

void thread_run( std::vector<aig_network> const& init_nets )
{
  std::vector<aig_network> nets( num_benchmarks );
  std::vector<std::pair<uint32_t, uint32_t>> cost( num_benchmarks );

  sop_factoring_params sop_ps;
  sop_ps.use_boolean_division = false;
  sop_ps.minimize_with_espresso = false;
  sop_factoring<aig_network> sop_resyn( sop_ps );

  aig_resyn resyn;
  exact_library_params eps;
  exact_library<aig_network, aig_resyn> exact_lib( resyn, eps );

  uint32_t id = move_id++;

  while ( id < std::pow( 8, steps ) )
  {
    /* clone initial designs */
    for ( uint32_t i = 0; i < num_benchmarks; ++i )
    {
      nets[i] = init_nets[i].clone();
      // rewards[i] = 0;
      cost[i] = std::make_pair( count_literals( nets[i] ), nets[i].num_gates() );
    }

    uint32_t local_steps = steps;
    uint32_t moves = id;
    bool skip = false;

    while ( local_steps-- > 0 )
    {
      /* fetch instruction */
      uint32_t move = moves & 7;
      moves = moves >> 3;

      /* execute */
      for ( uint32_t i = 0; i < num_benchmarks; ++i )
      {
        execute( nets[i], move, exact_lib, sop_resyn );
      }

      if ( skip )
        break;
    }

    if ( skip )
    {
      id = move_id++;
      continue;
    }

    /* compute total reward */
    double reward = 0;
    for ( uint32_t i = 0; i < num_benchmarks; ++i )
    {
      reward += 0.9 * ( cost[i].first - static_cast<double>( count_literals( nets[i] ) ) ) / cost[i].first + 0.1 * ( cost[i].second - static_cast<double>( nets[i].num_gates() ) ) / cost[i].second;
    }

    /* save reward */
    flow_reward[id] = reward;

    std::cout << fmt::format( "New reward {} %\n", reward / num_benchmarks * 100 );

    id = move_id++;
  }
}

int main( int argc, char **argv )
{
  using namespace experiments;

  std::vector<aig_network> initial_nets;
  initial_nets.reserve( num_benchmarks );

  bool use_starting_point = false;
  uint64_t starting_point = 0;

  if ( argc > 1 )
  {
    use_starting_point = true;
    starting_point = std::stoll( std::string( argv[1] ), nullptr, 0 );
  }

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    if ( benchmark == "hyp" || benchmark == "dec" || benchmark == "adder" || benchmark == "bar" )
      continue;

    fmt::print( "[i] adding {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    initial_nets.push_back( aig );
  }

  /* balance all the initial nets */
  // for ( auto& aig : initial_nets )
  //   aig_balance( aig );

  if ( use_starting_point )
  {
    goto_starting_point( initial_nets, starting_point, 10 );
  }

  move_id.store( 0 );
  std::vector<std::thread> threads;

  /* generate threads */
  const auto processor_count = std::thread::hardware_concurrency() - 4;

  fmt::print( "[i] Running on {} threads\n", processor_count );
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads.emplace_back( thread_run, initial_nets );
  }

  /* wait threads */
  for ( auto i = 0u; i < processor_count; ++i )
  {
    threads[i].join();
  }

  /* find the 10 best flows */
  std::array<std::pair<uint32_t, double>, 10> best_flows;
  for ( uint32_t i = 0; i < 10; ++i )
    best_flows[i] = std::make_pair( 0, 0 );

  double max10 = 0;
  for ( uint32_t i = 0; i < ( 1 << ( 3 * steps ) ); ++i )
  {
    std::cout << flow_reward[i] << "\n";

    if ( flow_reward[i] <= max10 )
      continue;

    for ( uint32_t j = 0; j < 10; ++j )
    {
      if ( flow_reward[i] > best_flows[j].second )
      {
        std::pair<uint32_t, double> tmp = best_flows[j];

        for ( uint32_t k = j + 1; k < 10; ++k )
        {
          std::swap( tmp, best_flows[k] );
        }

        best_flows[j] = std::make_pair( i, flow_reward[i] );

        break;
      }
    }

    max10 = best_flows[9].second;
  }

  for ( uint32_t i = 0; i < 10; ++i )
    std::cout << fmt::format( "{:2} : {:>10.8f}\t {}\n", i + 1, best_flows[i].second / num_benchmarks * 100, best_flows[i].first );

  return 0;
}
