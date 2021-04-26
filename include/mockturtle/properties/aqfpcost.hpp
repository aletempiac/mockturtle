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

/*!
  \file aqfpcost.hpp
  \brief Cost functions for AQFP networks

  \author Dewmini Marakkalage 
*/

#pragma once

#include <algorithm>
#include <limits>

#include "../utils/hash_functions.hpp"
#include "../views/fanout_view.hpp"

namespace mockturtle
{

/*! \brief Cost function for computing the best splitter and buffer cost for a fanout net with given relative levels. */
class balanced_fanout_net_cost
{
public:
  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

  balanced_fanout_net_cost( const std::unordered_map<uint32_t, double>& splitters )
      : buffer_cost( splitters.at( 1u ) ), splitters( remove_buffer( splitters ) )
  {
  }

  double operator()( const std::vector<uint32_t>& config )
  {
    return cost_for_config( config, false );
  }

  double operator()( const std::vector<uint32_t>& config, bool ignore_initial_buffers )
  {
    return cost_for_config( config, ignore_initial_buffers );
  }

private:
  double buffer_cost;
  std::unordered_map<uint32_t, double> splitters;
  std::unordered_map<std::tuple<bool, std::vector<uint32_t>>, double, hash<std::tuple<bool, std::vector<uint32_t>>>> cache;

  static std::unordered_map<uint32_t, double> remove_buffer( std::unordered_map<uint32_t, double> splitters )
  {
    splitters.erase( 1u );
    return splitters;
  }

  double cost_for_config( const std::vector<uint32_t> config, bool ignore_initial_buffers )
  {
    if ( config.size() == 1 )
    {
      if ( config[0] >= 1 )
      {
        return ignore_initial_buffers ? 0.0 : ( config[0] - 1 ) * buffer_cost;
      }
      else
      {
        return IMPOSSIBLE;
      }
    }

    std::tuple key = { ignore_initial_buffers, config };
    if ( cache.count( key ) )
    {
      return cache[key];
    }

    auto result = IMPOSSIBLE;

    for ( const auto& s : splitters )
    {
      for ( auto size = 2u; size <= std::min( s.first, uint32_t( config.size() ) ); size++ )
      {
        auto sp_lev = config[config.size() - size] - 1;
        if ( sp_lev == 0 )
        {
          continue;
        }

        auto temp = s.second;

        for ( auto i = config.size() - size; i < config.size(); i++ )
        {
          temp += ( config[i] - config[config.size() - size] ) * buffer_cost;
        }

        std::vector<uint32_t> new_config( config.begin(), config.begin() + ( config.size() - size ) );
        new_config.push_back( sp_lev );
        std::sort( new_config.begin(), new_config.end() );

        temp += cost_for_config( new_config, ignore_initial_buffers );

        if ( temp < result )
        {
          result = temp;
        }
      }
    }

    return ( cache[key] = result );
  }
};

/*! \brief Cost function for computing the cost of a path-balanced AQFP network with a given assignment of node levels. 
  *
  * Assumes no path balancing or splitters are needed for primary inputs or register outputs. 
  */
struct aqfp_network_cost
{
  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

  aqfp_network_cost( const std::unordered_map<uint32_t, double>& gate_costs, const std::unordered_map<uint32_t, double>& splitters,
                     bool pi_buffers = false, bool pi_splitters = false, bool po_buffers = true )
      : gate_costs( gate_costs ), fanout_cc( splitters ), pi_buffers( pi_buffers ), pi_splitters( pi_splitters ), po_buffers( po_buffers ) {}

  template<typename Ntk, typename LevelMap, typename PoLevelMap>
  double operator()( const Ntk& ntk, const LevelMap& level_of_node, const PoLevelMap& po_level_of_node )
  {
    fanout_view dest_fv{ ntk };
    auto gate_cost = 0.0;
    auto fanout_net_cost = 0.0;

    std::vector<node<Ntk>> nodes;
    if ( pi_splitters )
    {
      dest_fv.foreach_pi( [&]( auto n ) { nodes.push_back( n ); } );
    }

    dest_fv.foreach_gate( [&]( auto n ) { nodes.push_back( n ); } );

    // dest_fv.foreach_node( [&]( auto n ) {
    //   if ( dest_fv.is_constant( n ) )
    //   {
    //     return;
    //   }

    //   if ( dest_fv.is_pi( n ) && !pi_splitters ) {
    //     return;
    //   }

    //   if ( n > 0u && dest_fv.is_maj( n ) )
    //   {
    //     internal_nodes.push_back( n );
    //   }
    // } );

    assert( po_level_of_node.size() > 0 );
    size_t critical_po_level = std::max_element( po_level_of_node.begin(), po_level_of_node.end(), []( auto n1, auto n2 ) { return n1.second < n2.second; } )->second;

    for ( auto n : nodes )
    {

      if ( !dest_fv.is_pi( n ) )
      {
        gate_cost += gate_costs.at( ntk.fanin_size( n ) );
      }

      if (ntk.fanout_size(n) == 0) continue;

      std::vector<uint32_t> rellev;

      dest_fv.foreach_fanout( n, [&]( auto fo ) {
        assert( level_of_node.at( fo ) > level_of_node.at( n ) );
        rellev.push_back( level_of_node.at( fo ) - level_of_node.at( n ) );
      } );

      uint32_t pos = 0u;
      while ( rellev.size() < dest_fv.fanout_size( n ) )
      {
        pos++;
        if ( po_buffers )
        {
          rellev.push_back( critical_po_level + 1 - level_of_node.at( n ) );
        }
        else
        {
          assert(level_of_node.count(n) > 0);
          assert(po_level_of_node.count(n) > 0);
          rellev.push_back( po_level_of_node.at( n ) + 1 - level_of_node.at( n ) );
        }
      }

      if ( rellev.size() > 1u || ( rellev.size() == 1u && rellev[0] > 0 ) )
      {
        std::sort( rellev.begin(), rellev.end() );
        auto net_cost = fanout_cc( rellev, dest_fv.is_pi( n ) && pi_buffers );
        if (net_cost == std::numeric_limits<double>::infinity()) {
          std::cerr << fmt::format("[e] impossible to synthesize fanout net of node {} for relative levels [{}]\n", n, fmt::join(rellev, " "));
          std::abort();
        }
        fanout_net_cost += net_cost;
      } else {
        std::cerr << fmt::format("[e] invalid level assignment for node {} with levels [{}]\n", n, fmt::join(rellev, " "));
        std::abort();
      }
    }

    return gate_cost + fanout_net_cost;
  }

private:
  std::unordered_map<uint32_t, double> gate_costs;
  balanced_fanout_net_cost fanout_cc;
  bool pi_buffers, pi_splitters, po_buffers;
};

} // namespace mockturtle