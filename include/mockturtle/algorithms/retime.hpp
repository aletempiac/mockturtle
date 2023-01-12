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
  \file retime.hpp
  \brief Retiming

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <limits>

#include <fmt/format.h>
#include "../utils/stopwatch.hpp"
#include "../utils/node_map.hpp"
#include "../views/fanout_view.hpp"
#include "../views/topo_view.hpp"


namespace mockturtle
{
/*! \brief Parameters for retime.
 *
 * The data structure `retime_params` holds configurable parameters
 * with default arguments for `retime`.
 */
struct retime_params
{
  /*! \brief Do forward only retiming. */
  bool forward_only{ false };

  /*! \brief Do backward only retiming. */
  bool backward_only{ false };

  /*! \brief Retiming max iterations. */
  uint32_t iterations{ UINT32_MAX };

  /*! \brief Frontier based retiming. */
  bool frontier_retiming{ false };

  /*! \brief Be verbose */
  bool verbose { false };
};

/*! \brief Statistics for retiming.
 *
 * The data structure `retime_stats` provides data collected by running
 * `retime`.
 */
struct retime_stats
{
  /*! \brief Initial number of registers. */
  uint32_t registers_pre{ 0 };

  /*! \brief Number of registers after retime. */
  uint32_t registers_post{ 0 };

  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] Initial registers = {:7d}\t Final registers = {:7d}\n", registers_pre, registers_post );
    std::cout << fmt::format( "[i] Total runtime   = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

template<class Ntk>
class retime_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  static constexpr uint32_t sink_node = UINT32_MAX;

public:
  explicit retime_impl( Ntk& ntk, retime_params const& ps, retime_stats& st )
    : _ntk( ntk ),
      _ps( ps ),
      _st( st ),
      _flow_path( ntk )
  {}

public:
  void run()
  {
    stopwatch t( _st.time_total );

    _st.registers_pre =  _ntk.num_latches();

    if ( !_ps.backward_only )
    {
      if ( _ps.frontier_retiming )
      {
        init_frontiers<true>();
      }
      bool improvement = true;
      for ( auto i = 0; i < _ps.iterations && improvement == true; ++i )
      {
        improvement = retime_area<true>( i + 1 );
      }
    }

    if ( !_ps.forward_only )
    {
      if ( _ps.frontier_retiming )
      {
        init_frontiers<false>();
      }
      bool improvement = true;
      for ( auto i = 0; i < _ps.iterations && improvement == true; ++i )
      {
        improvement = retime_area<false>( i + 1 );
      }
    }

    _st.registers_post  = _ntk.num_latches();
  }

private:
  template<bool forward>
  bool retime_area( uint32_t iteration )
  {
    auto const num_latches_pre = _ntk.num_latches();

    init_values<forward>();

    auto min_cut = max_flow<forward>( iteration );

    if ( _ps.verbose )
    {
      float latch_improvement = ( (float) num_latches_pre - min_cut.size() ) / num_latches_pre * 100; 
      std::cout << fmt::format( "[i] Retiming {}\t pre = {:7d}\t post = {:7d}\t improvement = {:>5.2f}%\n",
                                 forward ? "forward" : "backward", num_latches_pre, min_cut.size(), latch_improvement );
    }

    if ( min_cut.size() >= num_latches_pre )
      return false;

    /* move latches */
    update_latches_position<forward>( min_cut, iteration );
  
    return true;
  }

  // template<bool forward>
  // void retime_area_frontiers()
  // {
  //   init_frontiers<forward>();

  //   uint32_t iterations = _frontiers.size();

  //   for ( auto i = 0u; i < iterations; ++i )
  //   {
  //     auto const num_latches_pre = _frontiers[i].size();
  //     init_values<forward>();

  //     auto min_frontier = max_flow_frontier<forward>( i );

  //     if ( _ps.verbose )
  //     {
  //       float latch_improvement = ( (float) num_latches_pre - min_frontier.size() ) / num_latches_pre * 100; 
  //       std::cout << fmt::format( "[i] Retiming {}\t pre = {:7d}\t post = {:7d}\t improvement = {:>5.2f}%\n",
  //                                 forward ? "forward" : "backward", num_latches_pre, min_frontier.size(), latch_improvement );
  //     }

  //     if ( min_frontier.size() < num_latches_pre )
  //     {
  //       /* move latches */
  //       update_latches_position_frontier<forward>( min_frontier, i );
  //     }
  //   }
  //   // std::cout << "--------------\n";
  // }

  template<bool forward>
  void init_frontiers()
  {
    _ntk.clear_values2();
    topo_view topo{_ntk};

    uint32_t max_frontiers = 0;

    /* assign frontiers to latches */
    if constexpr ( forward )
    {
      std::vector<node> topo_order;
      topo_order.reserve( topo.size() );

      topo.foreach_node( [&]( auto n ) {
        topo_order.push_back( n );
      } );

      /* iterate in reverse topological order */
      for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
      {
        if ( _ntk.is_pi( *it ) || _ntk.is_constant( *it ) )
          break;

        uint32_t level = _ntk.value( *it );
        if ( _ntk.is_latch( *it ) )
        {
          /* increase the frontier and set new value */
          ++level;
          _ntk.set_value2( *it, level );
        }

        max_frontiers = std::max( max_frontiers, level );

        _ntk.foreach_fanin( *it, [&]( auto const& f ) {
          uint32_t level_fanin = _ntk.value2( _ntk.get_node( f ) );
          _ntk.set_value2( _ntk.get_node( f ), std::min( level, level_fanin ) );
        } );
      }
    }
    else
    {
      topo.foreach_node( [&]( auto const& n ) {
        if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
          return true;

        uint32_t level = UINT32_MAX;
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          level = std::min( level, _ntk.value2( _ntk.get_node( f ) ) );
        } );

        if ( _ntk.is_latch( n ) )
        {
          /* increase the frontier */
          ++level;
        }

        max_frontiers = std::max( max_frontiers, level );

        _ntk.set_value2( n, level );
        return true;
      } );
    }

    // std::cout << "frontiers = " << max_frontiers << "\n";
  }

  template<bool forward>
  std::vector<node> max_flow( uint32_t iteration )
  {
    uint32_t flow = 0;

    _flow_path.reset();
    _ntk.incr_trav_id();

    /* run max flow from each register (capacity 1) */
    _ntk.foreach_latch( [&]( auto const& n ) {
      if ( _ps.frontier_retiming && ( _ntk.value2( n ) > iteration ) )
        return true;

      uint32_t local_flow;
      if constexpr ( forward )
      {
        local_flow = max_flow_forwards_compute_rec( _ntk.fanout( n )[0] );
      }
      else
      {
        node fanin = _ntk.get_node( _ntk.get_fanin0( n ) );
        local_flow = max_flow_backwards_compute_rec( fanin );
      }

      flow += local_flow;

      if ( local_flow )
        _ntk.incr_trav_id();

      return true;
    } );

    /* run reachability */
    _ntk.incr_trav_id();
    _ntk.foreach_latch( [&]( auto const& n ) {
      if ( _ps.frontier_retiming && ( _ntk.value2( n ) > iteration  ) )
        return true;

      uint32_t local_flow;
      if constexpr ( forward )
      {
        local_flow = max_flow_forwards_compute_rec( _ntk.fanout( n )[0] );
      }
      else
      {
        node fanin = _ntk.get_node( _ntk.get_fanin0( n ) );
        local_flow = max_flow_backwards_compute_rec( fanin );
      }

      assert( local_flow == 0 );
      return true;
    } );

    auto min_cut = get_min_cut();

    assert( check_min_cut<forward>( min_cut, iteration ) );

    // std::cout << fmt::format( "Initial latches: {}\t Max flow: {}\t Latches after: {}\n", _ntk.num_latches(), flow, min_cut.size() );

    legalize_retiming<forward>( min_cut, iteration );

    // std::cout << fmt::format( "Initial latches: {}\t Max flow: {}\t Latches after: {}\n\n", _ntk.num_latches(), flow, min_cut.size() );

    return min_cut;
  }

  uint32_t max_flow_forwards_compute_rec( node const& n )
  {
    uint32_t found_path = 0;

    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return 0;

    _ntk.set_visited( n, _ntk.trav_id() );

    /* node is not in a flow path */
    if ( _flow_path[n] == 0 )
    {
      /* cut boundary (sink) */
      if ( _ntk.value( n ) )
      {
        _flow_path[n] = sink_node;
        return 1;
      }

      _ntk.foreach_fanout( n, [&]( auto const& f ) {
        /* there is a path for flow */
        if ( max_flow_forwards_compute_rec( f ) )
        {
          _flow_path[n] = _ntk.node_to_index( f );
          found_path = 1;
          return false;
        }
        return true;
      } );

      return found_path;
    }

    /* path has flow already, find alternative path from fanin with flow */
    node fanin_flow = 0;
    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( _ntk.is_constant( _ntk.get_node( f ) ) )
        return true;
      if ( _flow_path[f] == _ntk.node_to_index( n ) )
      {
        fanin_flow = _ntk.get_node( f );
        return false;
      }
      return true;
    } );

    if ( fanin_flow == 0 )
      return 0;

    /* augment path */
    _ntk.foreach_fanout( fanin_flow, [&]( auto const& f ) {
      /* there is a path for flow */
      if ( max_flow_forwards_compute_rec( f ) )
      {
        _flow_path[fanin_flow] = _ntk.node_to_index( f );
        found_path = 1;
        return false;
      }
      return true;
    } );

    if ( found_path )
      return 1;

    if ( max_flow_forwards_compute_rec( fanin_flow ) )
    {
      _flow_path[fanin_flow] = 0;
      return 1;
    }

    return 0;
  }

  uint32_t max_flow_backwards_compute_rec( node const& n )
  {
    uint32_t found_path = 0;

    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return 0;

    _ntk.set_visited( n, _ntk.trav_id() );

    /* node is not in a flow path */
    if ( _flow_path[n] == 0 )
    {
      /* cut boundary (sink) */
      if ( _ntk.value( n ) )
      {
        _flow_path[n] = sink_node;
        return 1;
      }

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( _ntk.is_constant( _ntk.get_node( f ) ) )
          return true;
        /* there is a path for flow */
        if ( max_flow_backwards_compute_rec( _ntk.get_node( f ) ) )
        {
          _flow_path[n] = _ntk.node_to_index( _ntk.get_node( f ) );
          found_path = 1;
          return false;
        }
        return true;
      } );

      return found_path;
    }

    /* path has flow already, find alternative path from fanout with flow */
    node fanout_flow = 0;
    _ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( _flow_path[f] == _ntk.node_to_index( n ) )
      {
        fanout_flow = _ntk.get_node( f );
        return false;
      }
      return true;
    } );

    if ( fanout_flow == 0 )
      return 0;

    /* augment path */
    _ntk.foreach_fanin( fanout_flow, [&]( auto const& f ) {
      if ( _ntk.is_constant( _ntk.get_node( f ) ) )
        return true;
      /* there is a path for flow */
      if ( max_flow_backwards_compute_rec( _ntk.get_node( f ) ) )
      {
        _flow_path[fanout_flow] = _ntk.node_to_index( _ntk.get_node( f ) );
        found_path = 1;
        return false;
      }
      return true;
    } );

    if ( found_path )
      return 1;

    if ( max_flow_backwards_compute_rec( fanout_flow ) )
    {
      _flow_path[fanout_flow] = 0;
      return 1;
    }

    return 0;
  }

  std::vector<node> get_min_cut()
  {
    std::vector<node> min_cut;
    min_cut.reserve( _ntk.num_latches() );

    _ntk.foreach_node( [&]( auto const& n ) {
      if ( _flow_path[n] == 0 )
        return true;
      if ( _ntk.visited( n ) != _ntk.trav_id() )
        return true;

      if ( _ntk.value( n ) || _ntk.visited( _flow_path[n] ) != _ntk.trav_id() )
        min_cut.push_back( n );
      return true;
    } );

    return min_cut;
  }

  template<bool forward>
  void legalize_retiming( std::vector<node>& min_cut, uint32_t iteration )
  {
    _ntk.clear_values();

    _ntk.foreach_latch( [&]( auto const& n ) {
      if ( _ps.frontier_retiming && ( _ntk.value2( n ) > iteration  ) )
      {
        _ntk.set_value( _ntk.fanout( n )[0], 2 );
      }
      else
      {
        _ntk.set_value( _ntk.fanout( n )[0], 1 );
      }
    } );

    for ( auto const& n : min_cut )
    {
      rec_mark_tfi( n );
    }

    // _ntk.foreach_latch( [&]( auto const& n ) {
    //   if ( _ntk.value( n ) == 2 )
    //   {
    //     rec_mark_tfi( _ntk.get_fanin0( n ) );
    //   }
    // } );

    min_cut.clear();

    if constexpr ( forward )
    {
      _ntk.foreach_gate( [&]( auto const& n ) {
        if ( _ntk.value( n ) == 1 )
        {
          /* if is sink or before a register */
          _ntk.foreach_fanout( n, [&]( auto const& f ) {
            if ( _ntk.value( f ) != 1 )
            {
              min_cut.push_back( n );
              return false;
            }
            return true;
          } );
        }
      } );
    }
    else
    {
      _ntk.incr_trav_id();
      _ntk.foreach_latch( [&]( auto const& n ) {
        if ( _ps.frontier_retiming && ( _ntk.value2( n ) > iteration  ) )
          return true;

        node fanin = _ntk.get_node( _ntk.get_fanin0( n ) );
        collect_cut_nodes_tfi( fanin, min_cut );
        return true;
      } );
      _ntk.foreach_node( [&]( auto const& n ) {
        if ( _ntk.visited( n ) == _ntk.trav_id() )
          _ntk.set_value( n, 1 );
        else
          _ntk.set_value( n, 0 );
      } );
      for ( auto const& n : min_cut )
        _ntk.set_value( n, 0 );
    }
  }

  // template<bool forward>
  // void legalize_retiming_frontier( std::vector<node>& min_frontier, uint32_t iteration )
  // {
  //   _ntk.clear_values();

  //   /* set all latches to symbolic value of 2 to limit the marking recursion */
  //   _ntk.foreach_latch( [&]( auto const& n ) {
  //     _ntk.set_value( n, 2 );
  //   } );

  //   /* mark frontier */
  //   for ( auto n : _frontiers[iteration] )
  //   {
  //     _ntk.set_value( n, 1 );
  //   }

  //   for ( auto const& n : min_frontier )
  //   {
  //     rec_mark_tfi( n );
  //   }

  //   _ntk.foreach_latch( [&]( auto const& n ) {
  //     if ( _ntk.value( n ) == 2 )
  //       rec_mark_tfi( _ntk.get_node( _ntk.get_fanin0( n ) ) );
  //   } );

  //   min_frontier.clear();

  //   if constexpr ( forward )
  //   {
  //     _ntk.foreach_gate( [&]( auto const& n ) {
  //       if ( _ntk.value( n ) == 1 )
  //       {
  //         /* if is sink or before a register */
  //         _ntk.foreach_fanout( n, [&]( auto const& f ) {
  //           if ( _ntk.value( f ) == 0 )
  //           {
  //             min_frontier.push_back( n );
  //             return false;
  //           }
  //           return true;
  //         } );
  //       }
  //     } );
  //   }
  //   else
  //   {
  //     _ntk.incr_trav_id();
  //     for ( auto n : _frontiers[iteration] )
  //     {
  //       node fanin = _ntk.get_node( _ntk.get_fanin0( n ) );
  //       collect_cut_nodes_tfi( fanin, min_frontier );
  //     }
  //     _ntk.foreach_node( [&]( auto const& n ) {
  //       if ( _ntk.visited( n ) == _ntk.trav_id() )
  //         _ntk.set_value( n, 1 );
  //       else
  //         _ntk.set_value( n, 0 );
  //     } );
  //     for ( auto const& n : min_frontier )
  //       _ntk.set_value( n, 0 );
  //   }
  // }

  void collect_cut_nodes_tfi( node const& n, std::vector<node>& min_cut )
  {
    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return;

    _ntk.set_visited( n, _ntk.trav_id() );

    if ( _ntk.value( n ) )
    {
      min_cut.push_back( n );
      return;
    }

    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( _ntk.is_constant( _ntk.get_node( f ) ) )
        return;
      collect_cut_nodes_tfi( _ntk.get_node( f ), min_cut );
    } );
  }

  template<bool forward>
  void init_values()
  {
    _ntk.clear_values();

    /* marks the frontiers */
    if constexpr ( forward )
    {
      /* mark POs as sink */
      _ntk.foreach_po( [&]( auto const& f ) {
        _ntk.set_value( _ntk.get_node( f ), 1 );
      } );

      /* mark registers as sink */
      _ntk.foreach_latch( [&]( auto const& n ) {
        _ntk.set_value( n, 1 );
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          if ( _ntk.is_constant( _ntk.get_node( f ) ) )
            return;
          _ntk.set_value( _ntk.get_node( f ), 1 );
        } );
      } );

      /* exclude reachable nodes from PIs from retiming */
      _ntk.foreach_pi( [&]( auto const& n ) {
        rec_mark_tfo( n );
      } );

      /* mark childrens of marked nodes */
      std::vector<node> to_mark;
      to_mark.reserve( 200 );
      _ntk.foreach_gate( [&]( auto const& n ) {
        if ( _ntk.value( n ) == 1 )
        {
          _ntk.foreach_fanin( n, [&]( auto const& f ) {
            if ( _ntk.is_constant( _ntk.get_node( f ) ) )
              return;
            if ( _ntk.value( _ntk.get_node( f ) ) == 0 )
              to_mark.push_back( _ntk.get_node( f ) );
          } );
        }
      } );
      for ( auto const& n : to_mark )
      {
        _ntk.set_value( n, 1 );
      }
    }
    else
    {
      /* mark PIs as sink */
      _ntk.foreach_pi( [&]( auto const& n ) {
        _ntk.set_value( n, 1 );
      } );

      /* mark registers as sink */
      _ntk.foreach_latch( [&]( auto const& n ) {
        _ntk.set_value( n, 1 );
        _ntk.foreach_fanout( n, [&]( auto const& f ) {
          _ntk.set_value( f, 1 );
        } );
      } );

      /* exclude reachable nodes from POs from retiming */
      _ntk.foreach_po( [&]( auto const& f ) {
        rec_mark_tfi( _ntk.get_node( f ) );
      } );
    }
  }

  template<bool forward>
  void update_latches_position( std::vector<node> const& min_cut, uint32_t iteration )
  {
    _ntk.incr_trav_id();

    /* create new latches and mark the ones to reuse */
    for ( auto const& n : min_cut )
    {
      if constexpr ( forward )
      {
        if ( _ntk.is_box_output( n ) )
        {
          /* reuse the current latch */
          auto latch = _ntk.get_node( _ntk.get_fanin0( n ) );
          auto in_latch = _ntk.get_node( _ntk.get_fanin0( latch ) );
          auto in_in_latch = _ntk.get_node( _ntk.get_fanin0( in_latch ) );

          /* check for marked fanouts to connect to latch input */
           auto fanout = _ntk.fanout( n );
          for ( auto const& f : fanout )
          {
            if ( _ntk.value( f ) )
            {
              _ntk.replace_in_node( f, n, in_in_latch );
              _ntk.decr_fanout_size( n );
            }
          }

          _ntk.set_visited( latch, _ntk.trav_id() );
          // _ntk.foreach_fanin( n, [&]( auto const& f ) {
          //   if ( _ntk.is_latch( _ntk.get_node( f ) ) )
          //     _ntk.set_visited( _ntk.get_node( f ), _ntk.trav_id() );
          // } );

          /* check for not marked fanouts */
          // _ntk.foreach_fanout( n, [&]( auto const& f ) {
          //   if ( !_ntk.value( f ) )
          //     std::cout << "not marked fanout\n";
          // } );
        }
        else
        {
          /* create a new latch */
          auto const in_latch = _ntk.create_box_input( _ntk.make_signal( n ) );
          auto const latch = _ntk.create_latch( in_latch );
          auto const latch_out = _ntk.create_box_output( latch );

          /* replace in n fanout */
          auto fanout = _ntk.fanout( n );
          for ( auto const& f : fanout )
          {
            if ( f != _ntk.get_node( in_latch ) && !_ntk.value( f ) )
            {
              _ntk.replace_in_node( f, n, latch_out );
              _ntk.decr_fanout_size( n );
            }
          }

          _ntk.set_visited( _ntk.get_node( latch ), _ntk.trav_id() );
        }
      }
      else
      {
        if ( _ntk.is_box_input( n ) )
        {
          _ntk.foreach_fanout( n, [&]( auto const& f ) {
            _ntk.set_visited( f, _ntk.trav_id() );
          } );
        }
        else
        {
          /* create a new latch */
          auto const in_latch = _ntk.create_box_input( _ntk.make_signal( n ) );
          auto const latch = _ntk.create_latch( in_latch );
          auto const latch_out = _ntk.create_box_output( latch );

          /* replace in n fanout */
          auto fanout = _ntk.fanout( n );
          for ( auto const& f : fanout )
          {
            if ( f != _ntk.get_node( in_latch ) && _ntk.value( f ) )
            {
              _ntk.replace_in_node( f, n, latch_out );
              _ntk.decr_fanout_size( n );
            }
          }

          _ntk.set_visited( _ntk.get_node( latch ), _ntk.trav_id() );
        }
      }
    }

    /* remove retimed latches */
    _ntk.foreach_latch( [&]( auto const& n )
    {
      if ( _ps.frontier_retiming && ( _ntk.value2( n ) > iteration  ) )
        return true;

      if ( _ntk.visited( n ) == _ntk.trav_id() )
        return true;

      node latch_out;
      node latch_in = _ntk.get_node( _ntk.get_fanin0( n ) );
      signal latch_in_in = _ntk.get_fanin0( latch_in );

      _ntk.foreach_fanout( n, [&]( auto const& f )
      {
        latch_out = f;
      } );

      auto latch_fanout = _ntk.fanout_size( latch_out );
      auto fanin_fanout = _ntk.fanout_size( _ntk.get_node( latch_in_in ) );
      auto fanin_type = _ntk.is_box_output( _ntk.get_node( latch_in_in ) );

      _ntk.substitute_node( latch_out, latch_in_in );

      return true;
    } );
    

    /* cleanup dangling nodes */
    // _ntk = fanout_view( cleanup_dangling_generic<fanout_view<Ntk>, Ntk>( _ntk ) );
  }

  // template<bool forward>
  // void update_latches_position( std::vector<node> const& min_cut )
  // {
  //   _ntk.incr_trav_id();

  //   /* mark latches to reuse */
  //   for ( auto const& n : min_cut )
  //   {
  //     /* mark mincut nodes */
  //     _ntk.set_visited( n, _ntk.trav_id() );

  //     if constexpr ( forward )
  //     {
  //       if ( _ntk.is_box_output( n ) )
  //       {
  //         /* reuse the current latch */
  //         _ntk.foreach_fanin( n, [&]( auto const& f ) {
  //           if ( _ntk.is_latch( _ntk.get_node( f ) ) )
  //             _ntk.set_visited( _ntk.get_node( f ), _ntk.trav_id() );
  //         } );

  //         /* check for not marked fanouts */
  //         // _ntk.foreach_fanout( n, [&]( auto const& f ) {
  //         //   if ( !_ntk.value( f ) )
  //         //     std::cout << "not marked fanout\n";
  //         // } );
  //       }
  //     }
  //     else
  //     {
  //       if ( _ntk.is_box_input( n ) )
  //       {
  //         _ntk.foreach_fanout( n, [&]( auto const& f ) {
  //           if ( _ntk.is_latch( f ) )
  //             _ntk.set_visited( f, _ntk.trav_id() );
  //         } );
  //       }
  //     }
  //   }

  //   /* create a new network copy */
  //   node_map<signal, Ntk> old2new( _ntk );
  //   Ntk res = create_copy_retiming( old2new );

  //   _ntk.foreach_gate( [&]( auto const& n ) {
  //     if ( _ntk.is_po( n ) )
  //     {
  //       signal children;
  //       children = old2new[_ntk.get_fanin0( n )];
  //       res.create_po( children );
  //     }
  //     else if ( _ntk.is_box_input( n ) || _ntk.is_box_output( n ) )
  //     {
  //       /* link to children */
  //       signal children;
  //       children = old2new[_ntk.get_fanin0( n )];

  //       old2new[n] = children;
  //     }
  //     else if ( !_ntk.is_latch( n ) )
  //     {
  //       /* copy gate */
  //       std::vector<signal> children;

  //       _ntk.foreach_fanin( n, [&]( auto const& f ) {
  //         signal child = old2new[f];

  //         if ( _ntk.value( n ) ^ forward )
  //           child = res.get_fanin0
  //         children.push_back( child );
  //       } );

  //       const auto f = res.create_node( children, _ntk.node_function( n ) );
  //       old2new[n] = f;

  //       if constexpr ( has_add_binding_v<Ntk> )
  //       {
  //         res.add_binding( res.get_node( f ), _ntk.get_binding_index( n ) );
  //       }

  //       if ( _ntk.visited( n ) == _ntk.trav_id() )
  //       {
  //         /* create a new latch */
  //         auto const in_latch = res.create_box_input( f );
  //         auto const latch = res.create_latch( in_latch );
  //         auto const latch_out = res.create_box_output( latch );

  //         old2new[n] = latch_out;
  //       }
  //     }
  //     else
  //     {
  //       /* latch/register */
  //       if ( _ntk.visited( n ) == _ntk.trav_id() )
  //       {
  //         /* copy latch */
  //         signal children;
  //         children = old2new[_ntk.get_fanin0( n )];

  //         auto const in_latch = res.create_box_input( res.make_signal( children ) );
  //         auto const latch = res.create_latch( in_latch );
  //         auto const latch_out = res.create_box_output( latch );

  //         old2new[n] = latch_out;
  //       }
  //       else
  //       {
  //         old2new[n] = old2new[_ntk.get_fanin0( n )];
  //       }
  //     }
  //   } );

  //   _ntk = fanout_view( res );
  // }

  void rec_mark_tfo( node const& n )
  {
    if ( _ntk.value( n ) )
      return;

    _ntk.set_value( n, 1 );
    _ntk.foreach_fanout( n, [&]( auto const& f ) {
      rec_mark_tfo( f );
    } );
  }

  void rec_mark_tfi( node const& n )
  {
    if ( _ntk.value( n ) )
      return;

    _ntk.set_value( n, 1 );
    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( _ntk.is_constant( _ntk.get_node( f ) ) )
        return;
      rec_mark_tfi( _ntk.get_node( f ) );
    } );
  }

  template<bool forward>
  bool check_min_cut( std::vector<node> const& min_cut, uint32_t iteration )
  {
    _ntk.incr_trav_id();

    for ( node const& n : min_cut )
    {
      _ntk.set_visited( n, _ntk.trav_id() );
    }

    bool check = true;
    _ntk.foreach_latch( [&]( auto const& n ) {
      if ( _ps.frontier_retiming && _ntk.value2( n ) > iteration )
        return true;

      if constexpr ( forward )
      {
        if ( !check_min_cut_rec<forward>( _ntk.fanout( n )[0] ) )
          check = false;
      }
      else
      {
        node fanin =_ntk.get_node( _ntk.get_fanin0( n ) );
        if ( !check_min_cut_rec<forward>( fanin ) )
          check = false;
      }

      return check;
    } );

    return check;
  }

  template<bool forward>
  bool check_min_cut_rec( node const& n )
  {
    bool check = true;

    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return true;

    _ntk.set_visited( n, _ntk.trav_id() );

    if constexpr ( forward )
    {
      if ( _ntk.is_co( n ) )
      {
        check = false;
        return false;
      }

      _ntk.foreach_fanout( n, [&]( auto const& f ) {
        if ( !check_min_cut_rec<forward>( f ) )
        {
          check = false;
        }
      } );
    }
    else
    {
      if ( _ntk.is_ci( n ) )
      {
        check = false;
        return false;
      }

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( _ntk.is_constant( _ntk.get_node( f ) ) )
          return true;
        if ( !check_min_cut_rec<forward>( _ntk.get_node( f ) ) )
        {
          check = false;
        }
        return check;
      } );

      return check;
    }

    return check;
  }

private:
  Ntk& _ntk;
  retime_params const& _ps;
  retime_stats& _st;

  node_map<uint32_t, Ntk> _flow_path;
};

} /* namespace detail */

/*! \brief Retiming.
 *
 * This function implements a retiming algorithm for registers minimization.
 * The only supported network type is the `generic_network`.
 * The algorithm excecutes the retiming inplace.
 * 
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 *
 * \param ntk Network
 * \param ps Retiming params
 * \param pst Retiming statistics
 * 
 * The implementation of this algorithm was inspired by the
 * mapping command ``retime`` in ABC.
 */
template<class Ntk>
void retime( Ntk& ntk, retime_params const& ps = {}, retime_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_incr_value_v<Ntk>, "Ntk does not implement the incr_value method" );
  static_assert( has_decr_value_v<Ntk>, "Ntk does not implement the decr_value method" );
  static_assert( has_get_fanin0_v<Ntk>, "Ntk does not implement the get_fanin0 method" );

  retime_stats st;

  using fanout_view_t = fanout_view<Ntk>;
  fanout_view_t fanout_view{ntk};

  detail::retime_impl p( fanout_view, ps, st );
  p.run();

  if ( ps.verbose )
    st.report();

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
