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

#include "../views/fanout_view.hpp"


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

  /*! \brief Retiming iterations. */
  unsigned iterations{ 5 };
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
      bool improvement = true;
      for ( auto i = 0; i < _ps.iterations && improvement == true; ++i )
      {
        improvement = retime_area<true>();
      }
    }

    if ( !_ps.forward_only )
    {
      bool improvement = true;
      for ( auto i = 0; i < _ps.iterations && improvement == true; ++i )
      {
        improvement = retime_area<false>();
      }
    }

    _st.registers_post  = _ntk.num_latches();
  }

private:
  template<bool forward>
  bool retime_area()
  {
    init_values<forward>();

    auto min_cut = max_flow<forward>();

    if ( min_cut.size() >= ntk.num_latches() )
      return false;

    /* move latches */
    update_latches_position( min_cut );
  
    return false;
  }

  template<bool forward>
  std::vector<node> max_flow()
  {
    uint32_t flow = 0;

    _flow_path.reset();
    _ntk.incr_trav_id();

    /* run max flow from each register (capacity 1) */
    _ntk.foreach_latch( [&]( auto const& n ) {
      uint32_t local_flow;
      if constexpr ( forward )
      {
        local_flow = max_flow_forwards_compute_rec( _ntk.fanout( n )[0] );
      }
      else
      {
        node fanin = 0;
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = _ntk.get_node( f );
          return false;
        } );
        local_flow = max_flow_backwards_compute_rec( fanin );
      }

      flow += local_flow;

      if ( local_flow )
        _ntk.incr_trav_id();
    } );

    /* run reachability */
    _ntk.foreach_latch( [&]( auto const& n ) {
      uint32_t local_flow;
      if constexpr ( forward )
      {
        local_flow = max_flow_forwards_compute_rec( _ntk.fanout( n )[0] );
      }
      else
      {
        node fanin = 0;
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = _ntk.get_node( f );
          return false;
        } );
        local_flow = max_flow_backwards_compute_rec( fanin );
      }

      assert( local_flow == 0 );
    } );

    auto min_cut = get_min_cut();

    assert( check_min_cut<forward>( min_cut ) );

    std::cout << fmt::format( "Initial latches: {}\t Max flow: {}\t Latches after: {}\n", _st.registers_pre, flow, min_cut.size() );

    legalize_retiming<forward>( min_cut );

    std::cout << fmt::format( "Initial latches: {}\t Max flow: {}\t Latches after: {}\n\n", _st.registers_pre, flow, min_cut.size() );

    return min_cut;
  }

  // uint32_t max_flow_forwards_compute( node const& n )
  // {
  //   uint32_t found_path = 0;

  //   _ntk.set_visited( n, _ntk.trav_id() );

  //   /* node is not in a flow path */
  //   if ( _flow_path[n] == 0 )
  //   {
  //     _ntk.foreach_fanout( n, [&]( auto const& f ) {
  //       /* there is a path for flow */
  //       if ( max_flow_forwards_compute_rec( f ) )
  //       {
  //         _flow_path[n] = _ntk.node_to_index( f );
  //         found_path = 1;
  //         return false;
  //       }
  //       return true;
  //     } );

  //     if ( !found_path && _ntk.fanout_size( n ) != _ntk.fanout( n ).size() )
  //     {
  //       _flow_path[n] = sink_node;
  //       return 1;
  //     }
  //   }

  //   return found_path;
  // }

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

      // if ( _ntk.is_latch( n ) && _flow_path[n] != sink_node )
      // {
      //   if ( !_ntk.is_latch( _flow_path[n] ) && _ntk.visited( _flow_path[n] ) == _ntk.trav_id() )
      //     return true;
      // }

      if ( _ntk.value( n ) || _ntk.visited( _flow_path[n] ) != _ntk.trav_id() )
        min_cut.push_back( n );
      return true;
    } );

    return min_cut;
  }

  template<bool forward>
  void legalize_retiming( std::vector<node>& min_cut )
  {
    _ntk.clear_values();

    _ntk.foreach_latch( [&]( auto const& n ) {
      _ntk.set_value( n, 1 );
    } );

    for ( auto const& n : min_cut )
    {
      rec_mark_tfi( n );
    }

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
        node fanin = 0;
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = _ntk.get_node( f );
          return false;
        } );
        collect_cut_nodes_tfi( fanin, min_cut );
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
          _ntk.set_value( _ntk.get_node( f ), 1 );
        } );
      } );

      /* exclude reachable nodes from PIs from retiming */
      _ntk.foreach_pi( [&]( auto const& n ) {
        rec_mark_tfo( n );
      } );

      /* mark childrens of marked nodes */
      /* assume topological order */
      /* TODO: is it really needed? */
      _ntk.foreach_gate( [&]( auto const& n ) {
        if ( _ntk.value( n ) == 1 )
        {
          _ntk.foreach_fanin( n, [&]( auto const& f ) {
            _ntk.set_value( _ntk.get_node( f ), 1 );
          } );
        }
      } );
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

  void rec_mark_tfo( node const& n )
  {
    if ( _ntk.value( n ) == 1 )
      return;

    _ntk.set_value( n, 1 );
    _ntk.foreach_fanout( n, [&]( auto const& f ) {
      rec_mark_tfo( f );
    } );
  }

  void rec_mark_tfi( node const& n )
  {
    if ( _ntk.value( n ) == 1 )
      return;

    _ntk.set_value( n, 1 );
    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      rec_mark_tfi( _ntk.get_node( f ) );
    } );
  }

  template<bool forward>
  bool check_min_cut( std::vector<node> const& min_cut )
  {
    _ntk.incr_trav_id();

    for ( node const& n : min_cut )
    {
      _ntk.set_visited( n, _ntk.trav_id() );
    }

    bool check = true;
    _ntk.foreach_latch( [&]( auto const& n ) {
      if constexpr ( forward )
      {
        if ( !check_min_cut_rec<forward>( _ntk.fanout( n )[0] ) )
          check = false;
      }
      else
      {
        node fanin = 0;
        _ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = _ntk.get_node( f );
          return false;
        } );
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
        if ( !check_min_cut_rec<forward>( _ntk.get_node( f ) ) )
        {
          check = false;
        }
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
  uint32_t sink_node{UINT32_MAX};
};

} /* namespace detail */

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

  retime_stats st;

  using fanout_view_t = fanout_view<Ntk>;
  fanout_view_t fanout_view{ntk};

  detail::retime_impl p( fanout_view, ps, st );
  p.run();

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
