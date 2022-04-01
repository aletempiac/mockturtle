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
  \file sweep.hpp
  \brief Sweep utils for superconducting electronics

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>

#include <kitty/operations.hpp>

#include "../../utils/node_map.hpp"
#include "../../views/binding_view.hpp"
#include "../../views/depth_view.hpp"

namespace mockturtle
{

struct rsfq_balancing_sweep_params
{
  /*! \brief Maximum number of binate divisors to be considered. */
  bool allow_area_increase{ false };
};

namespace detail
{

template<class Ntk>
class rsfq_balancing_sweep_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit rsfq_balancing_sweep_impl( Ntk& ntk, rsfq_balancing_sweep_params const& ps )
      : _ntk( ntk ), _ps( ps )
  {
  }

  void run()
  {
    topo_view topo{ _ntk };
    topo.foreach_node( [this]( auto n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return;

      sweep_and_or( n );
      sweep_xor( n );
    } );
  }

private:
  void sweep_and_or( node const& n )
  {
    if ( _ntk.level( n ) <= 1 )
      return;

    bool is_and = true;
    if ( _ntk.node_function( n )._bits[0] == 0x8 )
    {
      /* AND */
      is_and = true;
    }
    else if ( _ntk.node_function( n )._bits[0] == 0xe )
    {
       /* OR */
      is_and = false;
    }
    else
    {
      return;
    }

    /* get children of top node, ordered by node level (ascending) */
    auto ocs = ordered_children( n );

    /* The output is an inverter too */
    auto const fanout = _ntk.fanout( n );
    bool out_inv = _ntk.fanout_size( n ) == 1 && fanout.size() == 1 && _ntk.node_function( fanout[0] )._bits[0] == 0x1;
    bool inv_0 = _ntk.node_function( _ntk.get_node( ocs[0] ) )._bits[0] == 0x1;

    /* depth of second child must be higher than depth of first child or area can be saved */
    bool depth_constraint = _ntk.level( _ntk.get_node( ocs[1] ) ) <= _ntk.level( _ntk.get_node( ocs[0] ) ) + 1;

    /* the critical children is not an inverter */
    if ( _ntk.node_function( _ntk.get_node( ocs[1] ) )._bits[0] != 0x1 )
      return;

    /* no optimization */
    if ( !out_inv && !inv_0 )
      return;

    signal child1 = get_child0( _ntk.get_node( ocs[1] ) );

    if ( out_inv && inv_0 )
    {
      signal child0 = get_child0( _ntk.get_node( ocs[0] ) );
      signal opt;
      if ( is_and )
        opt = _ntk.create_or( child0, child1 );
      else
        opt = _ntk.create_and( child0, child1 );
      _ntk.substitute_node( fanout[0], opt );
      _ntk.map_node( opt );
    }
    else if ( out_inv )
    {
      /* modification would increase the depth */
      if ( depth_constraint )
        return;

      signal opt;
      signal inv_s = search_inverter( _ntk.get_node( ocs[0] ) );
      
      /* modification would increase the area */
      if ( inv_s == 0 && _ntk.fanout_size( _ntk.get_node( ocs[1] ) ) != 1 )
        return;

      if ( inv_s == 0 )
      {
        inv_s = _ntk.create_not( _ntk.get_node( ocs[0] ) );
        _ntk.map_node( inv_s );
      }

      if ( is_and )
        opt = _ntk.create_or( inv_s, child1 );
      else
        opt = _ntk.create_and( inv_s, child1 );

      _ntk.substitute_node( fanout[0], opt );
      _ntk.map_node( opt );
    }
    else if ( _ntk.fanout_size( _ntk.get_node( ocs[0] ) ) == 1 && _ntk.fanout_size( _ntk.get_node( ocs[1] ) ) == 1 )
    {
      signal opt;
      signal child0 = get_child0( _ntk.get_node( ocs[0] ) );
      if ( is_and )
        opt = _ntk.create_or( child0, child1 );
      else
        opt = _ntk.create_and( child0, child1 );

      signal inv_s = _ntk.create_not( opt );
      _ntk.substitute_node( n, inv_s );
      _ntk.map_node( opt );
      _ntk.map_node( inv_s );
    }

    _ntk.update_levels();
  }

  void sweep_xor( node const& n )
  {
    if ( _ntk.level( n ) <= 1 )
      return;

    if ( _ntk.node_function( n )._bits[0] != 0x6 )
    {
      return;
    }

    /* get children of top node, ordered by node level (ascending) */
    auto ocs = ordered_children( n );

    /* The output is an inverter too */
    auto const fanout = _ntk.fanout( n );
    bool out_inv = _ntk.fanout_size( n ) == 1 && fanout.size() == 1 && _ntk.node_function( fanout[0] )._bits[0] == 0x1;
    bool inv_0 = _ntk.node_function( _ntk.get_node( ocs[0] ) )._bits[0] == 0x1;

    /* depth of second child must be higher than depth of first child or area can be saved */
    bool depth_constraint = _ntk.level( _ntk.get_node( ocs[1] ) ) <= _ntk.level( _ntk.get_node( ocs[0] ) ) + 1;

    /* the critical children is not an inverter */
    if ( _ntk.node_function( _ntk.get_node( ocs[1] ) )._bits[0] != 0x1 )
      return;

    /* no optimization */
    // if ( !out_inv && !inv_0 )
    //   return;

    signal child1 = get_child0( _ntk.get_node( ocs[1] ) );

    if ( out_inv )
    {
      /* remove inv_1 and out_inv */
      _ntk.replace_in_node( n, _ntk.get_node( ocs[1] ), child1 );

      if ( _ntk.decr_fanout_size( _ntk.get_node( ocs[1] ) ) )
        _ntk.take_out_node(  _ntk.get_node( ocs[1] ) );

      _ntk.substitute_node( fanout[0], n );
    }
    else if ( !out_inv && inv_0 )
    {
      signal child0 = get_child0( _ntk.get_node( ocs[0] ) );
      /* remove inv_0 and inv_1 */
      _ntk.replace_in_node( n, _ntk.get_node( ocs[0] ), child0 );
      _ntk.replace_in_node( n, _ntk.get_node( ocs[1] ), child1 );

      if ( _ntk.decr_fanout_size( _ntk.get_node( ocs[0] ) ) )
        _ntk.take_out_node(  _ntk.get_node( ocs[0] ) );
      if ( _ntk.decr_fanout_size( _ntk.get_node( ocs[1] ) ) )
        _ntk.take_out_node(  _ntk.get_node( ocs[1] ) );
    }
    else if ( !depth_constraint )
    {
      signal inv_s = search_inverter( _ntk.get_node( ocs[0] ) );

      /* modification would increase the area */
      if ( inv_s == 0 && _ntk.fanout_size( _ntk.get_node( ocs[1] ) ) != 1 )
        return;

      if ( inv_s == 0 )
      {
        inv_s = _ntk.create_not( _ntk.get_node( ocs[0] ) );
        _ntk.map_node( inv_s );
      }

      _ntk.replace_in_node( n, _ntk.get_node( ocs[0] ), inv_s );
      _ntk.replace_in_node( n, _ntk.get_node( ocs[1] ), child1 );

      _ntk.decr_fanout_size( _ntk.get_node( ocs[0] ) );
      if ( _ntk.decr_fanout_size( _ntk.get_node( ocs[1] ) ) )
        _ntk.take_out_node(  _ntk.get_node( ocs[1] ) );
    }
    else
    {
      return;
    }

    _ntk.update_levels();
  }

  inline std::array<signal, 2> ordered_children( node const& n ) const
  {
    std::array<signal, 2> children;
    _ntk.foreach_fanin( n, [&children]( auto const& f, auto i ) {
      children[i] = f;
    } );
    if ( _ntk.level( _ntk.get_node( children[0] ) ) > _ntk.level( _ntk.get_node( children[1] ) ) )
    {
      std::swap( children[0], children[1] );
    }
    return children;
  }

  inline signal get_child0( node const& n ) const
  {
    signal child;
    _ntk.foreach_fanin( n, [&child]( auto const& f ) {
      child = f;
      return false;
    } );
    return child;
  }

  inline signal search_inverter( node const& n ) const
  {
    signal inv_s = 0;
    _ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( _ntk.node_function( f )._bits[0] == 0x1 )
      {
        inv_s = _ntk.make_signal( f );
        return false;
      }
      return true;
    } );

    return inv_s;
  }

private:
  Ntk& _ntk;
  rsfq_balancing_sweep_params const& _ps;
};

} // namespace detail

template<class Ntk>
void rsfq_balancing_sweep( Ntk& ntk, rsfq_balancing_sweep_params const& ps = {} )
{
  detail::rsfq_balancing_sweep_impl p( ntk, ps );
  p.run();
}

}