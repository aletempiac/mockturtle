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
  \file network_conversion.hpp
  \brief Network conversion for generic network utils

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>

#include "../../networks/klut.hpp"
#include "../../networks/generic.hpp"
#include "../../utils/node_map.hpp"
#include "../../views/binding_view.hpp"
#include "../../views/topo_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class NtkSource, class NtkDest>
NtkDest create_copy_generic( NtkSource const& ntk, node_map<typename NtkDest::signal, NtkSource>& old2new )
{
  if constexpr ( has_add_binding_v<NtkSource> && has_add_binding_v<NtkDest> )
  {
    NtkDest res( ntk.get_library() );

    old2new[ntk.get_constant( false )] = res.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = res.get_constant( true );
    }
    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    return res;
  }
  else
  {
    NtkDest res;

    old2new[ntk.get_constant( false )] = res.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = res.get_constant( true );
    }
    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    return res;
  }
}

} /* namespace detail */

/*! \brief Network convertion to generic network.
 *
 * This function converts a network from a mapped network
 * generated from a technology mapper (binding_view<klut_network>)
 * to a mapped generic_network.
 * This function is essential to do registers retiming.
 *
 * \param ntk Mapped network
 */
template<class Ntk>
binding_view<generic_network> generic_network_create_from_mapped( Ntk const& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_is_complemented_v<Ntk>, "NtkDest does not implement the is_complemented method" );
  static_assert( has_add_binding_v<Ntk>, "Ntk does not implement the add_binding method" );
  static_assert( has_set_as_latch_v<Ntk>, "Ntk does not implement the set_as_latch method" );

  using signal = typename generic_network::signal;

  node_map<signal, Ntk> old2new( ntk );
  binding_view<generic_network> res( ntk.get_library() );

  old2new[ntk.get_constant( false )] = res.get_constant( false );
  if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
  {
    old2new[ntk.get_constant( true )] = res.get_constant( true );
  }
  ntk.foreach_pi( [&]( auto const& n ) {
    old2new[n] = res.create_pi();
  } );

  topo_view topo{ ntk };

  topo.foreach_node( [&] ( auto const& n ) {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return true;

    std::vector<signal> children;
    
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old2new[f] );
    } );

    if ( ntk.is_as_latch( n ) )
    {
      auto const in_latch = res.create_box_input( children[0] );
      auto const latch = res.create_latch( in_latch );
      auto const latch_out = res.create_box_output( latch );
      res.add_binding( res.get_node( latch ), ntk.get_binding_index( n ) );
      old2new[n] = latch_out;
    }
    else
    {
      const auto f = res.create_node( children, ntk.node_function( n ) );
      res.add_binding( res.get_node( f ), ntk.get_binding_index( n ) );
      old2new[n] = f;
    }

    return true;
  } );

  ntk.foreach_po( [&]( auto const& f ) {
    res.create_po( old2new[f] );
  } );

  return res;
}

/*! \brief Network convertion from generic network.
 *
 * This function converts a mapped generic_network to a
 * mapped network (binding_view<klut_network>).
 *
 * \param ntk Mapped generic network
 */
binding_view<klut_network> mapped_create_from_generic_network( binding_view<generic_network> const& ntk )
{
  using signal = typename klut_network::signal;

  node_map<signal, generic_network> old2new( ntk );
  binding_view<klut_network> res( ntk.get_library() );

  old2new[ntk.get_constant( false )] = res.get_constant( false );
  if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
  {
    old2new[ntk.get_constant( true )] = res.get_constant( true );
  }
  ntk.foreach_pi( [&]( auto const& n ) {
    old2new[n] = res.create_pi();
  } );

  /* get latch gate (buffer) */
  uint32_t buf_id;
  for ( auto const& gate : ntk.get_library() )
  {
    if ( gate.num_vars == 1 && kitty::is_const0( kitty::cofactor0( gate.function, 0 ) ) )
    {
      buf_id = gate.id;
      break;
    }
  }

  topo_view topo{ ntk };

  topo.foreach_node( [&] ( auto const& n ) {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return true;

    /* remove not represented nodes */
    if ( ntk.is_box_input( n ) ||  ntk.is_box_output( n ) || ntk.is_po( n ) )
    {
      signal children;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        children = old2new[f];
      } );

      old2new[n] = children;
      return true;
    }

    std::vector<signal> children;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old2new[f] );
    } );

    const auto f = res.create_node( children, ntk.node_function( n ) );

    if ( ntk.is_latch( n ) )
    {
      res.add_binding( res.get_node( f ), buf_id );
    }
    else if ( ntk.has_binding( n ) )
    {
      res.add_binding( res.get_node( f ), ntk.get_binding_index( n ) );
    }
    old2new[n] = f;

    return true;
  } );

  ntk.foreach_po( [&]( auto const& f ) {
    res.create_po( old2new[f] );
  } );

  return res;
}

/*! \brief Cleans up dangling nodes on generic network type.
 *
 * This method reconstructs a network and omits all dangling nodes.  The
 * network types of the source and destination network are the same.
 *
 */
template<class NtkSource, class NtkDest = NtkSource>
NtkDest cleanup_dangling_generic( NtkSource const& ntk )
{
  using signal = typename NtkDest::signal;

  /* create a new network copy */
  node_map<signal, NtkSource> old2new( ntk );
  NtkDest res = detail::create_copy_generic<NtkSource, NtkDest>( ntk, old2new );

  topo_view topo{ ntk };
  topo.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_po( n ) )
    {
      res.create_po( old2new[ntk.get_fanin0( n )] );
    }
    else if ( ntk.is_box_input( n ) || ntk.is_box_output( n ) )
    {
      /* link to children */
      signal children;
      children = old2new[ntk.get_fanin0( n )];

      old2new[n] = children;
    }
    else if ( ntk.is_latch( n ) )
    {
      /* copy latch */
      signal children;
      children = old2new[ntk.get_fanin0( n )];

      auto const in_latch = res.create_box_input( res.make_signal( children ) );
      auto const latch = res.create_latch( in_latch );
      auto const latch_out = res.create_box_output( latch );

      old2new[n] = latch_out;
    }
    else
    {
      /* copy gate */
      std::vector<signal> children;

      ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( old2new[f] );
      } );

      const auto f = res.create_node( children, ntk.node_function( n ) );
      old2new[n] = f;

      if constexpr ( has_add_binding_v<NtkDest> && has_add_binding_v<NtkSource> )
      {
        res.add_binding( res.get_node( f ), ntk.get_binding_index( n ) );
      }
    }
  } );

  return res;
}

} // namespace mockturtle