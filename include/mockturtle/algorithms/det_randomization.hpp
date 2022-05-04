/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file det_randomization.hpp
  \brief Randomizes the topological ordering of a network

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>
#include <random>
#include <algorithm>

#include "../traits.hpp"
#include "../utils/node_map.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class det_randomize_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  det_randomize_impl( Ntk& ntk, std::default_random_engine::result_type const& seed )
      : ntk( ntk ), seed( seed )
  {
  }

  Ntk run()
  {
    Ntk dest;
    node_map<signal, Ntk> old2new( ntk );

    ntk.incr_trav_id();
    ntk.incr_trav_id();

    old2new[ntk.get_constant( false )] = dest.get_constant( false );
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = dest.get_constant( true );
      ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
    }

    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = dest.create_pi();
      ntk.set_visited( n, ntk.trav_id() );
    } );

    /* To try: shuffle POs visiting order */
    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() )
        return;
      
      topo_rec( dest, old2new, ntk.get_node( f ) );
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
        dest.create_po( dest.create_not( old2new[f] ) );
      else
       dest.create_po( old2new[f] );
    } );

    return dest;
  }

private:
  void topo_rec( Ntk& dest, node_map<signal, Ntk>& old2new, node const& n )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    
    assert( ntk.visited( n ) != ntk.trav_id() - 1 );

    ntk.set_visited( n, ntk.trav_id() - 1 );

    std::vector<node> fanins;
    fanins.reserve( ntk.fanin_size( n ) );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      fanins.push_back( ntk.get_node( f ) );
    } );

    std::shuffle( fanins.begin(), fanins.end(), std::default_random_engine( seed++ ) );

    for ( node const& g : fanins )
    {
      topo_rec( dest, old2new, g );
    }

    std::vector<signal> children;
    children.reserve( ntk.fanin_size( n ) );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
        children.push_back( ntk.create_not( old2new[f] ) );
      else
        children.push_back( old2new[f] );
    } );

    old2new[n] = dest.clone_node( ntk, n, children );

    ntk.set_visited( n, ntk.trav_id() );
  }

private:
  Ntk& ntk;
  std::default_random_engine::result_type seed;
};

} /* namespace detail */

/*! \brief Topological ordering randomization.
 *
 * This method sorts the topological order of a network
 * using a deterministic random function and cleans up
 * dangling nodes.
 *
   \verbatim embed:rst

   .. note::

      This method returns the cleaned up network as a return value.  It does
      *not* modify the input network.
   \endverbatim
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_node`
 * - `foreach_pi`
 * - `foreach_po`
 * - `clone_node`
 * - `is_pi`
 * - `is_constant`
 */
template<class Ntk>
Ntk det_randomize( Ntk& ntk, std::default_random_engine::result_type const seed = 1 )
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
  static_assert( has_clone_node_v<Ntk>, "Ntk does not implement the clone_node method" );
  static_assert( has_create_pi_v<Ntk>, "Ntk does not implement the create_pi method" );
  static_assert( has_create_po_v<Ntk>, "Ntk does not implement the create_po method" );
  static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  detail::det_randomize_impl p( ntk, seed );
  return p.run();
}

} // namespace mockturtle