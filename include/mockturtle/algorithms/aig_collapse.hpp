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
  \file aig_collapse.hpp
  \brief Collapse nodes in an AIG

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>
#include <random>
#include <algorithm>

#include "../networks/aig.hpp"
#include "../networks/multi_aig.hpp"
#include "../traits.hpp"
#include "../utils/node_map.hpp"

namespace mockturtle
{

struct aig_collapse_params
{
  /*! \brief Limit on the number of collapsed nodes. */
  uint32_t collapse_limit{ 32u };
};

namespace detail
{

template<class Ntk>
class aig_collapse_impl
{
public:
  static constexpr size_t storage_init_size = 30;
  using node = typename Ntk::node;
  using storage_t = std::vector<std::vector<signal<Ntk>>>;
  using collapsed_map = node_map<signal<multi_aig_network>, Ntk>;

public:
  aig_collapse_impl( Ntk const& ntk, aig_collapse_params const& ps )
      : ntk( ntk ), ps( ps ), storage( storage_init_size )
  {
  }

  multi_aig_network run()
  {
    multi_aig_network res;
    collapsed_map old2new( ntk );

    initialize_multi_aig_network( res, old2new );
    ntk.clear_values();
    children.reserve( ps.collapse_limit );

    for ( auto i = 0; i < storage_init_size; ++i )
      storage[i].reserve( 10 );

    /* collapse in reverse topo order */
    ntk.foreach_co( [&]( auto const& f ) {
      collapse_rec( res, old2new, ntk.get_node( f ), 0 );
      res.create_po( old2new[f] ^ ntk.is_complemented( f ) );
    } );

    return res;
  }

private:
  void collapse_rec( multi_aig_network& res, collapsed_map& old2new, node const& n, uint32_t level )
  {
    if ( ntk.is_ci( n ) || ntk.value( n ) > 1 )
      return;
    
    assert( ntk.value( n ) != 1 );
    ntk.incr_value( n );

    if ( level >= storage.size() )
    {
      storage.emplace_back( std::vector<signal<Ntk>>() );
      storage.back().reserve( 10 );
    }

    /* collect leaves of the AND tree */
    collect_leaves( n, storage[level] );

    /* constant false */
    if ( storage[level].size() == 0 )
    {
      old2new[n] = res.get_constant( false );
      ntk.set_visited( n, ntk.trav_id() );
      return;
    }

    /* recur over the leaves */
    for ( auto& f : storage[level] )
    {
      collapse_rec( res, old2new, ntk.get_node( f ), level + 1 );
    }

    assert( storage[level].size() > 1 );
    ntk.incr_value( n );

    /* create the multi-input AND node */
    children.clear();
    std::transform( storage[level].begin(),
                    storage[level].begin() + std::min( static_cast<uint32_t>( storage[level].size() ), ps.collapse_limit ),
                    std::back_inserter( children ),
                    [&]( auto const& f ) {
                      return old2new[f] ^ ntk.is_complemented( f );
                    }
                  );

    old2new[n] = res.create_nary_and( children );
    
    /* leaves number exceeds the limit: create an AND chain */
    if ( storage[level].size() > ps.collapse_limit )
    {
      uint32_t l = ps.collapse_limit;
      while ( l < storage[level].size() )
      {
        children.clear();
        std::transform( storage[level].begin() + l,
                        storage[level].begin() + std::min( static_cast<uint32_t>( storage[level].size() ), l + ps.collapse_limit - 1 ),
                        std::back_inserter( children ),
                        [&]( auto const& f ) {
                          return old2new[f] ^ ntk.is_complemented( f );
                        }
                      );
        children.push_back( old2new[n] );
        old2new[n] = res.create_nary_and( children );
        l += ps.collapse_limit - 1;
      }
    }

    /* clean leaves storage */
    storage[level].clear();
}

  void collect_leaves( node const& n, std::vector<signal<Ntk>>& leaves )
  {
    ntk.incr_trav_id();

    int ret = collect_leaves_rec( ntk.make_signal( n ), leaves, true );

    /* check for constant false */
    if ( ret < 0 )
    {
      leaves.clear();
    }
  }

  int collect_leaves_rec( signal<Ntk> const& f, std::vector<signal<Ntk>>& leaves, bool is_root )
  {
    node n = ntk.get_node( f );

    /* check if already visited */
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      for ( signal<Ntk> const& s : leaves )
      {
        if ( ntk.get_node( s ) != n )
          continue;

        if ( s == f )
          return 1;   /* same polarity: duplicate */
        else
          return -1;  /* opposite polarity: const0 */
      }

      return 0;
    }

    /* set as leaf if signal is complemented or is a CI or has a multiple fanout */
    if ( !is_root && ( ntk.is_complemented( f ) || ntk.is_ci( n ) || ntk.fanout_size( n ) > 1 ) )
    {
      leaves.push_back( f );
      ntk.set_visited( n, ntk.trav_id() );
      return 0;
    }

    int ret = 0;
    ntk.foreach_fanin( n, [&]( auto const& child ) {
      ret |= collect_leaves_rec( child, leaves, false );
    } );

    return ret;
  }

  void insert_node_sorted( std::vector<signal<Ntk>>& leaves, signal<Ntk> const& f )
  {
    node n = ntk.get_node( f );

    /* check uniqueness */
    for ( auto const& s : leaves )
    {
      if ( s == f )
        return;
    }

    leaves.push_back( f );
    for ( size_t i = leaves.size() - 1; i > 0; --i )
    {
      auto& s2 = leaves[i - 1];

      if ( ntk.level( ntk.get_node( s2 ) ) < ntk.level( n ) )
      {
        std::swap( s2, leaves[i] );
      }
      else
      {
        break;
      }
    }
  }

  void initialize_multi_aig_network( multi_aig_network& dest, collapsed_map& old2new )
  {
    old2new[ntk.get_node( ntk.get_constant( false ) )] = dest.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
      old2new[ntk.get_node( ntk.get_constant( true ) )] = dest.get_constant( true );

    ntk.foreach_ci( [&]( auto const& n ) {
      old2new[n] = dest.create_pi();
    } );
  }

private:
  Ntk const& ntk;
  aig_collapse_params const& ps;

  storage_t storage;
  std::vector<signal<multi_aig_network>> children;
};

} /* namespace detail */

/*! \brief AIG collapse.
 *
 * This method collapses AND2 nodes in an AIG into
 * multi-input ANDs. The maximum number of inputs can
 * be limited using the parameter `collapse_limit`.
 * It returns the resulted network as a klut.
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
 * - `has_and`
 */
template<class Ntk>
multi_aig_network aig_collapse( Ntk const& ntk, aig_collapse_params const& ps = {} )
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

  detail::aig_collapse_impl p( ntk, ps );
  multi_aig_network res = p.run();

  return res;
}

} /* namespace mockturtle */
