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

template<class NtkSource, class NtkDest>
class aig_collapse_impl
{
public:
  static constexpr size_t storage_init_size = 30;
  using node = typename NtkSource::node;
  using storage_t = std::vector<std::vector<signal<NtkSource>>>;
  using collapsed_map = node_map<signal<NtkDest>, NtkSource>;

public:
  aig_collapse_impl( NtkSource const& ntk, aig_collapse_params const& ps )
      : ntk( ntk ), ps( ps ), storage( storage_init_size )
  {
  }

  NtkDest run()
  {
    NtkDest res;
    collapsed_map old2new( ntk );

    initialize_multi_aig_network( res, old2new );
    ntk.clear_values();
    children.reserve( ps.collapse_limit );

    for ( auto i = 0; i < storage_init_size; ++i )
      storage[i].reserve( 10 );

    /* collapse in reverse topo order */
    ntk.foreach_po( [&]( auto const& f ) {
      collapse_rec( res, old2new, ntk.get_node( f ), 0 );
      res.create_po( old2new[f] ^ ntk.is_complemented( f ) );
    } );

    if constexpr ( has_foreach_ri_v<NtkSource> )
    {
      static_assert( has_create_ri_v<NtkDest>, "NtkDest does not implement the create_ri method " );
      ntk.foreach_ri( [&]( auto const& f ) {
        collapse_rec( res, old2new, ntk.get_node( f ), 0 );
        res.create_ri( old2new[f] ^ ntk.is_complemented( f ) );
      } );
    }

    return res;
  }

private:
  void collapse_rec( NtkDest& res, collapsed_map& old2new, node const& n, uint32_t level )
  {
    if ( ntk.is_ci( n ) || ntk.value( n ) > 1 )
      return;
    
    assert( ntk.value( n ) != 1 );
    ntk.incr_value( n );

    if ( level >= storage.size() )
    {
      storage.emplace_back( std::vector<signal<NtkSource>>() );
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

    /* sort */
    std::sort( storage[level].begin(), storage[level].end() );

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

  void collect_leaves( node const& n, std::vector<signal<NtkSource>>& leaves )
  {
    ntk.incr_trav_id();

    int ret = collect_leaves_rec( ntk.make_signal( n ), leaves, true );

    /* check for constant false */
    if ( ret < 0 )
    {
      leaves.clear();
    }
  }

  int collect_leaves_rec( signal<NtkSource> const& f, std::vector<signal<NtkSource>>& leaves, bool is_root )
  {
    node n = ntk.get_node( f );

    /* check if already visited */
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      for ( signal<NtkSource> const& s : leaves )
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

  void initialize_multi_aig_network( NtkDest& dest, collapsed_map& old2new )
  {
    old2new[ntk.get_node( ntk.get_constant( false ) )] = dest.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
      old2new[ntk.get_node( ntk.get_constant( true ) )] = dest.get_constant( true );

    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = dest.create_pi();
    } );

    if constexpr ( has_foreach_ro_v<NtkSource> )
    {
      static_assert( has_create_ro_v<NtkDest>, "NtkDest does not implement the create_ro method" );
      ntk.foreach_ro( [&]( auto const& n ) {
        old2new[n] = dest.create_ro();
      } );
    }
  }

private:
  NtkSource const& ntk;
  aig_collapse_params const& ps;

  storage_t storage;
  std::vector<signal<NtkDest>> children;
};

} /* namespace detail */

/*! \brief AIG collapse.
 *
 * This method collapses AND2 nodes in an AIG into
 * multi-input ANDs. The maximum number of inputs can
 * be limited using the parameter `collapse_limit`.
 * It returns the resulted network as a NtkDest.
 * 
 * NtkDest is by default a `multi_aig`. In case of
 * sequential networks, it should be defined as
 * `sequential<multi_aig>`.
 *
 * **Required network functions:**
 * - `get_node`
 * - `node_to_index`
 * - `get_constant`
 * - `create_pi`
 * - `create_po`
 * - `create_ro`
 * - `create_ri`
 * - `create_not`
 * - `is_complemented`
 * - `foreach_po`
 * - `foreach_ri`
 * - `create_nary_and`
 * - `is_ci`
 * - `is_constant`
 * - `has_and`
 */
template<class NtkSource, class NtkDest = multi_aig_network>
NtkDest aig_collapse( NtkSource const& ntk, aig_collapse_params const& ps = {} )
{
  static_assert( is_network_type_v<NtkSource>, "NtkSource is not a network type" );
  static_assert( has_get_node_v<NtkSource>, "NtkSource does not implement the get_node method" );
  static_assert( has_node_to_index_v<NtkSource>, "NtkSource does not implement the node_to_index method" );
  static_assert( has_get_constant_v<NtkSource>, "NtkSource does not implement the get_constant method" );
  static_assert( has_foreach_node_v<NtkSource>, "NtkSource does not implement the foreach_node method" );
  static_assert( has_foreach_pi_v<NtkSource>, "NtkSource does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<NtkSource>, "NtkSource does not implement the foreach_po method" );
  static_assert( has_is_pi_v<NtkSource>, "NtkSource does not implement the is_pi method" );
  static_assert( has_is_constant_v<NtkSource>, "NtkSource does not implement the is_constant method" );
  static_assert( has_clone_node_v<NtkSource>, "NtkSource does not implement the clone_node method" );
  static_assert( has_create_pi_v<NtkSource>, "NtkSource does not implement the create_pi method" );
  static_assert( has_create_po_v<NtkSource>, "NtkSource does not implement the create_po method" );
  static_assert( has_create_not_v<NtkSource>, "NtkSource does not implement the create_not method" );
  static_assert( has_is_complemented_v<NtkSource>, "NtkSource does not implement the is_complemented method" );
  static_assert( has_create_nary_and_v<NtkDest>, "NtkDest does not implement the is_complemented method" );
  static_assert( has_create_pi_v<NtkDest>, "NtkDest does not implement the create_pi method" );
  static_assert( has_create_po_v<NtkDest>, "NtkDest does not implement the create_po method" );

  detail::aig_collapse_impl<NtkSource, NtkDest> p( ntk, ps );
  NtkDest res = p.run();

  return res;
}

} /* namespace mockturtle */
