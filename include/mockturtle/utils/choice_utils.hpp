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

/*!
  \file topo_view.hpp
  \brief Reimplements foreach_node to guarantee topological order

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cassert>
#include <vector>
#include <iostream>

#include "../algorithms/cleanup.hpp"
#include "../algorithms/functional_reduction.hpp"
#include "cost_functions.hpp"
#include "../networks/detail/foreach.hpp"
#include "../utils/node_map.hpp"
#include "../views/depth_choice_view.hpp"
#include "../views/choice_view.hpp"
#include "../views/topo_view.hpp"
#include "../traits.hpp"

namespace mockturtle::detail
{

template<typename Ntk>
bool check_choice_in_tfi_rec( choice_view<Ntk> const& ntk, node<Ntk> const& n, node<Ntk> const& choice )
{
  ntk.set_visited( n, ntk.trav_id() );

  if ( ntk.is_ci( n ) )
    return false;

  bool found = false;
  ntk.foreach_fanin( n, [&]( auto const& f ) {
    if ( ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() )
    {
      return true;
    }
    if ( ntk.get_choice_repr( ntk.get_node( f ) ) == choice )
    {
      found = true;
      return false;
    }

    found = check_choice_in_tfi_rec( ntk, ntk.get_node( f ), choice );
    return !found;
  } );
  return found;
}


/*! \brief Remove choices in TFI.
 *
 * Remove choices that contain the class representative in the
 * transitive fainin cone.
 */
template<typename Ntk>
void remove_choices_in_tfi( Ntk& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( !ntk.is_choice_repr( n ) ) {
      ntk.incr_trav_id();
      auto const& repr = ntk.get_choice_repr( n );
      if ( check_choice_in_tfi_rec( ntk, n, repr ) )
      {
        ntk.remove_choice( n );
      }
    }
  } );
}


/*! \brief Choice network consistency.
 *
 * Checks that the choice network is correctly structured.
 */
template<typename Ntk>
void check_consistency( choice_view<Ntk> const& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  ntk.foreach_node( [&] ( const auto& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( !ntk.is_choice_repr( n ) )
    {
      assert( ntk.fanout_size( n ) == 0u );
      assert( ntk.is_choice( n ) == true );
    }
    else
    {
      if ( ntk.is_choice( n ) )
      {
        assert( ntk.fanout_size( n ) == 0u );
      }
      else
      {
        assert( ntk.fanout_size( n ) != 0u );
      }
    }
  } );
}


/*! \brief Replace choices by the class representative.
 *
 * Each node in an equivalence class is substituted by the class
 * representative in the network.
 */
template<typename Ntk>
void replace_choices_by_repr( Ntk& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( !ntk.is_choice_repr( n ) && !ntk.is_choice( n ) ) {
      auto g = ntk.get_choice_repr_signal( n );
      ntk.substitute_node( n, g );
    }
  } );
}


/*! \brief Recursive deferencing for choice networks
 *
 * Recursive deferencing on the class representative
 */
template<typename Ntk, class NodeCostFn = unit_cost<Ntk>>
uint32_t choice_recursive_deref( Ntk const& ntk, node<Ntk> const& n )
{
  /* terminate? */
  if ( ntk.is_ci( n ) )
    return 0;

  /* recursively collect nodes */
  uint32_t value = NodeCostFn{}( ntk, n );
  ntk.foreach_fanin( n, [&]( auto const& child ) {
    auto s = ntk.get_choice_repr( ntk.get_node( child ) );
    if ( ntk.decr_value( s ) == 0 )
    {
      value += choice_recursive_deref<Ntk>( ntk, s );
    }
  } );
  return value;
}


/*! \brief Recursive referencing for choice networks
 *
 * Recursive referencing on the class representative
 */
template<typename Ntk, class NodeCostFn = unit_cost<Ntk>>
uint32_t choice_recursive_ref( Ntk const& ntk, node<Ntk> const& n )
{
  /* terminate? */
  if ( ntk.is_ci( n ) )
    return 0;

  /* recursively collect nodes */
  uint32_t value = NodeCostFn{}( ntk, n );
  ntk.foreach_fanin( n, [&]( auto const& child ) {
    auto s = ntk.get_choice_repr( ntk.get_node( child ) );
    if ( ntk.incr_value( s ) == 0 )
    {
      value += choice_recursive_ref<Ntk>( ntk, s );
    }
  } );
  return value;
}


/*! \brief Compute the required time in a choice network
 *
 * Compute the required time in a choice network given a max
 * required depth value
 */
template<typename Ntk, class DepthCostFn = unit_cost<Ntk>>
std::vector<int32_t> compute_required( Ntk const& ntk, uint32_t depth )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  std::vector<int32_t> required( ntk.size(), depth );

  auto i = ntk.size();
  while( i-- > 0 )
  {
    auto n = ntk.index_to_node( i );
    if ( ntk.is_ci( n ) )
      continue;
    // if ( ntk.is_choice_repr( n ) )
    if ( ntk.value( n ) && ntk.is_choice_repr( n ) )
    {
      ntk.foreach_fanin( n, [&]( auto const& child ) {
        auto child_repr = ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( child ) ) );
        int32_t cost = static_cast<int32_t>( DepthCostFn{}( ntk, child_repr ) );
        required[child_repr] = std::min( required[child_repr], required[i] - cost );
      } );
      ntk.foreach_choice( n, [&]( auto const& c ) {
        required[ntk.node_to_index( c )] = required[i];
        return true;
      } );
    }
  }
  return required;
}


template<typename Ntk, class NtkDest = Ntk>
void levelize_choice_network_rec( node<Ntk> const& root, choice_view<Ntk> const& src, choice_view<NtkDest>& dest, node_map<signal<NtkDest>, Ntk>& old_to_new )
{
  /* is permanently marked? */
  if ( src.visited( root ) == src.trav_id() )
    return;

  assert( src.is_choice_repr( root ) );

  src.foreach_choice( root, [&]( auto const& n ) {
    /* is permanently marked? */
    if ( src.visited( n ) == src.trav_id() )
      return true;

    /* ensure that the node is not temporarily marked */
    assert( src.visited( n ) != src.trav_id() - 1 );

    /* mark node temporarily */
    src.set_visited( n, src.trav_id() - 1 );

    /* mark children */
    src.foreach_fanin( n, [&]( auto const& child ) {
      auto const& repr = src.get_choice_repr( src.get_node( child ) );
      levelize_choice_network_rec( repr, src, dest, old_to_new );
    } );
    return true;
  } );

  node<Ntk> new_repr = 0;

  src.foreach_choice( root, [&]( auto const& n ) {
    /* mark node n permanently */
    if ( src.visited( n ) == src.trav_id() )
      return true;

    src.set_visited( n, src.trav_id() );

    /* collect children from equivalent class representative */
    std::vector<signal<Ntk>> children;
    src.foreach_fanin( n, [&]( auto child, auto ) {
      // assert( src.is_choice_repr( src.get_node( child ) ) );
      auto const& repr = src.get_choice_repr_signal( child );
      const auto f = old_to_new[repr];

      if ( src.is_complemented( repr ) )
      {
        children.push_back( dest.create_not( f ) );
      }
      else
      {
        children.push_back( f );
      }
    } );

    auto new_sig = dest.clone_node( src, n, children );
    old_to_new[n] = new_sig;

    if ( n != root )
    {
      dest.add_choice( new_repr, new_sig );
    }
    else
    {
      new_repr = dest.get_node( new_sig );
    }
    return true;
  } );
}

/*! \brief Compute area flow in a choice network
 *
 * Area flow, node based, for a choice network
 */
template<typename Ntk>
inline float area_flow( Ntk _ntk, node<Ntk> const& _root, std::vector<float> const& _area )
{
  float _res = 1.0;
  _ntk.foreach_fanin( _root, [&]( const auto sig ) {
    auto _node = _ntk.get_node( sig );
    float fanout_size = ( float ) _ntk.value( _node );

    if ( fanout_size == 0.0 )
      fanout_size = 1.0;

    _res += _area.at( _ntk.node_to_index( _node ) ) / fanout_size;
  } );

  return _res;
}

/*! \brief Initialize values with fanout size
 *
 * Initialize values in equivalence classes based on
 * the fanout size of the representative
 */
template<typename Ntk>
void init_value_with_fanout( Ntk& ntk )
{
  ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_choice_repr( n ) ) {
        ntk.set_value( n, ntk.fanout_size( n ) );
      } else {
        ntk.set_value( n, ntk.fanout_size( ntk.get_choice_repr( n ) ) );
      }
    } );
}

/*! \brief Update values w.r.t representatives
 *
 * Update values in classes w.r.t choice representatives
 */
template<typename Ntk>
void update_value_with_repr( Ntk& ntk )
{
  ntk.foreach_node( [&]( auto const& n ) {
    if ( !ntk.is_choice_repr( n ) )
      ntk.set_value( n, ntk.value( ntk.get_choice_repr( n ) ) );
  } );
}

} /* namespace mockturtle::detail */

namespace mockturtle
{

/*! \brief Add equivalences pairs as choices.
 *
 * Given the equivalence pairs, the function add the nodes
 * as choices.
 */
template<typename Ntk>
void insert_equivalences( choice_view<Ntk>& ntk, std::vector<std::pair<node<Ntk>, signal<Ntk>>> const& equivalences )
{  
  for ( auto const& pair : equivalences )
  {
    /* filter out possible dead nodes */
    if ( !( ntk.is_dead( std::get<0>( pair ) ) || ntk.is_dead( ntk.get_node( std::get<1>( pair ) ) ) ) )
      ntk.add_choice( std::get<0>( pair ), std::get<1>( pair ) );
  }
}


/*! \brief Reduce choice network given node equivalences
 *
 * Reduce a choice network given node equivalences:
 * replaces choices nodes by the class representative
 */
template<typename Ntk>
void reduce_choice_network( choice_view<Ntk>& ntk, std::vector<std::pair<node<Ntk>, signal<Ntk>>> const& equivalences )
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

  insert_equivalences( ntk, equivalences );
  detail::replace_choices_by_repr( ntk );
  detail::remove_choices_in_tfi( ntk );
}


// template<typename Ntk>
// choice_view<Ntk> cleanup_choice_network( Ntk& src, choice_view<Ntk> const& choice_src, node_map<signal<Ntk>, choice_view<Ntk> choice_map )


/*! \brief Update representatives using the choices currently used in the network.
 *
 * The representative in each class is updated with the choice currently
 * in use in the network. If more choices are in use, the one with the
 * highest index is used (last nodes added).
 */
template<typename Ntk>
void update_representatives( Ntk& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( !ntk.is_choice_repr( n ) && !ntk.is_choice( n ) )
    {
      ntk.update_choice_repr( n );
    }
  } );
  detail::replace_choices_by_repr( ntk );
  detail::remove_choices_in_tfi( ntk );
}


/*! \brief Improves the representative with a depth optimization strategy
 *
 * Improves the representative with a depth optimization strategy:
 * - depth optimization
 * - area recovery
 */
template<typename Ntk, class DepthCostFn = unit_cost<Ntk>, class NodeCostFn = unit_cost<Ntk>>
void improve_representatives( choice_view<Ntk>& ntk )
{
  std::vector<uint32_t> arrival( ntk.size(), 0 );
  uint32_t depth = 0u;

  detail::init_value_with_fanout( ntk );

  /* improve level */
  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    ntk.foreach_fanin( n, [&]( auto const& child ) {
      auto child_repr = ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( child ) ) );
      arrival[ntk.node_to_index( n )] = std::max( arrival[ntk.node_to_index( n )], arrival[child_repr] + DepthCostFn{}( ntk, child_repr ) );
    } );
    if ( ntk.is_choice_repr( n ) )
    {
      uint32_t min_level = UINT32_MAX;
      uint32_t min_mffc = UINT32_MAX;
      auto repr = n;

      if ( ntk.value( n ) )
        detail::choice_recursive_deref<choice_view<Ntk>, NodeCostFn>( ntk, n );

      ntk.foreach_choice( n, [&]( auto const& g ) {
        if ( arrival[ntk.node_to_index( g )] == 0u )
        {
          ntk.foreach_fanin( g, [&]( auto const& child ) {
            auto child_repr = ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( child ) ) );
            arrival[ntk.node_to_index( g )] = std::max( arrival[ntk.node_to_index( g )], arrival[child_repr] + DepthCostFn{}( ntk, child_repr ) );
          } );
        }
        auto level = arrival[ntk.node_to_index( g )];
        auto mffc = detail::choice_recursive_ref<choice_view<Ntk>, NodeCostFn>( ntk, g );
        auto v2 = detail::choice_recursive_deref<choice_view<Ntk>, NodeCostFn>( ntk, g );
        assert( mffc == v2 );
        if ( level < min_level || ( level == min_level && mffc < min_mffc ) )
        {
          min_level = level;
          min_mffc = mffc;
          repr = g;
        }
        return true;
      } );

      if ( ntk.value( n ) )
        detail::choice_recursive_ref<choice_view<Ntk>, NodeCostFn>( ntk, repr );

      ntk.update_choice_repr( repr );
    }
  } );

  ntk.foreach_po( [&]( auto const& po ) {
    depth = std::max( arrival[ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( po ) ) )], depth );
  } );

  for ( auto i = 0u; i < ntk.size(); arrival[i++] = 0u );

  detail::update_value_with_repr( ntk );
  auto const required = detail::compute_required<choice_view<Ntk>, DepthCostFn>( ntk, depth );

  /* area recovery */
  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    ntk.foreach_fanin( n, [&]( auto const& child ) {
      auto child_repr = ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( child ) ) );
      arrival[ntk.node_to_index( n )] = std::max( arrival[ntk.node_to_index( n )], arrival[child_repr] + DepthCostFn{}( ntk, child_repr ) );
    } );
    if ( ntk.is_choice_repr( n ) )
    {
      uint32_t min_mffc = UINT32_MAX;
      auto repr = n;

      if ( ntk.value( n ) )
        detail::choice_recursive_deref<choice_view<Ntk>, NodeCostFn>( ntk, n );

      ntk.foreach_choice( n, [&]( auto const& g ) {
        if ( arrival[ntk.node_to_index( g )] == 0u )
        {
          ntk.foreach_fanin( g, [&]( auto const& child ) {
            auto child_repr = ntk.node_to_index( ntk.get_choice_repr( ntk.get_node( child ) ) );
            arrival[ntk.node_to_index( g )] = std::max( arrival[ntk.node_to_index( g )], arrival[child_repr] + DepthCostFn{}( ntk, child_repr ) );
          } );
        }

        if ( required[ntk.node_to_index( g )] >= (int) arrival[ntk.node_to_index( g )] )
        {
          auto mffc = detail::choice_recursive_ref<choice_view<Ntk>, NodeCostFn>( ntk, g );
          auto v2 = detail::choice_recursive_deref<choice_view<Ntk>, NodeCostFn>( ntk, g );
          assert( mffc == v2 );
          if ( mffc < min_mffc )
          {
            min_mffc = mffc;
            repr = g;
          }
        }
        return true;
      } );

      if ( ntk.value( n ) )
        detail::choice_recursive_ref<choice_view<Ntk>, NodeCostFn>( ntk, repr );

      ntk.update_choice_repr( repr );
    }
  } );

  detail::replace_choices_by_repr( ntk );
}


/*! \brief Improves the representative using an area optimization strategy
 *
 * Improves the representative using an area optimization strategy
 */
template<typename Ntk>
void improve_representatives_area( choice_view<Ntk>& ntk )
{
  // std::vector<float> best_area( ntk.size() );
  ntk.clear_values();

  detail::init_value_with_fanout( ntk );

  ntk.foreach_node( [&]( auto const& n ) {
    // auto node_index = ntk.node_to_index( n );
    if ( ntk.is_ci( n ) )
    {
      // best_area[node_index] = 0u;
    }
    else if ( ntk.is_choice_repr( n ) )
    {
      // float min_aflow = std::numeric_limits<float>::max();
      uint32_t min_mffc = UINT32_MAX;
      auto repr = n;

      if ( ntk.value( n ) )
        detail::choice_recursive_deref( ntk, n );

      ntk.foreach_choice( n, [&]( auto const& g ) {
        // auto aflow = detail::area_flow( ntk, g, best_area );
        auto mffc = detail::choice_recursive_ref( ntk, g );
        auto v2 = detail::choice_recursive_deref( ntk, g );
        assert( mffc == v2 );
        // if ( aflow < min_aflow || ( aflow == min_aflow && mffc < min_mffc ) )
        if ( mffc < min_mffc )
        {
          // min_aflow = aflow;
          min_mffc = mffc;
          repr = g;
        }
        return true;
      } );
      // best_area[node_index] = min_aflow;
      if ( ntk.value( n ) )
        detail::choice_recursive_ref( ntk, repr );

      ntk.update_choice_repr( repr );
    }
  } );

  detail::replace_choices_by_repr( ntk );
}


/*! \brief Levelize and cleanup a choice network
 *
 * Rebuilds a choice network in a topological order levelizing
 * equivalence nodes. All the nodes in the equivalence class are
 * stored in the indices following the representative.
 * It cleans the dead nodes
 */
template<typename Ntk, class NtkDest = Ntk>
choice_view<NtkDest> levelize_choice_network( choice_view<Ntk> const& src )
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

  NtkDest dest;

  node_map<signal<NtkDest>, Ntk> old_to_new( src );

  src.incr_trav_id();
  src.incr_trav_id();
  
  const auto c0 = src.get_node( src.get_constant( false ) );
  src.set_visited( c0, src.trav_id() );

  old_to_new[src.get_constant( false )] = dest.get_constant( false );

  if ( src.get_node( src.get_constant( true ) ) != src.get_node( src.get_constant( false ) ) )
  {
    old_to_new[src.get_constant( true )] = dest.get_constant( true );
    src.set_visited( src.get_node( src.get_constant( true ) ), src.trav_id() );
  }

  src.foreach_pi( [&]( auto const& n ) {
    old_to_new[n] = dest.create_pi();
  } );

  src.foreach_ci( [&]( auto n ) {
    if ( src.visited( n ) != src.trav_id() )
    {
      src.set_visited( n, src.trav_id() );
    }
  } );

  choice_view<NtkDest> choice_net{dest};

  src.foreach_po( [&]( auto f ) {
    if ( src.visited( src.get_node( f ) ) == src.trav_id() )
      return;

    detail::levelize_choice_network_rec( src.get_node( f ), src, choice_net, old_to_new );
  } );

  /* create outputs in same order */
  src.foreach_po( [&]( auto po ) {
    const auto f = old_to_new[po];
    if ( src.is_complemented( po ) )
    {
      dest.create_po( dest.create_not( f ) );
    }
    else
    {
      dest.create_po( f );
    }
  } );

  choice_net.foreach_node( [&]( auto const& n ) {
    if ( !choice_net.is_choice_repr( n ) )
    {
      choice_net.take_out_choice( n );
    }
  } );

  return choice_net;
}


/*! \brief Creates a choice network starting from two equivalent networks
 *
 * Creates a choice networks starting from src1, it adds src2 nodes and
 * runs functional_reduction to find equivalent nodes.
 * A final levelized choice network is returned
 */
template<typename Ntk>
choice_view<Ntk> create_choice_network( Ntk const& src1, Ntk const& src2 )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_pi_at_v<Ntk>, "Ntk does not implement the pi_at method" );
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

  assert( src1.num_pis() == src2.num_pis() && src1.num_pos() == src2.num_pos() );

  Ntk dest = cleanup_dangling( src1 );

  node_map<signal<Ntk>, Ntk> old_to_new( src2 );

  old_to_new[src2.get_constant( false )] = dest.get_constant( false );

  if ( src2.get_node( src2.get_constant( true ) ) != src2.get_node( src2.get_constant( false ) ) )
  {
    old_to_new[src2.get_constant( true )] = dest.get_constant( true );
  }

  src2.foreach_pi( [&]( auto const& n, auto i ) {
    old_to_new[n] = dest.make_signal( dest.pi_at( i ) );
  } );

  topo_view topo{src2};
  topo.foreach_node( [&]( auto const& n ) {
    if ( src2.is_constant( n ) || src2.is_pi( n ) )
      return;

    std::vector<signal<Ntk>> children;
    src2.foreach_fanin( n, [&]( auto child, auto ) {
      const auto f = old_to_new[child];

      if ( src2.is_complemented( child ) )
      {
        children.push_back( dest.create_not( f ) );
      }
      else
      {
        children.push_back( f );
      }
    } );
    old_to_new[n] = dest.clone_node( src2, n, children );
  } );

  functional_reduction_params ps;
  functional_reduction_stats st;
  ps.compute_equivalence_classes = true;
  auto eqpairs = functional_reduction_eqclasses( dest, ps, &st );

  choice_view<Ntk> choice_dest{dest};

  /* create outputs in same order */
  src2.foreach_po( [&]( auto const& po, auto i ) {
    const auto f = old_to_new[po];
    bool inv = dest.is_complemented( dest.po_at( i ) ) ^ src2.is_complemented( po );
    choice_dest.add_choice( dest.get_node( dest.po_at( i ) ), f ^ inv );
  } );

  reduce_choice_network( choice_dest, eqpairs );
  
  return levelize_choice_network( choice_dest );
}

} // namespace mockturtle
