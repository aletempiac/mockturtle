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

#include "../networks/detail/foreach.hpp"
#include "../utils/node_map.hpp"
#include "../views/depth_choice_view.hpp"
#include "../views/choice_view.hpp"
#include "../algorithms/detail/mffc_utils.hpp"
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
    if ( ntk.get_node( f ) == choice )
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
void remove_choices_in_tfi( choice_view<Ntk>& ntk )
{
  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( !ntk.is_choice_repr( n ) ) {
      ntk.incr_trav_id();
      ntk.set_visited( n, ntk.trav_id() );
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
void replace_choices_by_repr( choice_view<Ntk>& ntk )
{
  // bool convergency = false;
  // while( !convergency )
  // {
    // convergency = true;
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_ci( n ) )
        return;
      if ( !ntk.is_choice_repr( n ) && !ntk.is_choice( n ) ) {
        auto g = ntk.get_choice_repr_signal( n );
        ntk.substitute_node( n, g );
        // convergency = false;
      }
    } );
  // }
  // check_consistency( ntk );
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


template<typename Ntk>
inline float area_flow( choice_view<Ntk> _ntk, node<Ntk> const& _root, std::vector<float> const& _area )
{
  float _res = 1;
  _ntk.foreach_fanin( _root, [&]( const auto sig ) {
    auto _node = _ntk.get_node( sig );
    float fanout_size = 1.0;
    if ( _ntk.fanout_size( _node ) > 0 ) {
      fanout_size = static_cast<float>( _ntk.fanout_size( _node ) );
    } else if ( !_ntk.is_choice_repr( _node ) ) {
      fanout_size = static_cast<float>( _ntk.fanout_size( _ntk.get_choice_repr( _node ) ) );
    }

    assert( fanout_size != 0 );
    _res += _area.at( _ntk.node_to_index( _node ) ) / fanout_size;
  } );

  return _res;
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
    ntk.add_choice( std::get<0>( pair ), std::get<1>( pair ) );
  }
}


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


template<typename Ntk>
void improve_representatives( choice_view<Ntk>& ntk )
{
  depth_choice_view depth_ntk{ntk};

  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_choice_repr( n ) ) {
      ntk.set_value( n, ntk.fanout_size( n ) );
    } else {
      ntk.set_value( n, ntk.fanout_size( ntk.get_choice_repr( n ) ) );
    }
    if ( ntk.value( n ) == 0u )
      ntk.incr_value( n );
  } );

  ntk.foreach_node( [&]( auto const& n ) {
    if ( ntk.is_ci( n ) )
      return;
    if ( ntk.is_choice_repr( n ) )
    {
      uint32_t min_level = UINT32_MAX;
      uint32_t min_mffc = UINT32_MAX;
      auto repr = n;
      ntk.foreach_choice( n, [&]( auto const& g ) {
        auto level = depth_ntk.level( g );
        auto mffc = detail::mffc_size( ntk, g );
        if ( level < min_level || ( level == min_level && mffc < min_mffc ) )
        // if ( mffc < min_mffc || ( mffc == min_mffc && level < min_level ) )
        {
          min_level = level;
          min_mffc = mffc;
          repr = g;
        }
        return true;
      } );
      ntk.update_choice_repr( repr );
    }
  } );

  detail::replace_choices_by_repr( ntk );
}


template<typename Ntk>
void improve_representatives_area( choice_view<Ntk>& ntk )
{
  std::vector<float> best_area( ntk.size() );

  ntk.foreach_node( [&]( auto const& n ) {
    auto node_index = ntk.node_to_index( n );
    if ( ntk.is_ci( n ) )
    {
      best_area[node_index] = 0u;
    }
    else if ( ntk.is_choice_repr( n ) )
    {
      float min_aflow = std::numeric_limits<float>::max();
      uint32_t min_mffc = UINT32_MAX;
      auto repr = n;
      ntk.foreach_choice( n, [&]( auto const& g ) {
        auto aflow = detail::area_flow( ntk, g, best_area );
        auto mffc = detail::mffc_size( ntk, g );
        if ( aflow < min_aflow || ( aflow == min_aflow && mffc < min_mffc ) )
        {
          min_aflow = aflow;
          min_mffc = mffc;
          repr = g;
        }
        return true;
      } );
      best_area[node_index] = min_aflow;
      ntk.update_choice_repr( repr );
    }
  } );

  detail::replace_choices_by_repr( ntk );
}


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

} // namespace mockturtle
