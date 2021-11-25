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

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>

#include "../algorithms/cleanup.hpp"
#include "../networks/klut.hpp"
#include "../networks/generic.hpp"
#include "../utils/node_map.hpp"
#include "../views/binding_view.hpp"
#include "../views/topo_view.hpp"
#include "../views/depth_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class buffering_impl
{
public:
  using signal = typename Ntk::signal;
  using buffer_map = node_map<std::vector<signal>, Ntk>;

public:
  explicit buffering_impl( Ntk const& ntk ) : _ntk( ntk )
  {}

  Ntk run()
  {
    auto [res, old2new] = initialize_copy_buf_network();

    load_dff_element();

    generate_buffered_network( res, old2new );

    return res;
  }

private:
  std::pair<Ntk, buffer_map> initialize_copy_buf_network()
  {
    buffer_map old2new( _ntk );
    Ntk res( _ntk.get_library() );

    old2new[_ntk.get_constant( false )].push_back( res.get_constant( false ) );
    if ( _ntk.get_node( _ntk.get_constant( true ) ) != _ntk.get_node( _ntk.get_constant( false ) ) )
    {
      old2new[_ntk.get_constant( true )].push_back( res.get_constant( true ) );
    }
    _ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n].push_back( res.create_pi() );
    } );
    return {res, old2new};
  }

  void load_dff_element()
  {
    for ( auto const& gate : _ntk.get_library() )
    {
      if ( gate.num_vars == 1 && kitty::is_const0( kitty::cofactor0( gate.function, 0 ) ) )
      {
        buf_id = gate.id;
        return;
      }
    }
  }

  void generate_buffered_network( Ntk& res, buffer_map& old2new )
  {
    depth_view res_d{res};

    /* network is supposed to be stored in topo order */
    _ntk.foreach_gate( [&]( auto const& n ) {
      uint32_t max_level = 0;

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        max_level = std::max( max_level, res_d.level( old2new[f][0] ) );
      } );

      std::vector<signal> children( _ntk.fanin_size( n ) );

      _ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        /* level indicates the number of padding buffers */
        uint32_t level = max_level - res_d.level( old2new[f][0] );
        auto& buffers = old2new[f];

        if ( buffers.size() <= level )
        {
          /* create buffers to pad up to the level */
          for ( auto i = buffers.size(); i <= level; ++i )
          {
            /* create a buffer */
            auto const buf = create_buffer( res_d, buffers[i - 1] );
            buffers.push_back( buf );
            res.add_binding( res.get_node( buf ), buf_id );
            res.set_as_latch( res.get_node( buf ) );
          }
        }

        children[i] = buffers[level];
      } );

      auto const new_node = res_d.clone_node( _ntk, n, children );
      old2new[n].push_back( new_node );

      res.add_binding( res.get_node( new_node ), _ntk.get_binding_index( n ) );
    } );

    /* no PO balancing for now */
    _ntk.foreach_po( [&]( auto const& f ) {
      res_d.create_po( old2new[f][0] );
    } );

    check_balancing( res_d );
  }

  inline signal create_buffer( depth_view<Ntk>& res_d, signal const& fanin )
  {
    static uint64_t _buf = 0x2;
    kitty::dynamic_truth_table tt_buf( 1 );
    kitty::create_from_words( tt_buf, &_buf, &_buf + 1 );
    return res_d.create_node( { fanin }, tt_buf );
  }

  void check_balancing( depth_view<Ntk>& res_d )
  {
    bool correct = true;
    res_d.foreach_gate( [&]( auto const& n ) {
      res_d.foreach_fanin( n, [&]( auto const& f ) {
        if ( res_d.level( res_d.get_node( f ) ) != res_d.level( n ) - 1 )
        {
          correct = false;
        }
      } );
    } );

    if ( !correct )
    {
      std::cout << "You messed up again!\n";
    }
  }

private:
  uint32_t buf_id{0};
  Ntk const& _ntk;
};

}

template<class Ntk>
Ntk buffering( Ntk const& ntk )
{
  detail::buffering_impl p( ntk );
  return p.run();
}

namespace detail
{

template<class Ntk>
class generic_network_convert_impl
{
public:
  using signal = typename generic_network::signal;
  using NtkDest = binding_view<generic_network>;

public:
  explicit generic_network_convert_impl( Ntk const& ntk ) : _ntk( ntk )
  {}

  binding_view<generic_network> run()
  {
    node_map<signal, Ntk> old2new( _ntk );
    NtkDest res( _ntk.get_library() );

    old2new[_ntk.get_constant( false )] = res.get_constant( false );
    if ( _ntk.get_node( _ntk.get_constant( true ) ) != _ntk.get_node( _ntk.get_constant( false ) ) )
    {
      old2new[_ntk.get_constant( true )] = res.get_constant( true );
    }
    _ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    topo_view topo{ _ntk };

    topo.foreach_node( [&] ( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return true;

      std::vector<signal> children;
      
      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( old2new[f] );
      } );

      if ( _ntk.is_as_latch( n ) )
      {
        const auto f = res.create_latch( children[0] );
        res.foreach_fanin( res.get_node( f ), [&]( auto const& latch ) {
          res.add_binding( res.get_node( latch ), _ntk.get_binding_index( n ) );
          return false;
        } );
        old2new[n] = f;
      }
      else
      {
        const auto f = res.create_node( children, _ntk.node_function( n ) );
        res.add_binding( res.get_node( f ), _ntk.get_binding_index( n ) );
        old2new[n] = f;
      }

      return true;
    } );

    _ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( old2new[f] );
    } );

    return res;
  }

private:
  Ntk const& _ntk;
};

}

template<class Ntk>
binding_view<generic_network> generic_network_convert( Ntk const& ntk )
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

  detail::generic_network_convert_impl p( ntk );
  return p.run();
}

}