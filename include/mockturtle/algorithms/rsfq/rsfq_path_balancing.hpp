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
  \file rsfq_path_balancing.hpp
  \brief Path balancing utils for superconducting electronics

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

namespace detail
{

template<class Ntk>
class rsfq_path_balancing_impl
{
public:
  using signal = typename Ntk::signal;
  using buffer_map = node_map<std::vector<signal>, Ntk>;

public:
  explicit rsfq_path_balancing_impl( Ntk const& ntk ) : _ntk( ntk )
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

    uint32_t depth = 0;

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
          for ( auto j = buffers.size(); j <= level; ++j )
          {
            /* create a buffer */
            auto const buf = create_buffer( res_d, buffers[j - 1] );
            buffers.push_back( buf );
            res.add_binding( res.get_node( buf ), buf_id );
            res.set_as_latch( res.get_node( buf ) );
          }
        }

        children[i] = buffers[level];
      } );

      auto const new_node = res_d.clone_node( _ntk, n, children );
      old2new[n].push_back( new_node );

      depth = std::max( depth, res_d.level( res.get_node( new_node) ) );

      res.add_binding( res.get_node( new_node ), _ntk.get_binding_index( n ) );
    } );

    /* buffer POs based on circuit depth */
    _ntk.foreach_po( [&]( auto const& f ) {
      /* don't buffer constant POs */
      if ( _ntk.is_constant( _ntk.get_node( f ) ) )
      {
        res_d.create_po( old2new[f][0] );
        return true;
      }

      uint32_t level = depth - res_d.level( old2new[f][0] );
      auto& buffers = old2new[f];
      /* create buffers to pad up to the depth */
      auto i = 0u;
      for ( i = buffers.size(); i <= level; ++i )
      {
        /* create a buffer */
        auto const buf = create_buffer( res_d, buffers[i - 1] );
        buffers.push_back( buf );
        res.add_binding( res.get_node( buf ), buf_id );
        res.set_as_latch( res.get_node( buf ) );
      }
      res_d.create_po( old2new[f][i - 1] );
      return true;
    } );

    res_d.set_depth( depth );

    assert( check_balancing( res_d ) );
  }

  inline signal create_buffer( depth_view<Ntk>& res_d, signal const& fanin )
  {
    return res_d._create_node( { fanin }, 0x2 );
  }

  bool check_balancing( depth_view<Ntk>& res_d )
  {
    /* check fanin balancing */
    bool correct = true;
    res_d.foreach_gate( [&]( auto const& n ) {
      res_d.foreach_fanin( n, [&]( auto const& f ) {
        if ( res_d.level( res_d.get_node( f ) ) != res_d.level( n ) - 1 )
        {
          correct = false;
        }
        return correct;
      } );
      return correct;
    } );

    /* check balanced POs */
    if ( correct )
    {
      auto depth = res_d.depth();
      res_d.foreach_po( [&]( auto const& f ) {
        if ( res_d.is_constant( res_d.get_node( f ) ) )
          return true;
        if ( res_d.level( res_d.get_node( f ) ) != depth )
        {
          correct = false;
        }
        return correct;
      } );
    }

    return correct;
  }

private:
  uint32_t buf_id{0};
  Ntk const& _ntk;
};

} // namespace detail

/*! \brief Path balancing for RSFQ.
 *
 * This function does path balancing according
 * to the RSFQ technology constraints:
 * - Inserts padding buffers (DFF) to balance nodes' fanin
 * - Insert padding buffers (DFF) to balance POs
 *
 * \param ntk Mapped network
 */
template<class Ntk>
Ntk rsfq_path_balancing( Ntk const& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_add_binding_v<Ntk>, "Ntk does not implement the add_binding method" );
  static_assert( has_set_as_latch_v<Ntk>, "Ntk does not implement the set_as_latch method" );

  detail::rsfq_path_balancing_impl p( ntk );
  return p.run();
}

/*! \brief Check path balancing for RSFQ.
 *
 * This function checks path balancing according
 * to the RSFQ technology constraints:
 * - Checks nodes are balanced
 * - Checks POs are balanced
 *
 * \param ntk Network
 */
template<class Ntk>
bool check_buffering( Ntk const& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );

  depth_view<Ntk> ntk_d{ ntk };
  bool result = true;

  /* check path balancing */
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk_d.level( ntk.get_node( f ) ) != ntk_d.level( n ) - 1 )
        result = false;
      return result;
    } );
    return result;
  } );

  /* check balanced POs */
  if ( result )
  {
    auto depth = ntk_d.depth();
    ntk_d.foreach_po( [&]( auto const& f ) {
      /* don't check constant POs */
      if ( ntk.is_constant( ntk.get_node( f ) ) )
        return true;
      if ( ntk_d.level( ntk_d.get_node( f ) ) != depth )
      {
        result = false;
      }
      return result;
    } );
  }

  return result;
}

} // namespace mockturtle