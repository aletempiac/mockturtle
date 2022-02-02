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
  \file aqfp_depth_optimization.hpp
  \brief AQFP depth optimization

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>
#include <random>

#include "aqfp_assumptions.hpp"
#include "../../networks/buffered.hpp"
#include "../../networks/generic.hpp"
#include "../../utils/node_map.hpp"
#include "../../views/topo_view.hpp"
#include "../../views/depth_view.hpp"
#include "../../views/fanout_view.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
class aqfp_optimize_depth_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit aqfp_optimize_depth_impl( Ntk& ntk )
    : _ntk( ntk )
    , _topo_order()
  {}

public:
  Ntk run()
  {
    _topo_order.reserve( _ntk.size() );
    topo_view<Ntk>( _ntk ).foreach_node( [&]( auto const& n ) {
      _topo_order.push_back( n );
    } );

    /* get real depth */
    auto achievable_depth = depth_view<Ntk, aqfp_depth_cost>( _ntk ).depth();
    auto current_depth = depth_view<Ntk>( _ntk ).depth();

    _ntk.clear_values();
    bool success = true;
    if ( achievable_depth < current_depth )
    {
      fanout_view<Ntk> f_ntk{ _ntk };
      success = run_cut_based_depth_reduction( f_ntk, current_depth - achievable_depth );
    }

    std::cout << achievable_depth << " " << success << "\n";

    /* create resulting network */
    node_map<signal, Ntk> old2new( _ntk );
    Ntk res;

    create_res_net( res, old2new );

    // if ( !success )
    // {
    //   /* find a new configuration of buffer splitters to reduce depth */
    //   run_critical_depth_reduction( res );
    //   return res;
    // }
    // else
    // {
      return res;
    // }
  }

private:
  bool run_cut_based_depth_reduction( fanout_view<Ntk>& f_ntk, uint32_t iterations )
  {
    /* find a cut of buffers and mark them as removable */
    uint32_t i;
    for ( i = 1; i <= iterations; ++i )
    {
      _ntk.incr_trav_id();
      uint32_t trav_id = _ntk.trav_id();

      _ntk.set_visited( _ntk.get_node( _ntk.get_constant( false ) ), trav_id );

      /* mark nodes to define a cut */
      _ntk.foreach_pi( [&]( auto const& n ) {
        mark_cut_rec( f_ntk, n );
      } );

      /* extract a cut if it exist */
       _ntk.incr_trav_id();
      bool legal_cut = true;
      _ntk.foreach_po( [&]( auto const& f ) {
        legal_cut = select_buf_cut_rec( _ntk.get_node( f ), i );
        return legal_cut;
      } );

      if ( !legal_cut )
      {
        /* depth reduction is not a cut, undo last iteration and exit */
        _ntk.foreach_node( [&]( auto const& n ) {
          if ( _ntk.value( n ) == i )
            _ntk.set_value( n, 0 );
        } );
        break;
      }

      assert( check_cut() );
      assert( check_balancing() );
    }

    std::cout << iterations << std::endl;

    /* last iteration is not a cut */
    if ( i != iterations + 1 )
    {
      return false;
    }

    return true;
  }

  void mark_cut_rec( fanout_view<Ntk>& f_ntk, node const& n )
  {
    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return;

    _ntk.set_visited( n, _ntk.trav_id() );

    /* recur towards TFI */
    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( _ntk.visited( _ntk.get_node( f ) ) != _ntk.trav_id() )
      {
        mark_cut_rec( f_ntk, _ntk.get_node( f ) );
      }
    } );

    /* found a new possible buffer cut */
    if ( _ntk.is_buf( n ) && _ntk.fanout_size( n ) == 1 && !_ntk.value( n ) )
      return;

    /* recur towards TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( _ntk.visited( f ) != _ntk.trav_id() )
      {
        mark_cut_rec( f_ntk, f );
      }
    } );
  }

  bool select_buf_cut_rec( node const& n, uint32_t value )
  {
    if ( _ntk.is_constant( n ) )
      return true;

    if ( _ntk.visited( n ) == _ntk.trav_id() )
      return true;

    /* if selected buffer, set as removable */
    if ( _ntk.visited( n ) == _ntk.trav_id() - 1 && _ntk.is_buf( n ) && _ntk.fanout_size( n ) == 1 )
    {
      _ntk.set_visited( n, _ntk.trav_id() );
      /* already selected in the past iterations */
      if ( _ntk.value( n ) != 0 && _ntk.value( n ) != value )
        return false;

      _ntk.set_value( n, value );
      return true;
    }

    /* check not a cut */
    if ( _ntk.visited( n ) == _ntk.trav_id() - 1 )
    {
      _ntk.set_visited( n, _ntk.trav_id() );
      return false;
    }

    _ntk.set_visited( n, _ntk.trav_id() );

    bool legal = true;
    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      legal = select_buf_cut_rec( _ntk.get_node( f ), value );
      return legal;
    } );

    return legal;
  }

  void run_critical_depth_reduction( Ntk& ntk )
  {
    fanout_view<Ntk> f_ntk{ ntk };

    ntk.incr_trav_id();
    uint32_t trav_id = ntk.trav_id();

    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), trav_id );

    /* mark nodes to define a cut */
    ntk.foreach_pi( [&]( auto const& n ) {
      mark_cut_rec2( f_ntk, n );
    } );

    /* extract a cut if it exist */
    bool legal_cut = true;
    ntk.foreach_po( [&]( auto const& f ) {
      legal_cut &= select_buf_cut_rec2( ntk, ntk.get_node( f ) );
    } );

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        return true;
      
      if ( ntk.value( n ) )
      {
        bool legal = true;
        f_ntk.foreach_fanout( n, [&]( auto const& f ) {
          if ( ntk.visited( f ) == trav_id )
            legal = false;
        } );

        if ( !legal )
        {
          auto fanouts = f_ntk.fanout( n );
          /* check fanout is a splitter */
          std::string fanout_type;
          if ( ntk.is_buf( fanouts[0] ) && ntk.fanout_size( fanouts[0] ) > 1 )
            fanout_type = "splitter";
          else if ( ntk.is_buf( fanouts[0] ) && ntk.fanout_size( fanouts[0] ) == 1 )
            fanout_type = "buffer";
          else
            fanout_type = "gate";
          std::cout << fanout_type << "\n";
          // for ( auto const& f : fanouts )
          // {
          //   if ( ntk.visited( f ) != trav_id )
          //     ntk.subsitute_node( )
          // }
        }
      }
      return true;
    } );
  }

  void mark_cut_rec2( fanout_view<Ntk>& f_ntk, node const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f_ntk.get_node( f ) ) != f_ntk.trav_id() )
      {
        mark_cut_rec2( f_ntk, _ntk.get_node( f ) );
      }
    } );

    /* found a new possible buffer cut */
    if ( f_ntk.is_buf( n ) && f_ntk.fanout_size( n ) == 1 && !f_ntk.value( n ) )
      return;

    /* recur towards TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f ) != f_ntk.trav_id() )
      {
        mark_cut_rec2( f_ntk, f );
      }
    } );
  }

  bool select_buf_cut_rec2( Ntk& ntk, node const& n )
  {
    if ( ntk.is_constant( n ) )
      return true;

    if ( ntk.is_pi( n ) )
      return false;

    /* if selected buffer, set as removable */
    if ( ntk.visited( n ) == ntk.trav_id() && ntk.is_buf( n ) && ntk.fanout_size( n ) == 1 )
    {
      ntk.set_value( n, 1 );
      return true;
    }

    bool legal = true;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      legal &= select_buf_cut_rec2( ntk, ntk.get_node( f ) );
    } );

    return legal;
  }

  bool check_cut()
  {
    bool correct = true;

    _ntk.foreach_po( [&]( auto const& f ) {
      correct = check_cut_rec( _ntk.get_node( f ), false );
      return correct;
    } );

    return correct;
  }

  bool check_cut_rec( node const& n, bool found )
  {
    if ( _ntk.is_constant( n ) )
      return true;

    if ( _ntk.is_pi( n ) )
      return found;

    bool correct = true;
    bool buf_in_cut = _ntk.visited( n ) == _ntk.trav_id() && _ntk.is_buf( n ) && _ntk.fanout_size( n ) == 1;
    auto value = _ntk.value( n );

    if ( ( found && _ntk.value( n ) ) || ( !found && buf_in_cut && !_ntk.value( n ) ) )
      return false;
    
    if (  _ntk.value( n )  )

    _ntk.foreach_fanin( n, [&]( auto const& f ) {
      correct = check_cut_rec( _ntk.get_node( f ),  _ntk.value( n ) | found );
      return correct;
    } );

    return correct;
  }

  bool check_balancing()
  {
    depth_view<Ntk, aqfp_depth_cost_balancing> d_ntk{ _ntk };
    bool balanced = true;

    for ( auto const& n : _topo_order )
    {
      if ( _ntk.value( n ) )
        continue;

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( !_ntk.is_constant( _ntk.get_node( f ) ) && d_ntk.level( _ntk.get_node( f ) ) != d_ntk.level( n ) - 1 )
        {
          balanced = false;
        }
      } );

      if ( !balanced )
        return false;
    }

    return true;
  }

  void create_res_net( Ntk& res, node_map<signal, Ntk>& old2new )
  {
    old2new[_ntk.get_constant( false )] = res.get_constant( false );
    if ( _ntk.get_node( _ntk.get_constant( true ) ) != _ntk.get_node( _ntk.get_constant( false ) ) )
    {
      old2new[_ntk.get_constant( true )] = res.get_constant( true );
    }
    _ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    for ( auto const n : _topo_order )
    {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        continue;

      std::vector<signal> children;

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( old2new[f] ^ _ntk.is_complemented( f ) );
      } );

      signal f;
      if ( _ntk.is_buf( n ) )
      {
        if ( !_ntk.value( n ) )
        {
          /* keep */
          f = res.create_buf( children[0] );
        }
        else
        {
          /* remove */
          f = children[0];
        }
      }
      else
      {
        f = res.create_maj( children[0], children[1], children[2] );
      }
      old2new[n] = f;
    }

    _ntk.foreach_po( [&]( auto const& f ) {
      if ( _ntk.is_complemented( f ) )
        res.create_po( res.create_not( old2new[f] ) );
      else
        res.create_po( old2new[f] );
    } );
  }

  struct aqfp_depth_cost
  {
    uint32_t operator()( Ntk const& ntk, node const& node ) const
    {
      if ( ntk.is_buf( node ) && ntk.fanout_size( node ) == 1 )
        return 0u;
      else
        return 1u;
    }
  };

  struct aqfp_depth_cost_balancing
  {
    uint32_t operator()( Ntk const& ntk, node const& node ) const
    {
      if ( ntk.is_buf( node ) && ntk.fanout_size( node ) == 1 && ntk.value( node ) )
        return 0u;
      else
        return 1u;
    }
  };

private:
  Ntk &_ntk;
  std::vector<node> _topo_order;
};

} /* namespace detail */

/*! \brief Depth optimization for AQFP networks.
 *
 * This function tries to reduce the depth of a mapped AQFP circuit
 *
 * \param ntk Mapped AQFP network
 */
template<class Ntk>
Ntk aqfp_optimize_depth( Ntk& ntk )
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
  static_assert( is_buffered_network_type_v<Ntk>, "BufNtk is not a buffered network type" );
  static_assert( has_is_buf_v<Ntk>, "BufNtk does not implement the is_buf method" );
  static_assert( has_create_buf_v<Ntk>, "BufNtk does not implement the create_buf method" );

  detail::aqfp_optimize_depth_impl p( ntk );
  return p.run();
}

} // namespace mockturtle