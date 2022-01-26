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
  \file buffer_network_conversion.hpp
  \brief Network conversion for generic network utils

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <vector>
#include <random>

#include "aqfp_assumptions.hpp"
#include "../../networks/klut.hpp"
#include "../../networks/buffered.hpp"
#include "../../networks/generic.hpp"
#include "../../utils/node_map.hpp"
#include "../../views/topo_view.hpp"
#include "../../views/fanout_view.hpp"

namespace mockturtle
{

struct aqfp_network_conversion_params
{
  /*! \brief AQFP technology assumptions. */
  aqfp_assumptions aqfp_assumptions_ps{};

  /*! \brief Random assignement. */
  bool use_random{ true };

  /*! \brief Random seed. */
  std::default_random_engine::result_type seed{ 1 };

  /*! \brief Direction of preferred registers. */
  bool forward{ true };
};

namespace detail
{

template<class Ntk>
class generic_network_create_from_buffered_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using signal_d = typename generic_network::signal;

public:
  explicit generic_network_create_from_buffered_impl( Ntk& ntk, aqfp_network_conversion_params const& ps )
    : _ntk( ntk )
    , _ps( ps )
  {}

public:
  generic_network run()
  {
    node_map<signal_d, Ntk> old2new( _ntk );
    generic_network res;

    old2new[_ntk.get_constant( false )] = res.get_constant( false );
    if ( _ntk.get_node( _ntk.get_constant( true ) ) != _ntk.get_node( _ntk.get_constant( false ) ) )
    {
      old2new[_ntk.get_constant( true )] = res.get_constant( true );
    }
    _ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    topo_view topo{ _ntk };

    select_retimeable_elements_random( topo );

    create_generic_network( topo, res, old2new );

    return res;
  }

private:
  void select_retimeable_elements_random( topo_view<Ntk>& topo )
  {
    fanout_view fntk{ _ntk };
    auto rseed = _ps.seed;

    _ntk.clear_values();

    /* select buffers and splitters to retime as soon as found some */
    _ntk.incr_trav_id();
    topo.foreach_node( [&] ( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return true;

      if ( _ntk.is_buf( n ) )
      {
        if ( _ntk.fanout_size( n ) == 1 )
        {
          _ntk.set_visited( n, _ntk.trav_id() );
        }
        else if ( _ntk.visited( n ) != _ntk.trav_id() ||  _ntk.value( n ) > 0 )
        {
          int free_spots;
          if ( _ntk.value( n ) > 0 )
          {
            free_spots = rec_fetch_root( n );
            if ( free_spots == 0 )
              return true;
          }
          else
          {
            free_spots = _ps.aqfp_assumptions_ps.splitter_capacity - _ntk.fanout_size( n );
          }

          int total_fanout = 0;
          std::vector<node> fanout_splitters;

          /* select retimeable splitters */
          fntk.foreach_fanout( n, [&]( auto const f ) {
            if ( _ntk.is_buf( f ) && _ntk.fanout_size( f ) > 1 && free_spots >= _ntk.fanout_size( f ) - 1 )
            {
              fanout_splitters.push_back( f );
              total_fanout += _ntk.fanout_size( f ) - 1;
            }
          } );

          /* check if they are all retimeable together */
          if ( free_spots >= total_fanout )
          {
            for ( auto f : fanout_splitters )
            {
              _ntk.set_value( f, free_spots - total_fanout );
              _ntk.set_visited( f, _ntk.trav_id() );
            }
            rec_update_root( n, free_spots - total_fanout );
            return true;
          }
          /* select one randomly */
          std::default_random_engine gen( rseed++ );
          std::uniform_int_distribution<uint32_t> dist( 0ul, fanout_splitters.size() - 1 );
          auto index = dist( gen );
          _ntk.set_value( fanout_splitters[index], free_spots - _ntk.fanout_size( fanout_splitters[index] ) + 1 );
          _ntk.set_visited( fanout_splitters[index], _ntk.trav_id() );
          rec_update_root( n, free_spots - _ntk.fanout_size( fanout_splitters[index] ) + 1 );
        }
      }
      return true;
    } );
  }

  // void select_retimeable_elements_simulate( topo_view<Ntk>& topo )
  // {
  //   fanout_view fntk{ _ntk };

  //   /* set the nodes to be retimed to visited */
  //   /* set value `index` for splitters to test */
  //   ntk.incr_trav_id();
  //   topo.foreach_node( [&] ( auto const& n ) {
  //     if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
  //       return true;

  //     if ( ntk.is_buf( n ) )
  //     {
  //       if ( ntk.fanout_size( n ) == 1 )
  //       {
  //         ntk.set_visited( n, ntk.trav_id() );
  //       }
  //       else if ( !ntk.value( n ) && ntk.fanout_size( n ) < ps.splitter_capacity )
  //       {
  //         int free_spots = ps.splitter_capacity - ntk.fanout_size( n );
  //         int total_fanout = 0;
  //         std::vector<typename Ntk::node> fanout_splitters;

  //         /* select retimeable splitters */
  //         fntk.foreach_fanout( n, [&]( auto const f ) {
  //           if ( ntk.is_buf( f ) && ntk.fanout_size( f ) > 1 && free_spots >= ntk.fanout_size( f ) - 1 )
  //           {
  //             fanout_splitters.push_back( f );
  //             total_fanout += ntk.fanout_size( f ) - 1;
  //           }
  //         } );

  //         /* check if they are all retimeable together */
  //         if ( free_spots >= total_fanout )
  //         {
  //           for ( auto f : fanout_splitters )
  //           {
  //             ntk.set_value( f, 1 );
  //           }
  //           return true;
  //         }
  //         /* select one randomly */
  //         std::default_random_engine gen( rseed++ );
  //         std::uniform_int_distribution<uint32_t> dist( 0ul, fanout_splitters.size() - 1 );
  //         auto index = dist( gen );
  //         ntk.set_value( fanout_splitters[index], 1 );
  //       }
  //     }
  //     return true;
  //   } );
  // }

  void create_generic_network( topo_view<Ntk>& topo, generic_network& res, node_map<signal_d, Ntk>& old2new )
  {
    topo.foreach_node( [&] ( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return true;

      std::vector<signal_d> children;
      
      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( _ntk.is_complemented( f ) )
          children.push_back( res.create_not( old2new[f] ) );
        else
          children.push_back( old2new[f] );
      } );

      if ( _ntk.is_buf( n ) && _ntk.visited( n ) == _ntk.trav_id() )
      {
        auto const in_latch = res.create_box_input( children[0] );
        auto const latch = res.create_latch( in_latch );
        auto const latch_out = res.create_box_output( latch );
        old2new[n] = latch_out;
      }
      else
      {
        const auto f = res.create_node( children, _ntk.node_function( n ) );
        old2new[n] = f;
      }

      return true;
    } );

    _ntk.foreach_po( [&]( auto const& f ) {
      if ( _ntk.is_complemented( f ) )
        res.create_po( res.create_not( old2new[f] ) );
      else
        res.create_po( old2new[f] );
    } );
  }

  uint32_t rec_fetch_root( node const n )
  {
    uint32_t value;

    _ntk.foreach_fanin( n, [&]( auto const f ) {
      auto g = _ntk.get_node( f );
      if ( !_ntk.is_buf( g ) || _ntk.fanout_size( g ) == 1 )
        value = _ntk.value( n );
      else
        value = rec_fetch_root( g );
    } );

    return value;
  }

  void rec_update_root( node const n, uint32_t const update )
  {
    if ( !_ntk.is_buf( n ) || _ntk.fanout_size( n ) == 1 )
      return;

    _ntk.set_value( n, update );
    _ntk.foreach_fanin( n, [&]( auto const f ) {
      rec_update_root( _ntk.get_node( f ), update );
    } );
  }

private:
  Ntk &_ntk;
  aqfp_network_conversion_params const& _ps;
};

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
generic_network generic_network_create_from_buffered( Ntk& ntk, aqfp_network_conversion_params const& ps = {} )
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

  detail::generic_network_create_from_buffered_impl p( ntk, ps );
  return p.run();

  // using signal = typename generic_network::signal;

  // node_map<signal, Ntk> old2new( ntk );
  // generic_network res;

  // old2new[ntk.get_constant( false )] = res.get_constant( false );
  // if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
  // {
  //   old2new[ntk.get_constant( true )] = res.get_constant( true );
  // }
  // ntk.foreach_pi( [&]( auto const& n ) {
  //   old2new[n] = res.create_pi();
  // } );

  // topo_view topo{ ntk };
  // fanout_view fntk{ ntk };

  // ntk.clear_values();
  // auto rseed = seed;

  /* set the nodes to be retimed to visited */
  /* set value `index` for splitters to test */
  // ntk.incr_trav_id();
  // topo.foreach_node( [&] ( auto const& n ) {
  //   if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
  //     return true;

  //   if ( ntk.is_buf( n ) )
  //   {
  //     if ( ntk.fanout_size( n ) == 1 )
  //     {
  //       ntk.set_visited( n, ntk.trav_id() );
  //     }
  //     else if ( !ntk.value( n ) && ntk.fanout_size( n ) < ps.splitter_capacity )
  //     {
  //       int free_spots = ps.splitter_capacity - ntk.fanout_size( n );
  //       int total_fanout = 0;
  //       std::vector<typename Ntk::node> fanout_splitters;

  //       /* select retimeable splitters */
  //       fntk.foreach_fanout( n, [&]( auto const f ) {
  //         if ( ntk.is_buf( f ) && ntk.fanout_size( f ) > 1 && free_spots >= ntk.fanout_size( f ) - 1 )
  //         {
  //           fanout_splitters.push_back( f );
  //           total_fanout += ntk.fanout_size( f ) - 1;
  //         }
  //       } );

  //       /* check if they are all retimeable together */
  //       if ( free_spots >= total_fanout )
  //       {
  //         for ( auto f : fanout_splitters )
  //         {
  //           ntk.set_value( f, 1 );
  //         }
  //         return true;
  //       }
  //       /* select one randomly */
  //       std::default_random_engine gen( rseed++ );
  //       std::uniform_int_distribution<uint32_t> dist( 0ul, fanout_splitters.size() - 1 );
  //       auto index = dist( gen );
  //       ntk.set_value( fanout_splitters[index], 1 );
  //     }
  //   }
  //   return true;
  // } );

  /* Alternative method: select buffers and splitters to retime */
  // ntk.incr_trav_id();
  // topo.foreach_node( [&] ( auto const& n ) {
  //   if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
  //     return true;

  //   if ( ntk.is_buf( n ) )
  //   {
  //     if ( ntk.fanout_size( n ) == 1 )
  //     {
  //       ntk.set_visited( n, ntk.trav_id() );
  //     }
  //     else if ( ntk.visited( n ) != ntk.trav_id() ||  ntk.value( n ) > 0 )
  //     {
  //       int free_spots;
  //       if ( ntk.value( n ) > 0 )
  //       {
  //         free_spots = detail::rec_fetch_root( ntk, n );
  //         if ( free_spots == 0 )
  //           return true;
  //       }
  //       else
  //       {
  //         free_spots = ps.splitter_capacity - ntk.fanout_size( n );
  //       }

  //       int total_fanout = 0;
  //       std::vector<typename Ntk::node> fanout_splitters;

  //       /* select retimeable splitters */
  //       fntk.foreach_fanout( n, [&]( auto const f ) {
  //         if ( ntk.is_buf( f ) && ntk.fanout_size( f ) > 1 && free_spots >= ntk.fanout_size( f ) - 1 )
  //         {
  //           fanout_splitters.push_back( f );
  //           total_fanout += ntk.fanout_size( f ) - 1;
  //         }
  //       } );

  //       /* check if they are all retimeable together */
  //       if ( free_spots >= total_fanout )
  //       {
  //         for ( auto f : fanout_splitters )
  //         {
  //           ntk.set_value( f, free_spots - total_fanout );
  //           ntk.set_visited( f, ntk.trav_id() );
  //         }
  //         detail::rec_update_root( ntk, n, free_spots - total_fanout );
  //         return true;
  //       }
  //       /* select one randomly */
  //       std::default_random_engine gen( rseed++ );
  //       std::uniform_int_distribution<uint32_t> dist( 0ul, fanout_splitters.size() - 1 );
  //       auto index = dist( gen );
  //       ntk.set_value( fanout_splitters[index], free_spots - ntk.fanout_size( fanout_splitters[index] ) + 1 );
  //       ntk.set_visited( fanout_splitters[index], ntk.trav_id() );
  //       detail::rec_update_root( ntk, n, free_spots - ntk.fanout_size( fanout_splitters[index] ) + 1 );
  //     }
  //   }
  //   return true;
  // } );

  // topo.foreach_node( [&] ( auto const& n ) {
  //   if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
  //     return true;

  //   std::vector<signal> children;
    
  //   ntk.foreach_fanin( n, [&]( auto const& f ) {
  //     if ( ntk.is_complemented( f ) )
  //       children.push_back( res.create_not( old2new[f] ) );
  //     else
  //       children.push_back( old2new[f] );
  //   } );

  //   if ( ntk.is_buf( n ) && ntk.visited( n ) == ntk.trav_id() )
  //   {
  //     auto const in_latch = res.create_box_input( children[0] );
  //     auto const latch = res.create_latch( in_latch );
  //     auto const latch_out = res.create_box_output( latch );
  //     old2new[n] = latch_out;
  //   }
  //   else
  //   {
  //     const auto f = res.create_node( children, ntk.node_function( n ) );
  //     old2new[n] = f;
  //   }

  //   return true;
  // } );

  // ntk.foreach_po( [&]( auto const& f ) {
  //   if ( ntk.is_complemented( f ) )
  //     res.create_po( res.create_not( old2new[f] ) );
  //   else
  //     res.create_po( old2new[f] );
  // } );

  // return res;
}

/*! \brief Network convertion from generic network.
 *
 * This function converts a generic_network to a
 * buffered mig network (buffered_mig_network).
 *
 * \param ntk Generic network
 */
buffered_mig_network buffered_mig_create_from_generic( generic_network const& ntk, aqfp_assumptions const& ps = {} )
{
  using signal = typename buffered_mig_network::signal;

  node_map<signal, generic_network> old2new( ntk );
  buffered_mig_network res;

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

    /* remove not represented nodes */
    if ( ntk.is_box_input( n ) ||  ntk.is_box_output( n ) || ntk.is_po( n ) )
    {
      signal children;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        children = old2new[f];
      } );

      // assert( ntk.fanout_size( n ) <= ps.splitter_capacity );

      // if ( ntk.is_box_output( n ) && fntk.fanout_size( n ) > 1 )
      // {
      //   std::cout << "Node " << n << " = ";
      //   kitty::print_hex( ntk.node_function( n ) );
      //   std::cout << " has fanout " <<  fntk.fanout( n ).size() << " and fanin " << ntk.get_fanin0( ntk.get_fanin0( ntk.get_fanin0 ( n ) ) ) << "\n";
      // }

      old2new[n] = children;
      return true;
    }

    std::vector<signal> children;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old2new[f] );
    } );

    signal f;

    if ( ntk.fanin_size( n ) == 3 )
    {
      /* majority */
      assert( children.size() == 3 );
      f = res.create_maj( children[0], children[1], children[2] );
    }
    else if ( ntk.fanin_size( n ) == 1 && ntk.node_function( n )._bits[0] == 0x1 )
    {
      /* not */
      assert( children.size() == 1 );
      f = !children[0];
    }
    else
    {
      /* buffer */
      assert( children.size() == 1 );
      assert( ntk.fanout_size( n ) <= ps.splitter_capacity );

      /* not balanced PIs */
      if ( !ps.balance_pis && ( res.is_pi( res.get_node( children[0] ) ) || res.is_constant( res.get_node( children[0] ) ) ) )
        f = children[0];
      else
        f = res.create_buf( children[0] );
    }

    old2new[n] = f;
    return true;
  } );

  ntk.foreach_po( [&]( auto const& f ) {
    res.create_po( old2new[f] );
  } );

  return res;
}

} // namespace mockturtle