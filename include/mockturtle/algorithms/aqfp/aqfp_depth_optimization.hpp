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

struct aqfp_optimize_depth_params
{
  /*! \brief AQFP technology assumptions. */
  aqfp_assumptions aqfp_assumptions_ps{};

  /*! \brief Critical depth reduction. */
  bool critical_depth_reduction{ true };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for AQFP depth optimization.
 *
 * The data structure `aqfp_optimize_depth_stats` provides data collected by running
 * `aqfp_optimize_depth`.
 */
struct aqfp_optimize_depth_stats
{
  /*! \brief Initial number of buffers/splitters. */
  // uint32_t buffers_pre{ 0 };

  /*! \brief Number of buffers/splitters after the algorithm. */
  // uint32_t buffers_post{ 0 };

  /*! \brief Initial depth. */
  uint32_t depth_pre{ 0 };

  /*! \brief Final depth. */
  uint32_t depth_post{ 0 };

  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  void report() const
  {
    // std::cout << fmt::format( "[i] Initial B/S   = {:7d}\t Final B/S   = {:7d}\n", buffers_pre, buffers_post );
    std::cout << fmt::format( "[i] Initial depth = {:7d}\t Final depth = {:7d}\n", depth_pre, depth_post );
    std::cout << fmt::format( "[i] Total runtime = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

template<class Ntk>
struct aqfp_depth_cost
{
  uint32_t operator()( Ntk const& ntk, node<Ntk> const& node ) const
  {
    if ( ntk.is_buf( node ) && ntk.fanout_size( node ) == 1 )
      return 0u;
    else
      return 1u;
  }
};

template<class Ntk>
class aqfp_optimize_depth_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using node_g = typename generic_network::node;
  using signal_g = typename generic_network::signal;
  using aqfp_level_t = depth_view<Ntk, aqfp_depth_cost<Ntk>>;

public:
  explicit aqfp_optimize_depth_impl( Ntk& ntk, aqfp_optimize_depth_params const& ps, aqfp_optimize_depth_stats& st )
    : _ntk( ntk )
    , _ps( ps )
    , _st( st )
    {}

public:
  Ntk run()
  {
    stopwatch t( _st.time_total );

    /* get real depth */
    auto achievable_depth = aqfp_level_t( _ntk ).depth();
    auto current_depth = depth_view<Ntk>( _ntk ).depth();

    _st.depth_pre = current_depth;
    _st.depth_post = current_depth;

    bool success = true;
    _ntk.clear_values();
    if ( achievable_depth < current_depth )
    {
      success = run_cut_based_depth_reduction( _ntk, current_depth - achievable_depth );
    }

    /* create resulting network */
    node_map<signal, Ntk> old2new( _ntk );
    Ntk res;

    create_res_net( _ntk, res, old2new );

    if ( !success )
    {
      /* push buffers over splitters to separate paths */
      auto current_depth2 = depth_view<Ntk>( res ).depth();
      Ntk ntk_forward = aqfp_push_buffers_forward( res );

      run_cut_based_depth_reduction( ntk_forward, current_depth2 - achievable_depth );

      Ntk res2;
      node_map<signal, Ntk> old2new2( ntk_forward );
      create_res_net( ntk_forward, res2, old2new2 );

      res = aqfp_push_buffers_backward( res2 );
    }

    if ( _ps.critical_depth_reduction )
    {
      res = run_critical_depth_reduction( res );
    }

    /* experiment */
    try_splitter_trees_repositioning( res );

    return res;
  }

private:
  bool run_cut_based_depth_reduction( Ntk& ntk, uint32_t iterations )
  {
    fanout_view<Ntk> f_ntk{ ntk };

    /* find a cut of buffers and mark them as removable */
    uint32_t i;
    for ( i = 1; i <= iterations; ++i )
    {
      ntk.incr_trav_id();
      uint32_t trav_id = ntk.trav_id();

      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), trav_id );

      /* mark nodes to define a cut */
      ntk.foreach_pi( [&]( auto const& n ) {
        mark_cut_rec( f_ntk, n );
      } );

      /* extract a cut if it exist */
      ntk.incr_trav_id();
      bool legal_cut = true;
      ntk.foreach_po( [&]( auto const& f ) {
        legal_cut = select_buf_cut_rec( ntk, ntk.get_node( f ), i );
        return legal_cut;
      } );

      if ( !legal_cut )
      {
        /* depth reduction is not a cut, undo last iteration and exit */
        ntk.foreach_node( [&]( auto const& n ) {
          if ( ntk.value( n ) == i )
            ntk.set_value( n, 0 );
        } );
        break;
      }

      --_st.depth_post;
    }

    /* last iteration is not a cut */
    if ( i != iterations + 1 )
    {
      return false;
    }

    return true;
  }

  void mark_cut_rec( fanout_view<Ntk>& f_ntk, node const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f_ntk.get_node( f ) ) != f_ntk.trav_id() )
      {
        mark_cut_rec( f_ntk, f_ntk.get_node( f ) );
      }
    } );

    /* found a new possible buffer cut */
    if ( f_ntk.is_buf( n ) && f_ntk.fanout_size( n ) == 1 && !f_ntk.value( n ) )
    {
      return;
    }

    /* recur towards TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f ) != f_ntk.trav_id() )
      {
        mark_cut_rec( f_ntk, f );
      }
    } );
  }

  bool select_buf_cut_rec( Ntk& ntk, node const& n, uint32_t value )
  {
    if ( ntk.is_constant( n ) )
      return true;

    if ( ntk.visited( n ) == ntk.trav_id() )
      return true;

    /* if selected buffer, set as removable */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 && ntk.is_buf( n ) && ntk.fanout_size( n ) == 1 )
    {
      ntk.set_visited( n, ntk.trav_id() );
      /* already selected in the past iterations */
      if ( ntk.value( n ) != 0 && ntk.value( n ) != value )
        return false;

      ntk.set_value( n, value );
      return true;
    }

    /* check not a cut */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 )
    {
      ntk.set_visited( n, ntk.trav_id() );
      return false;
    }

    ntk.set_visited( n, ntk.trav_id() );

    bool legal = true;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      legal = select_buf_cut_rec( ntk, ntk.get_node( f ), value );
      return legal;
    } );

    return legal;
  }

  Ntk run_critical_depth_reduction( Ntk& ntk_start )
  {
    bool success = false;

    Ntk ntk = aqfp_push_buffers_forward( ntk_start );
    Ntk res = ntk_start;

    bool additional_try = true;
    while ( true )
    {
      std::cout << "-------------------\n";
      aqfp_level_t d_ntk{ ntk };
      fanout_view<aqfp_level_t> f_ntk{ d_ntk };

      ntk.clear_values();
      ntk.incr_trav_id();
      uint32_t trav_id = ntk.trav_id();

      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), trav_id );

      /* set free spots foreach splitter */
      ntk.foreach_node( [&]( auto const& n ) {
        if ( ntk.is_buf( n ) )
          ntk.set_value( n, _ps.aqfp_assumptions_ps.splitter_capacity - ntk.fanout_size( n ) );
      } );

      /* A cut of buffers/splitters on the critical paths may exist */
      ntk.foreach_pi( [&]( auto const& n ) {
        if ( d_ntk.is_on_critical_path( n ) )
          mark_cut_critical_rec( f_ntk, n );
      } );

      /* search for the critical cut */
      bool legal_cut = true;
      ntk.clear_values();
      ntk.incr_trav_id();
      ntk.foreach_po( [&]( auto const& f ) {
        if ( d_ntk.is_on_critical_path( ntk.get_node( f ) ) )
          legal_cut = select_buf_cut_critical_rec( d_ntk, ntk.get_node( f ), 1 );
        return legal_cut;
      } );

      if ( !legal_cut )
      {
        /* critical path cannot be reduced */
        break;
      }

      /* mark critical splitters to reposition */
      ntk.foreach_node( [&]( auto const& n ) {
        if ( ntk.value( n ) == 1 )
        {
          if ( ntk.fanout_size( n ) > 1 )
          {
            f_ntk.foreach_fanout( n, [&]( auto const& f ) {
              if ( d_ntk.is_on_critical_path( f ) )
                d_ntk.set_value( f, 2 );
            } );
          }
          else
          {
            ntk.set_value( n, 0 );
          }
        }
      } );

      /* modify selected splitter trees and critical section */
      auto generic_net = convert_to_generic( d_ntk );
      fanout_view<generic_network> f_generic_net{ generic_net };
      change_splitter_trees( f_generic_net );
      lower_critical_section( f_generic_net );
      auto ntk_rewired = convert_to_buffered( generic_net );

      std::cout << aqfp_level_t( ntk_rewired ).depth() << std::endl;

      /* remove cut of buffers */
      auto result = run_cut_based_depth_reduction( ntk_rewired, 1 );

      if ( !result )
      {
        if ( !additional_try )
          break;
        else
          additional_try = false;
      }
      else
        success = true;

      /* create the new network */
      Ntk res_local;
      node_map<signal, Ntk> old2new( ntk_rewired );

      create_res_net( ntk_rewired, res_local, old2new );

      res = res_local;
      ntk = res_local;
    }

    if ( success )
      return aqfp_push_buffers_backward( res );
    else
      return res;
  }

  void mark_cut_critical_rec( fanout_view<aqfp_level_t>& f_ntk, node const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards critical TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto g = f_ntk.get_node( f );
      if ( f_ntk.visited( g ) != f_ntk.trav_id() && f_ntk.is_on_critical_path( g ) )
      {
        mark_cut_critical_rec( f_ntk, g );
      }
    } );

    /* find a cut */
    if ( f_ntk.is_buf( n ) )
    {
      if ( f_ntk.fanout_size( n ) == 1 || check_cut_critical_splitter( f_ntk, n ) )
      {
        return;
      }
    }

    /* recur towards critical TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f ) != f_ntk.trav_id() && f_ntk.is_on_critical_path( f ) )
      {
        mark_cut_critical_rec( f_ntk, f );
      }
    } );
  }

  inline bool check_cut_critical_splitter( fanout_view<aqfp_level_t>& f_ntk, node const& n )
  {
    /* check for input splitter */
    bool valid = false;
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      node g = f_ntk.get_node( f );
      if ( f_ntk.is_buf( g ) && f_ntk.value( g ) > 0 )
      {
        /* count current splitter critical signals */
        uint32_t count = 0;
        f_ntk.foreach_fanout( n, [&]( auto const& fanout ) {
          if ( f_ntk.is_on_critical_path( fanout ) )
            ++count;
        } );

        /* decrease if removable splitter */
        if ( count == f_ntk.fanout_size( n ) )
          --count;

        if ( f_ntk.value( g ) >= count )
        {
          f_ntk.set_value( g, f_ntk.value( g ) - count );
          valid = true;
        }
      }
    } );

    return valid;
  }

  bool select_buf_cut_critical_rec( aqfp_level_t& ntk, node const& n, uint32_t value )
  {
    if ( ntk.is_constant( n ) )
      return true;

    if ( ntk.visited( n ) == ntk.trav_id() )
      return true;

    /* if selected buffer, set as removable */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 && ntk.is_buf( n ) )
    {
      ntk.set_visited( n, ntk.trav_id() );
      ntk.set_value( n, value );
      return true;
    }

    /* check not a cut */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 )
    {
      ntk.set_visited( n, ntk.trav_id() );
      return false;
    }

    ntk.set_visited( n, ntk.trav_id() );

    bool legal = true;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_on_critical_path( ntk.get_node( f ) ) )
        legal = select_buf_cut_critical_rec( ntk, ntk.get_node( f ), value );
      return legal;
    } );

    return legal;
  }

  void change_splitter_trees( fanout_view<generic_network>& ntk )
  {
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.value( n ) == 1 )
      {
        auto fanin = ntk.get_fanin0( n );
        for ( auto const& f : ntk.fanout( n ) )
        {
          if ( ntk.value( f ) == 2 )
          {
            auto const buf = ntk.create_latch( fanin );
            ntk.replace_in_node( f, n, buf );
            ntk.decr_fanout_size( n );
            /* add to critical path */
            ntk.set_value2( ntk.get_node( buf ), 1 );
          }
        }

        if ( ntk.fanout_size( n ) == 0 )
        {
          ntk.take_out_node( n );
        }
        /* remove from critical path */
        ntk.set_value2( n, 0 );
      }
    } );
  }

  void lower_critical_section( fanout_view<generic_network>& ntk )
  {
    std::vector<node_g> critical_cut;

    /* get critical cut */
    ntk.incr_trav_id();
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
    ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
    ntk.foreach_pi( [&]( auto const& n ) {
      if ( ntk.value2( n ) )
        mark_cut_critical_generic_rec( ntk, n );
    } );

    ntk.incr_trav_id();
    ntk.foreach_po( [&]( auto const& f ) {
      auto n = ntk.get_node( ntk.get_fanin0( ntk.get_node( f ) ) );
      if ( ntk.value2( n ) )
        select_critical_cut_rec( ntk, n, critical_cut );
    } );

    /* remove TFI of critical cut from being critical */
    for ( auto n : critical_cut )
    {
      reset_value2_tfi( ntk, ntk.get_node( ntk.get_fanin0( n ) ) );
    }

    /* find blocking path buffers */
    ntk.incr_trav_id();
    ntk.clear_values();
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
    ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
    ntk.foreach_pi( [&]( auto const& n ) {
      visit_and_mark_tfo_buffer_rec( ntk, n );
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      auto n = ntk.get_node( ntk.get_fanin0( ntk.get_node( f ) ) );
      if ( ntk.visited( n ) == ntk.trav_id() && !ntk.is_latch( n ) )
        std::cout << n << " visited" << std::endl;
    } );

    /* find lower boundary (cut) of the critical section */
    ntk.incr_trav_id();
    bool incompatibilities = false;
    for ( auto n : critical_cut )
    {
      ntk.foreach_fanout( n, [&]( auto const& f ) {
        if ( ntk.visited( f ) == ntk.trav_id() - 1 )
        {
          incompatibilities = true;
          mark_critical_section_tfo( ntk, n );
        }
      } );
    }

    /* check validity */
    ntk.foreach_pi( [&]( auto const& n ) {
      if ( ntk.visited( n ) == ntk.trav_id() )
        incompatibilities = false;
    } );

    /* the cut is legal */
    if ( !incompatibilities )
    {
      return;
    }

    move_critical_section_down( ntk );
  }

  void mark_cut_critical_generic_rec( fanout_view<generic_network>& f_ntk, node const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards critical TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto const g = f_ntk.get_node( f );
      if ( f_ntk.visited( g ) != f_ntk.trav_id() && f_ntk.value2( g ) )
      {
        mark_cut_critical_generic_rec( f_ntk, g );
      }
    } );

    /* find a cut */
    if ( f_ntk.is_latch( n ) )
      return;

    /* recur towards critical TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f ) != f_ntk.trav_id() && f_ntk.value2( f ) )
      {
        mark_cut_critical_generic_rec( f_ntk, f );
      }
    } );
  }

  bool select_critical_cut_rec( fanout_view<generic_network>& ntk, node const& n, std::vector<node_g>& cut )
  {
    if ( ntk.is_constant( n ) )
      return true;

    if ( ntk.visited( n ) == ntk.trav_id() )
      return true;

    /* if selected buffer, set as removable */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 && ntk.is_latch( n ) )
    {
      ntk.set_visited( n, ntk.trav_id() );
      cut.push_back( n );
      return true;
    }

    /* check not a cut */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 )
    {
      std::cout << "Failed: " << n << "\n";
      ntk.set_visited( n, ntk.trav_id() );
      return false;
    }

    ntk.set_visited( n, ntk.trav_id() );

    bool legal = true;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.value2( ntk.get_node( f ) ) )
        legal = select_critical_cut_rec( ntk, ntk.get_node( f ), cut );
      return legal;
    } );

    return legal;
  }

  void reset_value2_tfi( fanout_view<generic_network>& ntk, node_g const& n )
  {
    ntk.set_value2( n, 0 );

    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto g = ntk.get_node( f );
      if ( !ntk.is_constant( g ) && ntk.value2( g ) )
        reset_value2_tfi( ntk, g );
    } );
  }

  void visit_and_mark_tfo_buffer_rec( fanout_view<generic_network>& f_ntk, node_g const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards critical TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto const g = f_ntk.get_node( f );
      if ( f_ntk.visited( g ) != f_ntk.trav_id() )
      {
        visit_and_mark_tfo_buffer_rec( f_ntk, g );
      }
    } );

    /* stop tfo recursion */
    if ( f_ntk.is_latch( n ) )
    {
      return;
    }

    /* recur towards TFO if not in critical section after the cut */
    if ( !f_ntk.value2( n ) )
    {
      f_ntk.foreach_fanout( n, [&]( auto const& f ) {
        if ( f_ntk.visited( f ) != f_ntk.trav_id() )
          visit_and_mark_tfo_buffer_rec( f_ntk, f );
      } );
    }
  }

  void move_critical_section_down( fanout_view<generic_network>& ntk )
  {
    ntk.clear_values();

    /* patch critical section fanout */
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.visited( n ) == ntk.trav_id() )
      {
        if ( ntk.is_latch( n ) && ntk.visited( ntk.get_node( ntk.get_fanin0( n ) ) ) != ntk.trav_id() )
        {
          ntk.set_value( n, 1 );
        }
        auto buf = ntk.get_constant( false );
        for ( auto const& f : ntk.fanout( n ) )
        {
          if ( ntk.visited( f ) != ntk.trav_id() )
          {
            if ( !ntk.value2( n ) )
            {
              if ( buf == ntk.get_constant( false ) )
              {
                buf = ntk.create_buf( n );
              }
              ntk.replace_in_node( f, n, buf );
              ntk.decr_fanout_size( n );
            }
            else
            {
              auto latch = ntk.create_latch( n );
              ntk.replace_in_node( f, n, latch );
              ntk.decr_fanout_size( n );
            }
          }
        }
      }
    } );

    /* remove lower buffer cut in critical section */
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.value( n ) )
      {
        auto fanin = ntk.get_fanin0( n );
        for ( auto const& f : ntk.fanout( n ) )
        {
          ntk.replace_in_node( f, n, fanin );
          ntk.take_out_node( n );
        }
      }
    } );
  }

  void mark_critical_section_tfo( fanout_view<generic_network>& f_ntk, node_g const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.value2( f ) && f_ntk.visited( f ) == f_ntk.trav_id() - 1 )
        mark_critical_section_tfo( f_ntk, f );
    } );

    if ( f_ntk.is_latch( n ) )
    {
      std::cout << n << std::endl;
      return;
    }
    
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( !f_ntk.is_constant( f_ntk.get_node( f ) ) )
          mark_critical_section_tfo( f_ntk, f );
    } );
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

  void try_splitter_trees_repositioning( Ntk& ntk )
  {
    aqfp_level_t d_ntk{ ntk };
    node_map<uint32_t, Ntk> req_time( ntk );

    std::vector<node> topo_order;
    topo_order.reserve( ntk.size() );

    ntk.foreach_node( [&]( auto const& n ) {
      topo_order.push_back( n );
      req_time[n] = UINT32_MAX;
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      req_time[f] = d_ntk.depth();
    } );

    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if ( ntk.is_pi( *it ) || ntk.is_constant( *it ) )
        continue;

      uint32_t update = req_time[*it];
      if ( !ntk.is_buf( *it ) || ntk.fanout_size( *it ) > 1 )
        --update;

      ntk.foreach_fanin( *it, [&]( auto const& f ) {
        req_time[f] = std::min( req_time[f], update );
      } );
    }

    fanout_view<aqfp_level_t> f_ntk{ d_ntk };

    /* set free spots foreach splitter */
    ntk.clear_values();
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_buf( n ) )
        ntk.set_value( n, _ps.aqfp_assumptions_ps.splitter_capacity - ntk.fanout_size( n ) );
    } );

    /* A cut of buffers/splitters on the critical paths may exist */
    ntk.incr_trav_id();
    ntk.foreach_pi( [&]( auto const& n ) {
      if ( d_ntk.is_on_critical_path( n ) )
        mark_cut_critical_rec_experiment( f_ntk, n, req_time );
    } );

    /* search for the critical cut */
    bool legal_cut = true;
    ntk.clear_values();
    ntk.incr_trav_id();
    ntk.foreach_po( [&]( auto const& f ) {
      if ( d_ntk.is_on_critical_path( ntk.get_node( f ) ) )
        legal_cut = select_buf_cut_critical_rec( d_ntk, ntk.get_node( f ), 1 );
      return legal_cut;
    } );

    std::cout << fmt::format( "Experiment legal cut : {}\n", legal_cut );
  }

  void mark_cut_critical_rec_experiment( fanout_view<aqfp_level_t>& f_ntk, node const& n, node_map<uint32_t, Ntk> const& req_time )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    /* recur towards critical TFI */
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto g = f_ntk.get_node( f );
      if ( f_ntk.visited( g ) != f_ntk.trav_id() && f_ntk.is_on_critical_path( g ) )
      {
        mark_cut_critical_rec_experiment( f_ntk, g, req_time );
      }
    } );

    /* find a cut */
    if ( f_ntk.is_buf( n ) )
    {
      if ( f_ntk.fanout_size( n ) == 1 || check_cut_critical_splitter_experiment( f_ntk, n, req_time ) )
      {
        return;
      }
    }

    /* recur towards critical TFO */
    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.visited( f ) != f_ntk.trav_id() && f_ntk.is_on_critical_path( f ) )
      {
        mark_cut_critical_rec_experiment( f_ntk, f, req_time );
      }
    } );
  }

  inline bool check_cut_critical_splitter_experiment( fanout_view<aqfp_level_t>& f_ntk, node const& n, node_map<uint32_t, Ntk> const& req_time )
  {
    /* check for input splitter */
    bool valid = false;
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      node g = f_ntk.get_node( f );
      if ( f_ntk.is_buf( g ) && f_ntk.value( g ) > 0 )
      {
        /* count current splitter critical signals */
        uint32_t count = 0;
        f_ntk.foreach_fanout( n, [&]( auto const& fanout ) {
          if ( f_ntk.is_on_critical_path( fanout ) )
            ++count;
        } );

        /* decrease if removable splitter */
        if ( count == f_ntk.fanout_size( n ) )
          --count;

        if ( f_ntk.value( g ) >= count )
        {
          f_ntk.set_value( g, f_ntk.value( g ) - count );
          valid = true;
        }
      }
    } );

    return valid;
  }

  // bool check_balancing()
  // {
  //   depth_view<Ntk, aqfp_depth_cost_balancing> d_ntk{ _ntk };
  //   bool balanced = true;

  //   for ( auto const& n : _topo_order )
  //   {
  //     if ( _ntk.value( n ) )
  //       continue;

  //     _ntk.foreach_fanin( n, [&]( auto const& f ) {
  //       if ( !_ntk.is_constant( _ntk.get_node( f ) ) && d_ntk.level( _ntk.get_node( f ) ) != d_ntk.level( n ) - 1 )
  //       {
  //         balanced = false;
  //       }
  //     } );

  //     if ( !balanced )
  //       return false;
  //   }

  //   return true;
  // }

  void create_res_net( Ntk& ntk, Ntk& res, node_map<signal, Ntk>& old2new )
  {
    old2new[ntk.get_constant( false )] = res.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = res.get_constant( true );
    }
    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    topo_view topo{ ntk };
    topo.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        return;

      std::vector<signal> children;

      ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( old2new[f] ^ ntk.is_complemented( f ) );
      } );

      signal f;
      if ( ntk.is_buf( n ) )
      {
        if ( !ntk.value( n ) )
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
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
        res.create_po( res.create_not( old2new[f] ) );
      else
        res.create_po( old2new[f] );
    } );
  }

  generic_network convert_to_generic( aqfp_level_t const& ntk )
  {
    node_map<signal_g, Ntk> old2new( ntk );
    generic_network res;

    old2new[ntk.get_constant( false )] = res.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = res.get_constant( true );
    }
    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();

      if ( ntk.is_on_critical_path( n ) )
        res.set_value2( res.get_node( old2new[n] ), 1 );
    } );

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        return;

      std::vector<signal_g> children;

      ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( ntk.is_complemented( f ) )
        {
          auto not_s = res.create_not( old2new[f] );

          if ( ntk.value( n ) == 2 )
            res.set_value( res.get_node( not_s ), 2 );

          if ( ntk.is_on_critical_path( n ) )
            res.set_value2( res.get_node( not_s ), 1 );

          children.push_back( not_s );
        }
        else
        {
          children.push_back( old2new[f] );
        }
      } );

      if ( ntk.is_buf( n ) && ntk.fanout_size( n ) == 1 )
      {
        /* create registers (without box-in and box-out) */
        auto const latch = res.create_latch( children[0] );
        old2new[n] = latch;
      }
      else
      {
        const auto f = res.create_node( children, ntk.node_function( n ) );
        old2new[n] = f;
      }

      /* mark splitter trees */
      if ( ntk.value( n ) )
      {
        res.set_value( res.get_node( old2new[n] ), ntk.value( n ) );
      }
      if ( ntk.is_on_critical_path( n ) )
      {
        res.set_value2( res.get_node( old2new[n] ), 1 );
      }
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
      {
        auto not_s = res.create_not( old2new[f] );

        if ( ntk.value( ntk.get_node( f ) ) == 2 )
          res.set_value( res.get_node( not_s ), 2 );
        
        if ( ntk.is_on_critical_path( ntk.get_node( f ) ) )
          res.set_value2( res.get_node( not_s ), 1 );

        res.create_po( not_s );
      }
      else
      {
        res.create_po( old2new[f] );
      }
    } );

    return res;
  }

  Ntk convert_to_buffered( generic_network& ntk )
  {
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
        return;

      /* remove not represented nodes */
      if ( ntk.is_po( n ) )
      {
        old2new[n] = old2new[ntk.get_fanin0( n )];
        return;
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
        f = res.create_buf( children[0] );
      }

      old2new[n] = f;
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( old2new[f] );
    } );

    return res;
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
  aqfp_optimize_depth_params const& _ps;
  aqfp_optimize_depth_stats& _st;
};

} /* namespace detail */

/*! \brief Depth optimization for AQFP networks.
 *
 * This function tries to reduce the depth of a mapped AQFP circuit
 *
 * \param ntk Mapped AQFP network
 */
template<class Ntk>
Ntk aqfp_optimize_depth( Ntk& ntk, aqfp_optimize_depth_params const& ps = {}, aqfp_optimize_depth_stats *pst = nullptr )
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

  aqfp_optimize_depth_stats st;
  detail::aqfp_optimize_depth_impl p( ntk, ps, st );
  Ntk res = p.run();

  if ( ps.verbose )
    st.report();

  if ( pst )
    *pst = st;

  return res;
}


namespace detail
{

template<class Ntk>
class aqfp_push_buffers_forward_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using node_g = typename generic_network::node;
  using signal_g = typename generic_network::signal;

public:
  aqfp_push_buffers_forward_impl( Ntk const& ntk )
    : _ntk( ntk ) {}

  Ntk run()
  {
    /* convert to generic (to deactivate strashing on node substitute) */
    auto res = to_generic();
    push_forward( res );
    /* convert back to buffered network */
    return to_buffered( res );
  }

private:
  void push_forward( generic_network& ntk )
  {
    /* collect the buffers (latches) */
    std::vector<node_g> buffers;
    buffers.reserve( ntk.num_latches() );

    ntk.foreach_latch( [&]( auto const& n ) {
      buffers.push_back( n );
    } );

    fanout_view f_ntk{ntk};

    /* move buffers (reverse order) */
    for ( auto it = buffers.rbegin(); it != buffers.rend(); ++it )
    {
      node_g g = f_ntk.fanout( *it )[0];
      if ( ntk.fanout_size( g ) != 1 || ntk.node_function( g )._bits[0] == 0x1 )
      {
        forward_push_rec( f_ntk, g );
        /* remove current register */
        signal_g fanin = ntk.get_fanin0( *it );
        f_ntk.substitute_node( *it, fanin );
      }
    }
  }

  void forward_push_rec( fanout_view<generic_network>& ntk, node_g const& g )
  {
    auto const fanouts = ntk.fanout( g );
    for ( auto const& f : fanouts )
    {
      if ( ntk.fanout_size( f ) == 1 && ntk.node_function( f )._bits[0] != 0x1 )
      {
        auto buf = ntk.create_latch( ntk.make_signal( g ) );
        ntk.replace_in_node( f, g, buf );
        ntk.decr_fanout_size( g );
      }
      else
      {
        forward_push_rec( ntk, f );
      }
    }
  }

  generic_network to_generic()
  {
    node_map<signal_g, Ntk> old2new( _ntk );
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

    topo.foreach_node( [&]( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return;

      std::vector<signal_g> children;

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( _ntk.is_complemented( f ) )
        {
          children.push_back( res.create_not( old2new[f] ) );
        }
        else
        {
          children.push_back( old2new[f] );
        }
      } );

      if ( _ntk.is_buf( n ) && _ntk.fanout_size( n ) == 1 )
      {
        /* create registers (without box-in and box-out) */
        auto const latch = res.create_latch( children[0] );
        old2new[n] = latch;
      }
      else
      {
        const auto f = res.create_node( children, topo.node_function( n ) );
        old2new[n] = f;
      }
    } );

    _ntk.foreach_po( [&]( auto const& f ) {
      if ( _ntk.is_complemented( f ) )
        res.create_po( res.create_not( old2new[f] ) );
      else
        res.create_po( old2new[f] );
    } );

    return res;
  }

  Ntk to_buffered( generic_network& ntk )
  {
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
        return;

      /* remove not represented nodes */
      if ( ntk.is_po( n ) )
      {
        old2new[n] = old2new[ntk.get_fanin0( n )];
        return;
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
        f = res.create_buf( children[0] );
      }

      old2new[n] = f;
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( old2new[f] );
    } );

    return res;
  }

private:
  Ntk const& _ntk;
};

} /* namespace detail */


/*! \brief Moves buffers forward in AQFP networks.
 *
 * This function pushes buffers forward over splitter trees
 *
 * \param ntk Mapped AQFP network
 */
template<class Ntk>
Ntk aqfp_push_buffers_forward( Ntk& ntk )
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

  detail::aqfp_push_buffers_forward_impl p( ntk );
  return p.run();
}

namespace detail
{

template<class Ntk>
class aqfp_push_buffers_backward_impl
{
public:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;
  using node_g = typename generic_network::node;
  using signal_g = typename generic_network::signal;

public:
  aqfp_push_buffers_backward_impl( Ntk const& ntk )
    : _ntk( ntk ) {}

  Ntk run()
  {
    /* convert to generic (to deactivate strashing on node substitute) */
    auto res = to_generic();
    push_backward( res );
    /* convert back to buffered network */
    return to_buffered( res );
  }

private:
  void push_backward( generic_network& ntk )
  {
    fanout_view f_ntk{ ntk };

    /* move buffers in topological order */
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        return;

      if ( ntk.fanout_size( n ) > 1 )
      {
        while ( check_push_backward( f_ntk, n ) )
        {
          backward_push_rec( f_ntk, n );
          /* add buffer */
          signal_g fanin = ntk.get_fanin0( n );
          auto buf = f_ntk.create_latch( fanin );
          f_ntk.replace_in_node( n, fanin, buf );
          f_ntk.decr_fanout_size( fanin );
        }
      }
    } );
  }

  bool check_push_backward( fanout_view<generic_network>& ntk, node_g const& n )
  {
    bool valid = true;

    ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( !ntk.is_latch( f ) )
      {
        if ( ntk.fanout_size( f ) == 1 && ntk.node_function( f )._bits[0] != 0x1 )
          valid = false;
        else
          valid &= check_push_backward( ntk, f );
      }

      return valid;
    } );

    return valid;
  }

  void backward_push_rec( fanout_view<generic_network>& ntk, node_g const& n )
  {
    auto const fanouts = ntk.fanout( n );
    for ( auto const& f : fanouts )
    {
      if ( !ntk.is_latch( f ) )
      {
        backward_push_rec( ntk, f );
      }
      else
      {
        /* remove */
        ntk.substitute_node( f, ntk.make_signal( n ) );
      }
    }
  }

  generic_network to_generic()
  {
    node_map<signal_g, Ntk> old2new( _ntk );
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

    topo.foreach_node( [&]( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return;

      std::vector<signal_g> children;

      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        if ( _ntk.is_complemented( f ) )
        {
          children.push_back( res.create_not( old2new[f] ) );
        }
        else
        {
          children.push_back( old2new[f] );
        }
      } );

      if ( _ntk.is_buf( n ) && _ntk.fanout_size( n ) == 1 )
      {
        /* create registers (without box-in and box-out) */
        auto const latch = res.create_latch( children[0] );
        old2new[n] = latch;
      }
      else
      {
        const auto f = res.create_node( children, topo.node_function( n ) );
        old2new[n] = f;
      }
    } );

    _ntk.foreach_po( [&]( auto const& f ) {
      if ( _ntk.is_complemented( f ) )
        res.create_po( res.create_not( old2new[f] ) );
      else
        res.create_po( old2new[f] );
    } );

    return res;
  }

  Ntk to_buffered( generic_network& ntk )
  {
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
        return;

      /* remove not represented nodes */
      if ( ntk.is_po( n ) )
      {
        old2new[n] = old2new[ntk.get_fanin0( n )];
        return;
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
        f = res.create_buf( children[0] );
      }

      old2new[n] = f;
    } );

    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( old2new[f] );
    } );

    return res;
  }

private:
  Ntk const& _ntk;
};

} /* namespace detail */

/*! \brief Moves buffers backward in AQFP networks.
 *
 * This function pushes buffers backward before splitter trees
 *
 * \param ntk Mapped AQFP network
 */
template<class Ntk>
Ntk aqfp_push_buffers_backward( Ntk& ntk )
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

  detail::aqfp_push_buffers_backward_impl p( ntk );
  return p.run();
}

} // namespace mockturtle