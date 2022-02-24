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
#include <list>
#include <algorithm>

#include "aqfp_assumptions.hpp"
#include "buffer_insertion.hpp"
#include "aqfp_network_convertion.hpp"
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
  using splitter_tuple = std::tuple<signal, node, int>;

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

    Ntk ntk = cleanup_dangling_buffered( _ntk );

    /* reposition buffers */
    push_buffers_forward( ntk );

    bool success = false;
    ntk.clear_values();
    if ( achievable_depth < current_depth )
    {
      fanout_view<Ntk> f_ntk{ ntk };
      success = run_cut_based_depth_reduction( f_ntk, current_depth - achievable_depth );
    }

    if ( success )
    {
      node_map<signal, Ntk> old2new( ntk );
      Ntk res;
      create_res_net( ntk, res, old2new );
      ntk = res;
    }

    if ( _ps.critical_depth_reduction )
    {
      run_critical_depth_reduction( ntk );
    }

    /* splitter trees reconstruction params */
    buffer_insertion_params buf_ps;
    buf_ps.assume = _ps.aqfp_assumptions_ps;
    buf_ps.scheduling = buffer_insertion_params::provided;
    buf_ps.optimization_effort = buffer_insertion_params::none;
    return aqfp_reconstruct_splitter_trees( ntk, buf_ps );
  }

private:
  template<class FNtk>
  bool run_cut_based_depth_reduction( FNtk& ntk, uint32_t iterations )
  {
    static_assert( has_foreach_fanout_v<FNtk>, "Ntk does not implement the foreach_fanout method" );

    /* find a cut of buffers and mark them as removable */
    uint32_t i;
    for ( i = 1; i <= iterations; ++i )
    {
      ntk.incr_trav_id();
      uint32_t trav_id = ntk.trav_id();

      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), trav_id );

      /* mark nodes to define a cut */
      ntk.foreach_pi( [&]( auto const& n ) {
        mark_cut_rec( ntk, n );
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

      if ( _ps.verbose )
      {
        std::cout << fmt::format( "[i] Initial depth = {:7d}\t Final depth = {:7d}\r", _st.depth_pre, _st.depth_post ) << std::flush;
      }
    }

    /* no cut found */
    if ( i == 1 )
    {
      return false;
    }

    return true;
  }

  void run_critical_depth_reduction( Ntk& ntk )
  {
    aqfp_level_t d_ntk{ ntk };
    fanout_view<aqfp_level_t> f_ntk{ d_ntk };

    bool success = false;
    while ( true )
    {
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
        {
          mark_cut_critical_rec( f_ntk, n );
        }
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

      if ( legal_cut )
      {
        /* TODO: a more "in depth" study */
        ntk.foreach_po( [&]( auto const& f ) {
          if ( ntk.value( ntk.get_node( f ) ) && ntk.fanout_size( ntk.get_node( f ) ) > 1 )
          {
            /* check validity */
            legal_cut = false;
          }
          return legal_cut;
        } );
      }

      if ( !legal_cut )
      {
        /* critical path cannot be reduced */
        break;
      }

      /* modify selected splitter trees and critical section */
      std::vector<node> critical_cut;
      change_splitter_trees2( f_ntk, critical_cut );

      lower_critical_section( f_ntk, critical_cut );

      /* remove cut of buffers */
      auto result = run_cut_based_depth_reduction( f_ntk, 1 );

      if ( !result )
        break;
      else
        success = true;

      remove_buffers_inplace( f_ntk );
      f_ntk.update_levels();
    }
  }

  template<class FNtk>
  void mark_cut_rec( fanout_view<FNtk>& f_ntk, node const& n )
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

  template<class FNtk>
  bool select_buf_cut_rec( FNtk& ntk, node const& n, uint32_t value )
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
      if ( f_ntk.fanout_size( n ) == 1 || check_cut_critical_splitter2( f_ntk, n ) )
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

  /* old version: fast, but finds less optimization opportunities that the newer version */
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

  /* new version: slightly slower, but finds more optimization opportunities that the old version */
  inline bool check_cut_critical_splitter2( fanout_view<aqfp_level_t>& ntk, node const& n )
  {
    /* TODO: extend using required time */
    node fanin;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      fanin = ntk.get_node( f );
    } );

    /* return if not a splitter tree root */
    if ( ntk.is_buf( fanin ) )
      return false;

    std::vector<int> level_assignment;

    bool modify = collect_splitter_tree_leaves_levels( ntk, n, 0, level_assignment );

    /* no need to rewrite the splitter tree */
    if ( modify == false )
      return true;

    /* sort vector by level in decreasing order */
    std::sort( level_assignment.begin(), level_assignment.end(), std::greater<int>() );

    /* check if negative level (not valid) */
    if ( level_assignment.empty() || level_assignment.back() < 0 )
      return false;

    /* see if the new level assignment has a solution */
    uint32_t nodes_in_level = 0;
    uint32_t last_level = level_assignment.front();
    for ( int const l : level_assignment )
    {
      if ( l == last_level )
      {
        ++nodes_in_level;
      }
      else
      {
        /* update splitters */
        for ( auto i = 0; ( i < last_level - l ) && ( nodes_in_level != 1 ); ++i )
          nodes_in_level = std::ceil( float( nodes_in_level ) / float( _ps.aqfp_assumptions_ps.splitter_capacity ) );

        ++nodes_in_level;
        last_level = l;
      }
    }
    for ( auto i = 0; i < last_level; ++i )
      nodes_in_level = std::ceil( float( nodes_in_level ) / float( _ps.aqfp_assumptions_ps.splitter_capacity ) );

    if ( nodes_in_level > _ps.aqfp_assumptions_ps.splitter_capacity )
      return false;
    else
      return true;
  }

  bool collect_splitter_tree_leaves_levels( fanout_view<aqfp_level_t> const& ntk, node const& n, int level, std::vector<int>& level_assignment )
  {
    bool modify = false;
    ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( ntk.is_buf( f ) && ntk.fanout_size( f ) > 1 )
      {
        modify |= collect_splitter_tree_leaves_levels( ntk, f, level + 1, level_assignment );
      }
      else
      {
        /* lower critical signal by one ( if not a buffer ) */
        if ( ntk.is_on_critical_path( f ) && !ntk.is_buf( f ) )
        {
          level_assignment.push_back( level - 1 );
          modify = true;
        }
        else
        {
          level_assignment.push_back( level );
        }
      }
    } );

    /* TODO: consider POs, for now POs are considered as balanced */
    return modify;
  }

  bool collect_splitter_tree_leaves( fanout_view<aqfp_level_t> const& ntk, node const& n, int level, std::vector<splitter_tuple>& signal_assignment, bool phase )
  {
    bool modify = false;
    ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( ntk.is_buf( f ) && ntk.fanout_size( f ) > 1 )
      {
        bool phase_s = phase;
        ntk.foreach_fanin( n, [&]( auto const& fanin ) {
          phase_s ^= ntk.is_complemented( fanin );
        } );
        modify |= collect_splitter_tree_leaves( ntk, f, level + 1, signal_assignment, phase_s );
      }
      else
      {
        /* lower critical signal by one ( if not a buffer ) */
        if ( ntk.is_on_critical_path( f ) && !ntk.is_buf( f ) )
        {
          signal_assignment.emplace_back( std::make_tuple( ntk.make_signal( f ) ^ phase, n, level - 1 ) );
          modify = true;
        }
        else
        {
          signal_assignment.emplace_back( std::make_tuple( ntk.make_signal( f ) ^ phase, n, level ) );
        }
      }
    } );

    /* TODO: consider POs, for now POs are considered as balanced */
    return modify;
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

  /* old version: fast, but finds less optimization opportunities that the newer version */
  void change_splitter_trees( fanout_view<aqfp_level_t>& ntk, std::vector<node>& critical_cut )
  {
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.value( n ) == 1 )
      {
        if ( ntk.fanout_size( n ) > 1 )
        {
          signal fanin;
          ntk.foreach_fanin( n, [&]( auto const& f ) {
            fanin = f;
          } );
          for ( auto const& f : ntk.fanout( n ) )
          {
            if ( ntk.is_on_critical_path( f ) )
            {
              auto const buf = ntk.create_buf( fanin );
              ntk.replace_in_node( f, n, buf );
              ntk.decr_fanout_size( n );
              /* add to critical path */
              ntk.set_level( ntk.get_node( buf ), ntk.level( n ) );
              ntk.set_on_critical_path( ntk.get_node( buf ), true );
              /* add to critical cut */
              critical_cut.push_back( ntk.get_node( buf ) );
            }
          }

          /* remove n from critiacl path */
          ntk.set_on_critical_path( n, false );

          if ( ntk.fanout_size( n ) == 0 )
          {
            ntk.take_out_node( n );
          }
        }
        else
        {
          critical_cut.push_back( n );
        }
      }
    } );
  }

  /* new version: slightly slower, but finds more optimization opportunities that the old version */
  void change_splitter_trees2( fanout_view<aqfp_level_t>& ntk, std::vector<node>& critical_cut )
  {
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.value( n ) == 1 )
      {
        if ( ntk.fanout_size( n ) > 1 )
        {
          /* reconstruct splitter tree lowering the critical paths */
          std::vector<splitter_tuple> signal_assignment;

          bool modify = collect_splitter_tree_leaves( ntk, n, 0, signal_assignment, false );

          if ( modify == false )
          {
            /* no need to rewrite the splitter tree, just collect the critical cut */
            for ( auto const& t : signal_assignment )
            {
              if ( ntk.is_on_critical_path( ntk.get_node( std::get<0>( t ) ) ) )
              {
                critical_cut.push_back( ntk.get_node( std::get<0>( t ) ) );
              }
            }
            return;
          }

          std::sort( signal_assignment.begin(), signal_assignment.end(), []( auto const& a, auto const& b ) {
            return std::get<2>( a ) > std::get<2>( b );
          } );

          uint32_t max_level = std::get<2>( signal_assignment.front() );
          std::vector<uint32_t> splitters_per_level( max_level, 0 );
          uint32_t nodes_in_level = 0;
          uint32_t last_level = max_level;

          for ( auto const& t : signal_assignment )
          {
            auto l = std::get<2>( t );
            if ( l == max_level )
            {
              ++nodes_in_level;
            }
            else
            {
              /* update splitters */
              for ( auto i = last_level; i > l; --i )
              {
                splitters_per_level[i - 1] = std::ceil( float( nodes_in_level ) / float( _ps.aqfp_assumptions_ps.splitter_capacity ) );
                nodes_in_level = splitters_per_level[i - 1];
              }

              ++nodes_in_level;
              last_level = l;
            }
          }
          for ( auto i = last_level; i > 0; --i )
          {
            splitters_per_level[i - 1] = std::ceil( float( nodes_in_level ) / float( _ps.aqfp_assumptions_ps.splitter_capacity ) );
            nodes_in_level = splitters_per_level[i - 1];
          }

          /* get root node */
          signal root_s;
          ntk.foreach_fanin( n, [&]( auto const& f ) {
            root_s = f;
          } );

          std::vector<std::list<signal>> splitters( max_level + 1 );
          splitters[0].push_back( ntk.create_buf( root_s ) );

          /* create splitter tree */
          for ( auto i = 0; i < splitters_per_level.size(); ++i )
          {
            auto it = splitters[i].begin();
            for ( auto j = 0; j < splitters_per_level[i]; ++j )
            {
              splitters[i + 1].push_back( ntk.create_buf( *it ) );
              if ( ntk.fanout_size( ntk.get_node( *it ) ) == _ps.aqfp_assumptions_ps.splitter_capacity )
                ++it;
            }
          }

          /* assign signals from splitter trees */
          last_level = max_level;
          auto it = splitters[max_level].begin();
          for ( auto const& t : signal_assignment )
          {
            if ( std::get<2>( t ) != last_level )
            {
              it = splitters[std::get<2>( t )].begin();
              while ( ntk.fanout_size( ntk.get_node( *it ) ) >= _ps.aqfp_assumptions_ps.splitter_capacity )
                ++it;
              last_level = std::get<2>( t );
            }
            signal f = std::get<0>( t );
            if ( ntk.is_on_critical_path( ntk.get_node( f ) ) && !ntk.is_buf( ntk.get_node( f ) ) )
            {
              signal buf = ntk.create_buf( *it ) ^ ntk.is_complemented( f );
              ntk.replace_in_node( ntk.get_node( f ), std::get<1>( t ), buf );
              critical_cut.push_back( ntk.get_node( buf ) );
              set_critical_path_fanin_rec( ntk, ntk.get_node( buf ) );
            }
            else
            {
              ntk.replace_in_node( ntk.get_node( f ), std::get<1>( t ), *it ^ ntk.is_complemented( f ) );
              if ( ntk.is_on_critical_path( ntk.get_node( f ) ) )
              {
                critical_cut.push_back( ntk.get_node( f ) );
                set_critical_path_fanin_rec( ntk, ntk.get_node( *it ) );
              }
            }
            if ( ntk.fanout_size( ntk.get_node( *it ) ) >= _ps.aqfp_assumptions_ps.splitter_capacity )
              ++it;
          }

          /* take out nodes */
          for ( auto const& t : signal_assignment )
          {
            if ( !ntk.is_dead( std::get<1>( t ) ) )
              ntk.take_out_node( std::get<1>( t ) );
          }
        }
        else
        {
          critical_cut.push_back( n );
        }
      }
    } );
  }

  void set_critical_path_fanin_rec( fanout_view<aqfp_level_t>& ntk, node const& n )
  {
    if ( ntk.is_on_critical_path( n ) )
      return;

    ntk.set_on_critical_path( n, true );
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      set_critical_path_fanin_rec( ntk, ntk.get_node( f ) );
    } );
  }

  void lower_critical_section( fanout_view<aqfp_level_t>& ntk, std::vector<node> const& critical_cut )
  {
    /* remove TFI of critical cut from being critical */
    ntk.incr_trav_id();
    for ( auto n : critical_cut )
    {
      node g;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        g = ntk.get_node( f );
      } );
      reset_on_critical_path_tfi( ntk, g );
    }

    /* find blocking path buffers */
    ntk.incr_trav_id();
    ntk.clear_values();
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
    ntk.foreach_pi( [&]( auto const& n ) {
      visit_and_mark_tfo_buffer_rec( ntk, n );
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
      return incompatibilities;
    } );

    /* the cut is legal or configuration is not valid */
    if ( !incompatibilities )
    {
      return;
    }

    move_critical_section_down( ntk );
  }

  void reset_on_critical_path_tfi( fanout_view<aqfp_level_t>& ntk, node const& n )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.set_visited( n, ntk.trav_id() );
    ntk.set_on_critical_path( n, false );

    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto g = ntk.get_node( f );
      if ( !ntk.is_constant( g ) && ntk.is_on_critical_path( g ) )
        reset_on_critical_path_tfi( ntk, g );
    } );
  }

  void visit_and_mark_tfo_buffer_rec( fanout_view<aqfp_level_t>& f_ntk, node_g const& n )
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
    if ( f_ntk.is_buf( n ) && f_ntk.fanout_size( n ) == 1 )
    {
      return;
    }

    /* recur towards TFO if not in critical section after the cut */
    if ( !f_ntk.is_on_critical_path( n ) )
    {
      f_ntk.foreach_fanout( n, [&]( auto const& f ) {
        if ( f_ntk.visited( f ) != f_ntk.trav_id() )
          visit_and_mark_tfo_buffer_rec( f_ntk, f );
      } );
    }
  }

  void move_critical_section_down( fanout_view<aqfp_level_t>& ntk )
  {
    ntk.clear_values();

    /* patch critical section fanout */
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.visited( n ) == ntk.trav_id() )
      {
        if ( ntk.is_buf( n ) && ntk.fanout_size( n ) == 1 )
        {
          node fanin;
          ntk.foreach_fanin( n, [&]( auto const& f ) {
            fanin = ntk.get_node( f );
          } );
          if ( ntk.visited( fanin ) != ntk.trav_id() )
          {
            ntk.set_value( n, 1 );
          }
        }
        auto splitter = ntk.get_constant( false );
        for ( auto const& f : ntk.fanout( n ) )
        {
          if ( ntk.visited( f ) != ntk.trav_id() )
          {
            if ( !ntk.is_on_critical_path( n ) )
            {
              if ( splitter == ntk.get_constant( false ) )
              {
                splitter = ntk.create_buf( ntk.make_signal( n ) );
              }
              ntk.replace_in_node( f, n, splitter );
              ntk.decr_fanout_size( n );
            }
            else
            {
              auto buf = ntk.create_buf( ntk.make_signal( n ) );
              ntk.replace_in_node( f, n, buf );
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
        signal fanin;
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = f;
        } );
        for ( auto const& f : ntk.fanout( n ) )
        {
          ntk.replace_in_node( f, n, fanin );
          ntk.take_out_node( n );
        }
      }
    } );
  }

  void mark_critical_section_tfo( fanout_view<aqfp_level_t>& f_ntk, node_g const& n )
  {
    if ( f_ntk.visited( n ) == f_ntk.trav_id() )
      return;

    f_ntk.set_visited( n, f_ntk.trav_id() );

    f_ntk.foreach_fanout( n, [&]( auto const& f ) {
      if ( f_ntk.is_on_critical_path( f ) && f_ntk.visited( f ) == f_ntk.trav_id() - 1 )
        mark_critical_section_tfo( f_ntk, f );
    } );

    if ( f_ntk.is_buf( n ) && f_ntk.fanout_size( n ) == 1 )
    {
      return;
    }
    
    f_ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( !f_ntk.is_constant( f_ntk.get_node( f ) ) )
        mark_critical_section_tfo( f_ntk, f_ntk.get_node( f ) );
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

    if ( legal_cut )
    {
      ntk.foreach_po( [&]( auto const& f ) {
        if ( ntk.value( ntk.get_node( f ) ) && ntk.fanout_size( ntk.get_node( f ) ) > 1 )
          legal_cut = false;
        return legal_cut;
      } );
    }

    if ( !legal_cut )
      return;

    std::vector<node> critical_cut;

    change_splitter_trees2( f_ntk, critical_cut );
    lower_critical_section( f_ntk, critical_cut );

    /* remove cut of buffers */
    auto result = run_cut_based_depth_reduction( ntk, 1 );

    /* create the new network */
    Ntk res_local;
    node_map<signal, Ntk> old2new( ntk );

    create_res_net( ntk, res_local, old2new );

    ntk = res_local;
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
      if ( f_ntk.fanout_size( n ) == 1 || check_cut_critical_splitter2( f_ntk, n ) )
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

      assert( children.size() > 0 );

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
        f = res.clone_node( ntk, n, children );
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

  template<class FNtk>
  void remove_buffers_inplace( FNtk& ntk )
  {
    static_assert( has_foreach_fanout_v<FNtk>, "Ntk does not implement the foreach_fanout method" );

    ntk.foreach_node( [&]( auto const& n ) {
      /* remove selected buffers */
      if ( ntk.value( n ) )
      {
        signal fanin;
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          fanin = f;
        } );
        
        assert( ntk.fanout_size( n ) == 1 );

        std::vector<node> fanout = ntk.fanout( n );
        
        if ( fanout.empty() )
        {
          /* PO */
          ntk.replace_in_outputs( n, fanin );
        }
        else
        {
          ntk.replace_in_node( fanout[0], n, fanin );
        }
        ntk.take_out_node( n );
      }
    } );
  }

  void push_buffers_forward( Ntk& ntk )
  {
    /* ntk must be topologically sorted */

    /* collect the buffers (latches) */
    std::vector<node> buffers;
    buffers.reserve( 100 );

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_buf( n ) && ntk.fanout_size( n ) == 1 )
        buffers.push_back( n );
    } );

    fanout_view<Ntk> f_ntk{ ntk };

    /* reverse topological order */
    for ( auto it = buffers.rbegin(); it != buffers.rend(); ++it )
    {
      for ( auto const g : f_ntk.fanout( *it ) )
      {
        /* output splitter */
        if ( ntk.fanout_size( g ) != 1 )
        {
          forward_push_rec( f_ntk, g );
          /* remove current buffer */
          signal fanin;
          ntk.foreach_fanin( *it, [&]( auto const& f ) {
            fanin = f;
          } );
          f_ntk.substitute_node( *it, fanin );
        }
      }
    }
  }

  bool forward_push_rec( fanout_view<Ntk>& ntk, node const n )
  {
    auto const fanouts = ntk.fanout( n );
    for ( auto const& f : fanouts )
    {
      if ( ntk.fanout_size( f ) == 1 )
      {
        auto buf = ntk.create_buf( ntk.make_signal( n ) );
        ntk.replace_in_node( f, n, buf );
        ntk.decr_fanout_size( n );
      }
      else
      {
        forward_push_rec( ntk, f );
      }
    }
    /* PO */
    if ( fanouts.size() == 0 )
    {
      /* set it as a fanin */
      signal fanin;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        fanin = f ^ ntk.is_complemented( f );
      } );
      auto buf = ntk.create_buf( fanin );
      ntk.replace_in_node( n, ntk.get_node( fanin ), buf );
      ntk.decr_fanout_size( ntk.get_node( fanin ) );
    }
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
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_is_complemented_v<Ntk>, "NtkDest does not implement the is_complemented method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the substitute_node method" );
  static_assert( has_replace_in_node_v<Ntk>, "Ntk does not implement the replace_in_node method" );
  static_assert( has_take_out_node_v<Ntk>, "Ntk does not implement the take_out_node method" );
  static_assert( is_buffered_network_type_v<Ntk>, "Ntk is not a buffered network type" );
  static_assert( has_is_buf_v<Ntk>, "Ntk does not implement the is_buf method" );

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

class aqfp_reconstruct_splitter_trees_impl
{
public:
  using node = typename aqfp_network::node;
  using signal = typename aqfp_network::signal;

public:
  aqfp_reconstruct_splitter_trees_impl( buffered_aqfp_network const& ntk, buffer_insertion_params const& ps )
    : _ntk( ntk )
    , _ps( ps )
    {}

  buffered_aqfp_network run()
  {
    /* save the level of each node */
    depth_view ntk_level{ _ntk };

    /* create a network removing the splitter trees */
    aqfp_network clean_ntk;
    node_map<signal, buffered_aqfp_network> old2new( _ntk );
    remove_splitter_trees( clean_ntk, old2new );

    /* compute the node level on the new network */
    node_map<uint32_t, aqfp_network> levels( clean_ntk );
    _ntk.foreach_gate( [&]( auto const& n ) {
      levels[old2new[n]] = ntk_level.level( n );
    } );

    /* recompute splitter trees and return the new buffered network */
    buffered_aqfp_network res;
    buffer_insertion buf_inst( clean_ntk, levels, _ps );
    buf_inst.run( res );
    return res;
  }

private:
  void remove_splitter_trees( aqfp_network& res, node_map<signal, buffered_aqfp_network>& old2new )
  {
    topo_view topo{ _ntk };

    old2new[_ntk.get_constant( false )] = res.get_constant( false );
 
    _ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    topo.foreach_node( [&]( auto const& n ) {
      if ( _ntk.is_pi( n ) || _ntk.is_constant( n ) )
        return;

      std::vector<signal> children;
      _ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( old2new[f] ^ _ntk.is_complemented( f ) );
      } );

      if ( _ntk.is_buf( n ) )
      {
        old2new[n] = children[0];
      }
      else if ( children.size() == 3 )
      {
        old2new[n] = res.create_maj( children[0], children[1], children[2] );
      }
      else
      {
        old2new[n] = res.create_maj( children );
      }
    } );

    _ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( old2new[f] ^ _ntk.is_complemented( f ) );
    } );
  }

private:
  buffered_aqfp_network const& _ntk;
  buffer_insertion_params const& _ps;
};

} /* namespace detail */

/*! \brief Rebuilds buffer trees in AQFP network.
 *
 * This function rebuilds buffer trees in AQFP network.
 *
 * \param ntk Buffered AQFP network
 */
buffered_aqfp_network aqfp_reconstruct_splitter_trees( buffered_aqfp_network& ntk, buffer_insertion_params const& ps = {} )
{
  detail::aqfp_reconstruct_splitter_trees_impl p( ntk, ps );
  return p.run();
}

} // namespace mockturtle