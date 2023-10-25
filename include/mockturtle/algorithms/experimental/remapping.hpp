/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file remapping.hpp
  \brief Remapping engine

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include "detail/resyn_opt.hpp"
#include "detail/resyn_windowing.hpp"
#include "../../networks/block.hpp"
#include "../../networks/klut.hpp"
#include "../../traits.hpp"
#include "../../utils/algorithm.hpp"
#include "../../utils/cuts.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/tech_library.hpp"
#include "../../views/binding_view.hpp"
#include "../../views/cell_view.hpp"
#include "../../views/depth_view.hpp"
#include "../../views/fanout_view.hpp"
#include "../../views/topo_view.hpp"
#include "../../views/window_view.hpp"
#include "../cleanup.hpp"
#include "../cut_enumeration.hpp"
#include "../detail/mffc_utils.hpp"
#include "../detail/switching_activity.hpp"

namespace mockturtle
{

/*! \brief Parameters for remap.
 *
 * The data structure `remap_params` holds configurable parameters
 * with default arguments for `remap`.
 */
struct remap_params
{
  /*! \brief Do area-oriented remapping. */
  bool area_oriented_remapping{ false };

  /*! \brief Required time for delay optimization. */
  double required_time{ 0.0f };

  /*! \brief Maps using multi-output gates */
  bool use_multioutput{ false };

  /*! \brief Window number of PIs */
  uint32_t num_pis{ 12 };

  /*! \brief Window number of POs */
  uint32_t num_pos{ 12 };

  remap_params()
  {
    cut_enumeration_ps.cut_limit = 16;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 16.
   * The maximum cut limit is 15.
   * By default, truth table minimization
   * is performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for remap.
 *
 * The data structure `remap_stats` provides data collected by running
 * `remap`.
 */
struct remap_stats
{
  /*! \brief Recovered area. */
  double area_save{ 0 };
  /*! \brief Successful remappings. */
  uint32_t num_success{ 0 };
  /*! \brief Failed remappings. */
  uint32_t num_fail{ 0 };

  /*! \brief Area result. */
  double area{ 0 };
  /*! \brief Worst delay result. */
  double delay{ 0 };
  /*! \brief Power result. */
  double power{ 0 };

  /*! \brief Mapped multi-output gates. */
  uint32_t multioutput_gates{ 0 };

  /*! \brief Total runtime. */
  /*! \brief Total runtime. */
  stopwatch<>::duration time_windowing{ 0 };
  stopwatch<>::duration time_opt{ 0 };
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Rempping error. */
  bool remapping_error{ false };

  void report() const
  {
    std::cout << fmt::format( "[i] Save = {:>5.2f}; Area = {:>5.2f}; Delay = {:>5.2f};", area_save, area, delay );
    std::cout << fmt::format( "[i] Time W = {:>5.2f}s; Time O = {:>5.2f}s; Total = {:>5.2f}s\n", to_seconds( time_windowing ), to_seconds( time_opt ), to_seconds( time_total ) );
  }
};

namespace detail
{

template<class Ntk, typename WindowEngine, typename DecompFn, typename OptScript, typename MapperFn, unsigned NInputs, classification_type Configuration>
class remap_impl
{
private:
public:
  explicit remap_impl( Ntk& ntk, tech_library<NInputs, Configuration> const& library, remap_params const& ps, remap_stats& st )
      : ntk( ntk ),
        library( library ),
        ps( ps ),
        st( st )
  {}

  void run()
  {
    stopwatch t( st.time_total );

    resyn_windowing_params win_ps = ps.win_ps;
    WindowEngine win( ntk, win_ps );
    node_map<signal<Ntk>, Ntk> ntk_to_win( ntk );
    OptScript opt;

    ntk.foreach_gate( [&]( auto const& n ) {
      /* extract window */
      call_with_stopwatch( st.time_windowing, [&]() { win.compute_window( n ); } );
      if ( visited_window( win.get_hash() ) )
        return;
      // win.report_info();

      /* create new network instance */
      // Ntk win_ntk;
      // win_copy( win_ntk, win, ntk_to_win );

      /* decompose */


      /* optimize */
      uint32_t size_before = win_ntk.num_gates();
      call_with_stopwatch( st.time_opt, [&]() { opt( win_ntk ); } );

      /* map */
      // auto res = MapperFn( ).run();

      /* evaluate */
      if ( !evaluate( win_ntk, win ) )
      {
        st.num_fail++;
        return;
      }

      /* replace */
      replace( win_ntk, win, ntk_to_win );
      ++st.num_success;
      return;
    } );
  }

private:
  void win_copy( Ntk& win_ntk, WindowEngine const& win, node_map<signal<Ntk>, Ntk>& ntk_to_win )
  {
    std::vector<signal<Ntk>> children;
    children.reserve( Ntk::max_fanin_size );

    if ( ntk_to_win.size() != ntk.size() )
      ntk_to_win.resize();

    /* create PIs */
    for ( auto const& n : win.get_leaves() )
    {
      ntk_to_win[n] = win_ntk.create_pi();
    }

    /* create gates */
    for ( auto const& n : win.get_gates() )
    {
      children.clear();
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        signal<Ntk> s = ntk.is_complemented( f ) ? ntk.create_not( ntk_to_win[f] ) : ntk_to_win[f];
        children.push_back( s );
      } );
      ntk_to_win[n] = win_ntk.clone_node( ntk, n, children );
    }

    /* create POs */
    for ( auto const& f : win.get_roots() )
    {
      if ( ntk.is_complemented( f ) )
      {
        win_ntk.create_po( ntk.create_not( ntk_to_win[f] ) );
      }
      else
      {
        win_ntk.create_po( ntk_to_win[f] );
      }
    }
  }

  bool evaluate( Ntk const& win_ntk, WindowEngine const& win )
  {
    if ( win_ntk.num_gates() < win.num_gates() )
    {
      st.size_save += win.num_gates() - win_ntk.num_gates();
      return true;
    }

    return false;
  }

  void replace( Ntk& win_ntk, WindowEngine const& win, node_map<signal<Ntk>, Ntk>& ntk_to_win )
  {
    topo_view topo{ win_ntk };

    std::vector<signal<Ntk>> children;
    children.reserve( Ntk::max_fanin_size );

    /* Get PIs */
    uint32_t i = 0;
    for ( auto const& n : win.get_leaves() )
    {
      ntk_to_win[win_ntk.pi_at( i++ )] = ntk.make_signal( n );
    }

    /* add the new gates */
    topo.foreach_gate( [&]( auto const& n ) {
      children.clear();
      win_ntk.foreach_fanin( n, [&]( auto const& f ) {
        signal<Ntk> s = win_ntk.is_complemented( f ) ? ntk.create_not( ntk_to_win[f] ) : ntk_to_win[f];
        children.push_back( s );
      } );
      ntk_to_win[n] = ntk.clone_node( win_ntk, n, children );
    } );

    /* substitute POs */
    auto const& roots = win.get_roots();
    win_ntk.foreach_po( [&]( auto const& f, uint32_t index ) {
      signal<Ntk> s = ntk.is_complemented( f ) ? ntk.create_not( ntk_to_win[f] ) : ntk_to_win[f];
      ntk.substitute_node( ntk.get_node( roots.at( index ) ), s );
    } );
  }

  bool visited_window( std::array<uint64_t, 2> const& hash )
  {
    return !window_cache.insert( hash ).second;
  }

private:
  Ntk& ntk;
  tech_library<NInputs, Configuration> const& library; /* TODO: replace with &&Lib? */
  remap_params const& ps;
  remap_stats& st;

  phmap::flat_hash_set<std::array<uint64_t, 2>, window_hash_fn> window_cache;
};

} /* namespace detail */

/*! \brief Remapping.
 *
 * This function implements a simple technology mapping algorithm.
 * The algorithm maps each node to the first implementation in the technology library.
 *
 * The input must be a binding_view with the gates correctly loaded.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `has_binding`
 *
 * \param ntk Network
 *
 */
/* TODO: add resynFn */
template<class Ntk, typename OptScript = detail::resyn_aig_size<Ntk>, unsigned NInputs, classification_type Configuration>
void remap( Ntk& ntk, tech_library<NInputs, Configuration> const& library, remap_params const& ps = {}, remap_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  // static_assert( has_has_binding_v<Ntk> || has_has_cell_v<Ntk>, "Ntk does not implement the has_binding or has_cell method" );

  using remap_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> d_ntk{ ntk };
  remap_view_t remap_ntk{ d_ntk };

  using WindowEngine = detail::remap_windowing<remap_view_t>;
  using MapperFn = detail::emap_impl<remap_view_t, 6, NInputs, Configuration>;

  remap_stats st;
  detail::remap_impl<remap_view_t, WindowEngine, OptScript, MapperFn, NInputs, Configuration> p( remap_ntk, library, ps, st );
  p.run();

  if ( ps.verbose && !st.remapping_error )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */