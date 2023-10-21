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
  \file rremapping.hpp
  \brief A rremapping engine

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
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Rempping error. */
  bool remapping_error{ false };

  void report() const
  {
    std::cout << fmt::format( "[i] Area = {:>5.2f}; Delay = {:>5.2f};", area, delay );
    if ( power != 0 )
      std::cout << fmt::format( " Power = {:>5.2f};\n", power );
    else
      std::cout << "\n";
    if ( multioutput_gates )
    {
      std::cout << fmt::format( "[i] Multi-output gates   = {:>5}\n", multioutput_gates );
    }
    std::cout << fmt::format( "[i] Total runtime        = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

struct remap_windowing_params
{
  /*! \brief Window number of PIs */
  uint32_t num_pis{ 12 };

  /*! \brief Window number of POs */
  uint32_t num_pos{ 12 };

  /*! \brief TFI max levels */
  uint32_t tfi_levels{ 4 };

  /*! \brief TFO max levels */
  uint32_t tfo_levels{ 3 };

  /*! \brief Maximum number of gates to include in a window. */
  uint32_t max_gates{ 10 };

  /*! \brief Maximum fanout of a node to expand. */
  uint32_t skip_fanout_limit{ 5 };
};

template<class Ntk>
class remap_windowing
{
private:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit remap_windowing( Ntk const& ntk, remap_windowing_params const& ps )
      : ntk( ntk ), ps( ps )
  {
    leaves.reserve( ps.max_gates );
    roots.reserve( ps.max_gates );
    gates.reserve( ps.max_gates );
    candidates.reserve( ps.max_gates );
  }

  void compute_window( node const& pivot )
  {
    leaves.clear();
    roots.clear();
    gates.clear();

    /* add pivot to gates */
    ntk.incr_trav_id();
    ntk.set_visited( pivot, ntk.trav_id() );
    gates.push_back( pivot );

    if ( ps.max_gates < 2 )
      return;

    /* decrement fanout size of leaves */
    ntk.foreach_fanin( pivot, [&]( auto const& f ) {
      ntk.decr_fanout_size( ntk.get_node( f ) );
    } );

    /* increment traverse ID */
    ntk.incr_trav_id();

    /* add iteratively nodes to the window */
    std::optional<node> next;
    while ( gates.size() < ps.max_gates && ( next = find_next_pivot() ) )
    {
      assert( ntk.visited( *next ) < ntk.trav_id() - 1 );
      gates.push_back( *next );
      ntk.set_visited( *next, ntk.trav_id() - 1 );

      /* decrement fanout size of leaves */
      ntk.foreach_fanin( *next, [&]( auto const& f ) {
        ntk.decr_fanout_size( ntk.get_node( f ) );
      } );
    }

    /* restore fanout counts */
    for ( node const& n : gates )
    {
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        ntk.incr_fanout_size( ntk.get_node( f ) );
      } );
    }

    /* collect roots, leaves, and gates in topo order */
    collect_roots();
    collect_nodes();

    assert( gates.size() <= ps.max_gates );
    return;
  }

  std::vector<node> const& get_gates() const
  {
    return gates;
  }

  void report_info( std::ostream& out = std::cout )
  {
    out << fmt::format( "[i] W: I = {};\t G = {};\t O = {}\n", leaves.size(), gates.size(), roots.size() );
  }

private:
  std::optional<node> find_next_pivot()
  {
    candidates.clear();
    auto best = candidates.cend();

    /* marks meaning */
    /* visited == trav_id - 1 --> in the window */
    /* visited == trav_id --> candidate in the TFI of the window */

    /* look for MFFC nodes */
    do
    {
      for ( node const& n : gates )
      {
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          node const& g = ntk.get_node( f );
          if ( ntk.visited( g ) < ntk.trav_id() - 1 && ntk.fanout_size( g ) == 0 && !ntk.is_ci( g ) )
          {
            candidates.push_back( g );
            ntk.set_visited( g, ntk.trav_id() );
          }
        } );
      }

      if ( !candidates.empty() )
      {
        /* select the best candidate */
        best = max_element_unary(
            candidates.begin(), candidates.end(),
            [&]( auto const& cand ) {
              auto cnt{ 0 };
              ntk.foreach_fanin( cand, [&]( auto const& f ) {
                cnt += ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() ? 1 : 0;
              } );
              return cnt;
            },
            -1 );
        break;
      }

      /* add all the input candidates */
      for ( node const& n : gates )
      {
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          node const& g = ntk.get_node( f );
          if ( ntk.visited( g ) < ntk.trav_id() - 1 && !ntk.is_ci( g ) )
          {
            candidates.push_back( g );
            ntk.set_visited( g, ntk.trav_id() );
          }
        } );
      }

      /* add all the output candidates */
      for ( node const& n : gates )
      {
        if ( ntk.fanout_size( n ) == 0 || ntk.fanout_size( n ) > ps.skip_fanout_limit )
          continue;

        /* Output candidate member of the MFFC */
        auto const& fanout_v = ntk.fanout( n );
        if ( ntk.fanout_size( n ) == 1 && fanout_v.size() == 1 && ntk.visited( fanout_v.front() ) < ntk.trav_id() - 1 )
        {
          candidates.push_back( fanout_v.front() );
          best = candidates.cend() - 1;
          break;
        }

        /* does not mark the nodes (TFI has the priority to avoid reconvergences passing from outside the window) */
        std::copy_if( fanout_v.begin(), fanout_v.end(),
                      std::back_inserter( candidates ),
                      [&]( auto const& g ) {
                        return ntk.visited( g ) < ntk.trav_id() - 1;
                      } );
      }

      if ( !candidates.empty() )
      {
        /* select the best candidate */
        best = max_element_unary(
            candidates.begin(), candidates.end(),
            [&]( auto const& cand ) {
              auto cnt{ 0 };
              ntk.foreach_fanin( cand, [&]( auto const& f ) {
                cnt += ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() ? 1 : 0;
              } );
              return cnt;
            },
            -1 );
      }
    } while ( false );

    if ( candidates.empty() )
    {
      return std::nullopt;
    }

    assert( best != candidates.end() );

    /* reset the marks */
    for ( auto const& n : candidates )
    {
      ntk.set_visited( n, ntk.trav_id() - 2 );
    }

    return *best;
  }

  void collect_roots()
  {
    for ( node const& n : gates )
    {
      /* equivalent to is_po() */
      if ( ntk.fanout_size( n ) != ntk.fanout( n ).size() )
      {
        roots.push_back( n );
        continue;
      }

      ntk.foreach_fanout( n, [&]( auto const& g ) {
        if ( ntk.visited( g ) < ntk.trav_id() - 1 )
        {
          roots.push_back( n );
          return false;
        }
        return true;
      } );
    }
  }

  void collect_nodes()
  {
    auto prev_size = gates.size();
    gates.clear();

    /* collects gates and leaves */
    for ( node const& n : roots )
    {
      if ( ntk.visited( n ) == ntk.trav_id() )
        continue;

      collect_nodes_rec( n );
    }

    assert( gates.size() == prev_size );
  }

  void collect_nodes_rec( node const& n )
  {
    assert( ntk.visited( n ) != ntk.trav_id() );
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      node const& g = ntk.get_node( f );

      if ( ntk.visited( g ) < ntk.trav_id() - 1 )
      {
        /* leaf */
        leaves.push_back( g );
        ntk.set_visited( g, ntk.trav_id() );
      }
      else if ( ntk.visited( g ) == ntk.trav_id() - 1 )
      {
        /* gate */
        collect_nodes_rec( g );
      }
    } );

    gates.push_back( n );
  }

private:
  Ntk const& ntk;
  remap_windowing_params const& ps;

  uint32_t num_pis{ 0 };

  std::vector<node> leaves;
  std::vector<node> roots;
  std::vector<node> gates;
  std::vector<node> candidates;
};

template<class Ntk, typename WindowEngine, typename MapperFn, unsigned NInputs, classification_type Configuration>
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

    remap_windowing_params win_ps;
    remap_windowing<Ntk> win( ntk, win_ps );

    ntk.foreach_gate( [&]( auto const& n ) {
      /* extract window */
      win.compute_window( n );
      win.report_info();
      ++st.num_success;
      // window_view<Ntk> win{ ntk };

      /* optimize */

      /* remap */

      /* evaluate */
      // st.num_fail++;
    } );
  }

private:
  Ntk& ntk;
  // WindowsEngine engine;
  // ResynFn resyn;
  tech_library<NInputs, Configuration> const& library; /* TODO: replace with &&Lib? */
  remap_params const& ps;
  remap_stats& st;
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
template<class Ntk, unsigned CutSize = 6u, unsigned NInputs, classification_type Configuration>
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
  static_assert( has_has_binding_v<Ntk> || has_has_cell_v<Ntk>, "Ntk does not implement the has_binding or has_cell method" );

  using remap_view_t = fanout_view<depth_view<Ntk>>;
  depth_view<Ntk> d_ntk{ ntk };
  remap_view_t remap_ntk{ d_ntk };

  using WindowEngine = detail::remap_windowing<remap_view_t>;
  using MapperFn = detail::emap_impl<remap_view_t, CutSize, NInputs, Configuration>;

  remap_stats st;
  detail::remap_impl<remap_view_t, WindowEngine, MapperFn, NInputs, Configuration> p( remap_ntk, library, ps, st );
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