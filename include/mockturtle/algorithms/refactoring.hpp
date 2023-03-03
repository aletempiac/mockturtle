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
  \file refactoring.hpp
  \brief Refactoring

  \author Alessandro Tempia Calvino
  \author Eleonora Testa
  \author Heinz Riener
  \author Mathias Soeken
  \author Siang-Yun (Sonia) Lee
*/
#pragma once

#include "../networks/mig.hpp"
#include "../traits.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/node_map.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/cut_view.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "../views/mffc_view.hpp"
#include "../views/topo_view.hpp"
#include "../views/window_view.hpp"
#include "../views/color_view.hpp"
#include "cleanup.hpp"
#include "detail/mffc_utils.hpp"
#include "dont_cares.hpp"
#include "simulation.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>

namespace mockturtle
{

/*! \brief Parameters for refactoring.
 *
 * The data structure `refactoring_params` holds configurable parameters with
 * default arguments for `refactoring`.
 */
struct refactoring_params
{
  /*! \brief Maximum number of PIs of the MFFC or window. */
  uint32_t max_pis{ 6 };

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{ false };

  /*! \brief Use don't cares for optimization. */
  bool use_dont_cares{ false };

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{ false };

  /*! \brief Show progress. */
  bool progress{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for refactoring.
 *
 * The data structure `refactoring_stats` provides data collected by running
 * `refactoring`.
 */
struct refactoring_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Accumulated runtime for computing MFFCs. */
  stopwatch<>::duration time_mffc{ 0 };

  /*! \brief Accumulated runtime for rewriting. */
  stopwatch<>::duration time_refactoring{ 0 };

  /*! \brief Accumulated runtime for simulating MFFCs. */
  stopwatch<>::duration time_simulation{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
    std::cout << fmt::format( "[i] MFFC time        = {:>5.2f} secs\n", to_seconds( time_mffc ) );
    std::cout << fmt::format( "[i] refactoring time = {:>5.2f} secs\n", to_seconds( time_refactoring ) );
    std::cout << fmt::format( "[i] simulation time  = {:>5.2f} secs\n", to_seconds( time_simulation ) );
  }
};

namespace detail
{

template<class Ntk, class RefactoringFn, class Iterator, class = void>
struct has_refactoring_with_dont_cares : std::false_type
{
};

template<class Ntk, class RefactoringFn, class Iterator>
struct has_refactoring_with_dont_cares<Ntk,
                                       RefactoringFn, Iterator,
                                       std::void_t<decltype( std::declval<RefactoringFn>()( std::declval<Ntk&>(),
                                                                                            std::declval<kitty::dynamic_truth_table>(),
                                                                                            std::declval<kitty::dynamic_truth_table>(),
                                                                                            std::declval<Iterator const&>(),
                                                                                            std::declval<Iterator const&>(),
                                                                                            std::declval<void( signal<Ntk> )>() ) )>> : std::true_type
{
};

template<class Ntk, class RefactoringFn, class Iterator>
inline constexpr bool has_refactoring_with_dont_cares_v = has_refactoring_with_dont_cares<Ntk, RefactoringFn, Iterator>::value;

template<class Ntk, class RefactoringFn, class NodeCostFn>
class refactoring_impl
{
public:
  refactoring_impl( Ntk& ntk, RefactoringFn&& refactoring_fn, refactoring_params const& ps, refactoring_stats& st, NodeCostFn const& cost_fn )
      : ntk( ntk ), refactoring_fn( refactoring_fn ), ps( ps ), st( st ), cost_fn( cost_fn ), required( ntk, UINT32_MAX )
  {
    if constexpr ( has_level_v<Ntk> )
    {
      ntk.events().release_add_event( add_event );
      ntk.events().release_modified_event( modified_event );
      ntk.events().release_delete_event( delete_event );
    }
  }

  ~refactoring_impl()
  {
    if constexpr ( has_level_v<Ntk> )
    {
      ntk.events().release_add_event( add_event );
      ntk.events().release_modified_event( modified_event );
      ntk.events().release_delete_event( delete_event );
    }
  }

  void run()
  {
    progress_bar pbar{ ntk.size(), "refactoring |{0}| node = {1:>4}   cand = {2:>4}   est. reduction = {3:>5}", ps.progress };

    stopwatch t( st.time_total );

    ntk.clear_visited();

    if ( ps.preserve_depth )
      compute_required();

    reconvergence_driven_cut_parameters rps;
    rps.max_leaves = ps.max_pis;
    reconvergence_driven_cut_statistics rst;
    detail::reconvergence_driven_cut_impl<Ntk, false, false> reconv_cuts( ntk, rps, rst );

    color_view<Ntk> color_ntk{ ntk };

    const auto size = ntk.size();
    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      // if ( i >= size )
      // {
      //   return false;
      // }
      if ( ntk.fanout_size( n ) == 0u )
      {
        return true;
      }

      /* update level for node */
      // if constexpr ( has_level_v<Ntk> )
      // {
      //   if ( ps.preserve_depth )
      //   {
      //     uint32_t level = 0;
      //     ntk.foreach_fanin( n, [&]( auto const& f ) {
      //       level = std::max( level, ntk.level( ntk.get_node( f ) ) );
      //     } );
      //     ntk.set_level( n, level + 1 );
      //   }
      // }

      const auto mffc = make_with_stopwatch<mffc_view<Ntk>>( st.time_mffc, ntk, n );

      pbar( i, i, _candidates, _estimated_gain );

      if ( mffc.num_pos() == 0 || mffc.size() < 4 )
      {
        return true;
      }

      kitty::dynamic_truth_table tt;
      std::vector<signal<Ntk>> leaves( ps.max_pis );
      uint32_t num_leaves = 0;

      if ( mffc.num_pis() <= ps.max_pis )
      {
        /* use MFFC */
        mffc.foreach_pi( [&]( auto const& m, auto j ) {
          leaves[j] = ntk.make_signal( m );
        } );

        num_leaves = mffc.num_pis();

        default_simulator<kitty::dynamic_truth_table> sim( mffc.num_pis() );
        tt = call_with_stopwatch( st.time_simulation,
                                  [&]() { return simulate<kitty::dynamic_truth_table>( mffc, sim )[0]; } );
      }
      else
      {
        /* compute a reconvergent-driven cut */
        std::vector<node<Ntk>> roots = { n };
        auto const extended_leaves = reconv_cuts.run( roots ).first;

        num_leaves = extended_leaves.size();
        assert( num_leaves <= ps.max_pis );

        for ( auto j = 0u; j < num_leaves; ++j )
        {
          leaves[j] = ntk.make_signal( extended_leaves[j] );
        }

        cut_view<Ntk> cut( ntk, extended_leaves, ntk.make_signal( n ) );
        default_simulator<kitty::dynamic_truth_table> sim( num_leaves );
        tt = call_with_stopwatch( st.time_simulation,
                                  [&]() { return simulate<kitty::dynamic_truth_table>( cut, sim )[0]; } );
      }

      /* mark cut boundaries */
      for ( auto j = 0u; j < num_leaves; ++j )
      {
        ntk.incr_fanout_size( ntk.get_node( leaves[j] ) );
      }

      signal<Ntk> new_f;
      bool resynthesized{ false };

      ntk.incr_trav_id();
      int32_t gain = recursive_deref_mark( n );

      /* unmark cut boundaries */
      for ( auto j = 0u; j < num_leaves; ++j )
      {
        ntk.decr_fanout_size( ntk.get_node( leaves[j] ) );
      }

      {
        if ( ps.use_dont_cares )
        {
          if constexpr ( has_refactoring_with_dont_cares_v<Ntk, RefactoringFn, decltype( leaves.begin() )> )
          {
            std::vector<node<Ntk>> pivots;
            for ( auto const& c : leaves )
            {
              pivots.push_back( ntk.get_node( c ) );
            }
            stopwatch t( st.time_refactoring );

            refactoring_fn( ntk, tt, satisfiability_dont_cares( ntk, pivots, 16u ), leaves.begin(), leaves.begin() + num_leaves, [&]( auto const& f ) { new_f = f; resynthesized = true; return false; } );
          }
          else
          {
            stopwatch t( st.time_refactoring );
            refactoring_fn( ntk, tt, leaves.begin(), leaves.begin() + num_leaves, [&]( auto const& f ) { new_f = f; resynthesized = true; return false; } );
          }
        }
        else
        {
          stopwatch t( st.time_refactoring );
          refactoring_fn( ntk, tt, leaves.begin(), leaves.begin() + num_leaves, [&]( auto const& f ) { new_f = f; resynthesized = true; return false; } );
        }
      }

     
      if ( !resynthesized || n == ntk.get_node( new_f ) )
      {
        /* mark cut boundaries */
        for ( auto j = 0u; j < num_leaves; ++j )
        {
          ntk.incr_fanout_size( ntk.get_node( leaves[j] ) );
        }

        recursive_ref( n );

        /* unmark cut boundaries */
        for ( auto j = 0u; j < num_leaves; ++j )
        {
          ntk.decr_fanout_size( ntk.get_node( leaves[j] ) );
        }

        return true;
      }

      /* ref only if it is a new node */
      if ( ntk.fanout_size( ntk.get_node( new_f ) ) == 0 )
      {
        /* mark cut boundaries */
        for ( auto j = 0u; j < num_leaves; ++j )
        {
          ntk.incr_fanout_size( ntk.get_node( leaves[j] ) );
        }

        recursive_deref_check_mark( ntk.get_node( new_f ) );
        gain -= recursive_ref( ntk.get_node( new_f ) );

        /* unmark cut boundaries */
        for ( auto j = 0u; j < num_leaves; ++j )
        {
          ntk.decr_fanout_size( ntk.get_node( leaves[j] ) );
        }
      }

      /* mark cut boundaries */
      for ( auto j = 0u; j < num_leaves; ++j )
      {
        ntk.incr_fanout_size( ntk.get_node( leaves[j] ) );
      }

      recursive_ref( n );

      /* unmark cut boundaries */
      for ( auto j = 0u; j < num_leaves; ++j )
      {
        ntk.decr_fanout_size( ntk.get_node( leaves[j] ) );
      }

      if constexpr ( has_level_v<Ntk> )
      {
        if ( ps.preserve_depth && ntk.level( ntk.get_node( new_f ) ) > ntk.required( n ) )
        {
          /* remove */
          if ( ntk.fanout_size( ntk.get_node( new_f ) ) == 0 )
            ntk.take_out_node( ntk.get_node( new_f ) );

          return true;
        }
      }

      if ( gain > 0 || ( ps.allow_zero_gain && gain == 0 ) )
      {
        ++_candidates;
        _estimated_gain += gain;
        ntk.substitute_node( n, new_f );
      }
      else
      {
        /* remove */
        if ( ntk.fanout_size( ntk.get_node( new_f ) ) == 0 )
          ntk.take_out_node( ntk.get_node( new_f ) );
      }

      if constexpr ( has_level_v<Ntk> )
      {
        if ( ps.preserve_depth )
          ntk.update_levels();
        // for ( auto i = 0; i < size; ++i )
        //   required[ntk.index_to_node( i )] = ntk.required( ntk.index_to_node( i ) );
        // ntk.update_required();
        // propagate_required_rec( ntk.node_to_index( n ), ntk.get_node( new_f ), size, required[n] );
      }
      return true;
    } );
  }

private:
  uint32_t recursive_deref_mark( node<Ntk> const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;
    
    ntk.set_visited( n, ntk.trav_id() );

    /* recursively collect nodes */
    uint32_t value{ cost_fn( ntk, n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref_mark( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  uint32_t recursive_deref_check_mark( node<Ntk> const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;
    
    if ( ntk.visited( n ) == ntk.trav_id() )
      // return recursive_deref_after_mark( n );
      return 0;

    /* recursively collect nodes */
    uint32_t value{ cost_fn( ntk, n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref_check_mark( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  uint32_t recursive_deref_after_mark( node<Ntk> const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    uint32_t value{ cost_fn( ntk, n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref_after_mark( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  uint32_t recursive_ref( node<Ntk> const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    uint32_t value{ cost_fn( ntk, n ) };
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  void compute_required()
  {
    if constexpr ( has_level_v<Ntk> )
    {
      ntk.foreach_po( [&]( auto const& f ) {
        required[f] = ntk.depth();
      } );

      for ( uint32_t index = ntk.size() - 1; index > ntk.num_pis(); index-- )
      {
        node<Ntk> n = ntk.index_to_node( index );
        uint32_t req = required[n];

        ntk.foreach_fanin( n, [&]( auto const& f ) {
          required[f] = std::min( required[f], req - 1 );
        } );
      }
    }
  }

  void propagate_required_rec( uint32_t root, node<Ntk> const& n, uint32_t size, uint32_t req )
  {
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return;

    /* recursively update required time */
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto const g = ntk.get_node( f );
      
      /* recur if it is still a node to explore and to update */
      if ( ntk.node_to_index( g ) > root && ( ntk.node_to_index( g ) > size || required[g] > req ) )
        propagate_required_rec( root, g, size, req - 1 );
      
      /* update the required time */
      if ( ntk.node_to_index( g ) < size )
        required[g] = std::min( required[g], req - 1 );
    } );
  }

private:
  void register_events()
  {
    if constexpr ( has_level_v<Ntk> )
    {
      auto const update_level_of_new_node = [&]( const auto& n ) {
        ntk.resize_levels();
        update_node_level( n );
      };

      auto const update_level_of_existing_node = [&]( node<Ntk> const& n, const auto& old_children ) {
        (void)old_children;
        ntk.resize_levels();
        update_node_level( n );
      };

      auto const update_level_of_deleted_node = [&]( node<Ntk> const& n ) {
        ntk.set_level( n, -1 );
      };

      // add_event = ntk.events().register_add_event( update_level_of_new_node );
      modified_event = ntk.events().register_modified_event( update_level_of_existing_node );
      delete_event = ntk.events().register_delete_event( update_level_of_deleted_node );
    }
  }

  /* maybe should move to depth_view */
  void update_node_level( node<Ntk> const& n, bool top_most = true )
  {
    if constexpr ( has_level_v<Ntk> )
    {
      uint32_t curr_level = ntk.level( n );

      uint32_t max_level = 0;
      ntk.foreach_fanin( n, [&]( const auto& f ) {
        auto const p = ntk.get_node( f );
        auto const fanin_level = ntk.level( p );
        if ( fanin_level > max_level )
        {
          max_level = fanin_level;
        }
      } );
      ++max_level;

      if ( curr_level != max_level )
      {
        ntk.set_level( n, max_level );

        /* update only one more level */
        if ( top_most )
        {
          ntk.foreach_fanout( n, [&]( const auto& p ) {
            update_node_level( p, false );
          } );
        }
      }
    }
  }

private:
  Ntk& ntk;
  RefactoringFn&& refactoring_fn;
  refactoring_params const& ps;
  refactoring_stats& st;
  NodeCostFn cost_fn;

  node_map<uint32_t, Ntk> required;

  uint32_t _candidates{ 0 };
  uint32_t _estimated_gain{ 0 };

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;
};

} /* namespace detail */

/*! \brief Boolean refactoring.
 *
 * This algorithm performs refactoring by collapsing maximal fanout-free cones
 * (MFFCs) into truth tables and recreating a new network structure from it.
 * If the MFFC is too large a reconvergent-driven cut is extracted.
 * The algorithm performs changes directly in the input network and keeps the
 * substituted structures dangling in the network.  They can be cleaned up using
 * the `cleanup_dangling` algorithm.
 *
 * The refactoring function must be of type `NtkDest::signal(NtkDest&,
 * kitty::dynamic_truth_table const&, LeavesIterator, LeavesIterator)` where
 * `LeavesIterator` can be dereferenced to a `NtkDest::signal`.  The last two
 * parameters compose an iterator pair where the distance matches the number of
 * variables of the truth table that is passed as second parameter.  There are
 * some refactoring algorithms in the folder
 * `mockturtle/algorithms/node_resyntesis`, since the resynthesis functions
 * have the same signature.
 *
 * **Required network functions:**
 * - `get_node`
 * - `size`
 * - `make_signal`
 * - `foreach_gate`
 * - `substitute_node`
 * - `clear_visited`
 * - `clear_values`
 * - `fanout_size`
 * - `set_value`
 * - `foreach_node`
 *
 * \param ntk Input network (will be changed in-place)
 * \param refactoring_fn Refactoring function
 * \param ps Refactoring params
 * \param pst Refactoring statistics
 * \param cost_fn Node cost function (a functor with signature `uint32_t(Ntk const&, node<Ntk> const&)`)
 */
template<class Ntk, class RefactoringFn, class NodeCostFn = unit_cost<Ntk>>
void refactoring( Ntk& ntk, RefactoringFn&& refactoring_fn, refactoring_params const& ps = {}, refactoring_stats* pst = nullptr, NodeCostFn const& cost_fn = {} )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the substitute_node method" );
  static_assert( has_clear_visited_v<Ntk>, "Ntk does not implement the clear_visited method" );
  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );

  refactoring_stats st;

  if ( ps.preserve_depth )
  {
    depth_view<Ntk> d_ntk{ ntk };
    fanout_view<depth_view<Ntk>> f_ntk{ d_ntk };

    detail::refactoring_impl<fanout_view<depth_view<Ntk>>, RefactoringFn, NodeCostFn> p( f_ntk, refactoring_fn, ps, st, cost_fn );
    p.run();
  }
  else
  {
    fanout_view<Ntk> f_ntk{ ntk };

    detail::refactoring_impl<fanout_view<Ntk>, RefactoringFn, NodeCostFn> p( f_ntk, refactoring_fn, ps, st, cost_fn );
    p.run();
  }

  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */