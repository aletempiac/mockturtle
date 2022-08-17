/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file rewrite.hpp
  \brief Rewrite

  \author Alessandro Tempia Calvino
*/
#pragma once

#include "../traits.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/cut_view.hpp"
#include "../views/fanout_view.hpp"
#include "../views/mffc_view.hpp"
#include "cleanup.hpp"
#include "cut_enumeration/rewrite_cut.hpp"
#include "cut_enumeration.hpp"
#include "detail/mffc_utils.hpp"
#include "dont_cares.hpp"
#include "simulation.hpp"

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>

namespace mockturtle
{

/*! \brief Parameters for Rewrite.
 *
 * The data structure `rewrite_params` holds configurable parameters with
 * default arguments for `rewrite`.
 */
struct rewrite_params
{
  rewrite_params()
  {
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Cut enumeration parameters. */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Rewrite using MFFC instead of cuts */
  bool use_mffc{false};

  /*! \brief If true, candidates are only accepted if they do not increase logic level of node. */
  bool preserve_depth{false};

  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{true};

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{false};

  /*! \brief Use don't cares for optimization. */
  bool use_dont_cares{false};

  /*! \brief Show progress. */
  bool progress{false};

  /*! \brief Be verbose. */
  bool verbose{false};
};

/*! \brief Statistics for rewrite.
 *
 * The data structure `rewrite_stats` provides data collected by running
 * `rewrite`.
 */
struct rewrite_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  /*! \brief Accumulated runtime for computing MFFCs. */
  stopwatch<>::duration time_cuts{0};

  /*! \brief Accumulated runtime for rewrite. */
  stopwatch<>::duration time_matching{0};

  /*! \brief Accumulated runtime for rewrite. */
  stopwatch<>::duration time_rewrite{0};

  /*! \brief Accumulated runtime for simulating MFFCs. */
  stopwatch<>::duration time_simulation{0};

  /*! \brief Expected gain. */
  uint32_t estimated_gain{0};

  /*! \brief Candidates */
  uint32_t candidates{0};

  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
    std::cout << fmt::format( "[i] cuts time        = {:>5.2f} secs\n", to_seconds( time_cuts ) );
    std::cout << fmt::format( "[i] matching time    = {:>5.2f} secs\n", to_seconds( time_matching ) );
    std::cout << fmt::format( "[i] rewrite time     = {:>5.2f} secs\n", to_seconds( time_rewrite ) );
    std::cout << fmt::format( "[i] simulation time  = {:>5.2f} secs\n", to_seconds( time_simulation ) );
  }
};

namespace detail
{

template<class Ntk, class Library, class NodeCostFn>
class rewrite_impl
{
  static constexpr uint32_t num_vars = 4u;
  using network_cuts_t = dynamic_network_cuts<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_manager_t = detail::dynamic_cut_enumeration_impl<Ntk, num_vars, true, cut_enumeration_rewrite_cut>;
  using cut_t = typename network_cuts_t::cut_t;
  using node_data = typename Ntk::storage::element_type::node_type;

public:
  rewrite_impl( Ntk& ntk, Library&& library, rewrite_params const& ps, rewrite_stats& st, NodeCostFn const& cost_fn )
      : ntk( ntk ), library( library ), ps( ps ), st( st ), cost_fn( cost_fn ) {}

  void run()
  {
    progress_bar pbar{ntk.size(), "rewrite |{0}| node = {1:>4}   cand = {2:>4}   est. reduction = {3:>5}", ps.progress};

    stopwatch t( st.time_total );

    /* for cost estimation we use reference counters initialized by the fanout size */
    initialize_values_with_fanout( ntk );
    ntk.incr_trav_id();

    /* initialize cut manager */
    cut_enumeration_stats cst;
    network_cuts_t cuts( ps.use_mffc ? 0 : ntk.size() + ( ntk.size() >> 1 ) );
    cut_manager_t cut_manager( ntk, ps.cut_enumeration_ps, cst, cuts );

    /* initialize cuts for constant nodes and PIs */
    if ( !ps.use_mffc )
      cut_manager.init_cuts();

    auto& db = library.get_database();

    const auto size = ntk.num_gates();
    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      if ( ntk.fanout_size( n ) == 0u )
      {
        return true;
      }

      pbar( i, i, _candidates, _estimated_gain );

      int32_t best_gain = -1;
      uint32_t best_cut = 0;
      signal<Ntk> best_signal;
      std::vector<signal<Ntk>> best_leaves;
      bool best_phase;
      std::vector<signal<Ntk>> leaves( num_vars, ntk.get_constant( false ) );

      if ( ps.use_mffc )
      {
        const auto mffc = make_with_stopwatch<mffc_view<Ntk>>( st.time_cuts, ntk, n );

        if ( mffc.num_pos() == 0 || mffc.num_pis() > num_vars || mffc.size() < num_vars + 1 )
        {
          return true;
        }

        /* disconnect n */
        ntk.set_value( n, 0 );

        default_simulator<kitty::static_truth_table<num_vars>> sim;
        const auto tt = call_with_stopwatch( st.time_simulation,
                                            [&]() { return simulate<kitty::static_truth_table<num_vars>>( mffc, sim )[0]; } );

        auto config = kitty::exact_npn_canonization( tt );
        auto tt_npn = std::get<0>( config );
        auto neg = std::get<1>( config );
        auto perm = std::get<2>( config );

        auto const structures = call_with_stopwatch( st.time_matching, [&]() {
            if ( ps.use_dont_cares )
            {
              std::vector<node<Ntk>> pivots;
              mffc.foreach_pi( [&]( auto const& m, auto j ) {
                pivots.push_back( m );
              } );

              const auto sdc = satisfiability_dont_cares<Ntk, num_vars>( ntk, pivots, 12u );
              const auto dc_npn = apply_npn_transformation( sdc, neg & ~( 1 << num_vars ), perm );

              return library.get_supergates( tt_npn, dc_npn, neg, perm );
            }
            else
            {
              return library.get_supergates( tt_npn );
            }
          }
        );

        if ( structures == nullptr )
        {
          return true;
        }

        /* dereference n */
        int32_t mffc_size = recursive_deref( n );

        uint32_t negation = 0;
        std::array<uint8_t, num_vars> permutation;

        for ( auto j = 0u; j < num_vars; ++j )
        {
          permutation[perm[j]] = j;
          negation |= ( ( neg >> perm[j] ) & 1 ) << j;
        }

        /* save output negation to apply */
        bool phase = ( neg >> num_vars == 1 ) ? true : false;

        mffc.foreach_pi( [&]( auto const& m, auto j ) {
          leaves[permutation[j]] = ntk.make_signal( m );
        } );

        for ( auto j = 0u; j < num_vars; ++j )
        {
          if ( ( negation >> j ) & 1 )
          {
            leaves[j] = !leaves[j];
          }
        }

        {
          stopwatch t( st.time_rewrite );

          for ( auto const& dag : *structures )
          {
            int32_t nodes_added = evaluate_entry( db, db.get_node( dag.root ), leaves );

            /* discard if dag.root and n are the same */
            if ( ntk.node_to_index( n ) == db.value( db.get_node( dag.root ) ) )
              continue;
            
            int32_t gain = mffc_size - nodes_added;

            if ( ( gain > 0 || ( ps.allow_zero_gain && gain == 0 ) ) && gain > best_gain )
            {
              ++_candidates;
              best_gain = gain;
              best_phase = phase;
              best_signal = dag.root;
              best_leaves = leaves;
            }

            if ( !ps.allow_multiple_structures )
              break;
          }
        }
      }
      else
      {
        /* use cuts */
        cut_manager.clear_cuts( n );
        cut_manager.compute_cuts( n );

        /* disconnect n */
        ntk.set_value( n, 0 );

        uint32_t cut_index = 0;
        for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
        {
          /* skip trivial cut */
          if ( ( cut->size() == 1 && *cut->begin() == ntk.node_to_index( n ) ) )
          {
            ++cut_index;
            continue;
          }

          /* Boolean matching */
          auto config = kitty::exact_npn_canonization( cuts.truth_table( *cut ) );
          auto tt_npn = std::get<0>( config );
          auto neg = std::get<1>( config );
          auto perm = std::get<2>( config );

          auto const structures = call_with_stopwatch( st.time_matching, [&]() {
              if ( ps.use_dont_cares )
              {
                std::vector<node<Ntk>> pivots;
                for ( auto leaf : *cut )
                {
                  pivots.push_back( ntk.index_to_node( leaf ) );
                }

                const auto sdc = satisfiability_dont_cares<Ntk, num_vars>( ntk, pivots, 8u );
                const auto dc_npn = apply_npn_transformation( sdc, neg & ~( 1 << num_vars ), perm );

                return library.get_supergates( tt_npn, dc_npn, neg, perm );
              }
              else
              {
                return library.get_supergates( tt_npn );
              }
            }
          );

          if ( structures == nullptr )
          {
            ++cut_index;
            continue;
          }

          uint32_t negation = 0;
          std::array<uint8_t, num_vars> permutation;

          for ( auto j = 0u; j < num_vars; ++j )
          {
            permutation[perm[j]] = j;
            negation |= ( ( neg >> perm[j] ) & 1 ) << j;
          }

          /* save output negation to apply */
          bool phase = ( neg >> num_vars == 1 ) ? true : false;

          {
            auto j = 0u;
            for ( auto const leaf : *cut )
            {
              leaves[permutation[j++]] = ntk.make_signal( ntk.index_to_node( leaf ) );
            }
          }

          for ( auto j = 0u; j < num_vars; ++j )
          {
            if ( ( negation >> j ) & 1 )
            {
              leaves[j] = !leaves[j];
            }
          }

          {
            stopwatch t( st.time_rewrite );

            /* measure the MFFC contained in the cut */
            int32_t mffc_size = measure_mffc_deref( n, cut );

            for ( auto const& dag : *structures )
            {
              int32_t nodes_added = evaluate_entry( db, db.get_node( dag.root ), leaves );

              /* discard if dag.root and n are the same */
              if ( ntk.node_to_index( n ) == db.value( db.get_node( dag.root ) ) )
                continue;

              int32_t gain = mffc_size - nodes_added;

              if ( ( gain > 0 || ( ps.allow_zero_gain && gain == 0 ) ) && gain > best_gain )
              {
                ++_candidates;
                best_gain = gain;
                best_signal = dag.root;
                best_leaves = leaves;
                best_phase = phase;
                best_cut = cut_index;
              }

              if ( !ps.allow_multiple_structures )
                break;
            }

            /* restore contained MFFC */
            measure_mffc_ref( n, cut );
            ++cut_index;

            if ( cut->size() == 0 || ( cut->size() == 1 && *cut->begin() != ntk.node_to_index( n ) ) )
              break;
          }
        }
      }

      /* reconnect n */
      ntk.set_value( n, ntk.fanout_size( n ) );

      if ( best_gain > 0 || ( ps.allow_zero_gain && best_gain == 0 ) )
      {
        if ( !ps.use_mffc )
        {
          measure_mffc_deref( n, &cuts.cuts( ntk.node_to_index( n ) )[best_cut] );
          clear_cuts_fanout_rec( cuts, cut_manager, n );
        }

        topo_view topo{db, best_signal};
        auto new_f = cleanup_dangling( topo, ntk, best_leaves.begin(), best_leaves.end() ).front();

        assert( n != ntk.get_node( new_f ) );
        recursive_ref( ntk.get_node( new_f ) );

        _estimated_gain += best_gain;
        ntk.substitute_node( n, new_f ^ best_phase );
        ntk.set_value( n, 0 );
        ntk.set_value( ntk.get_node( new_f ), ntk.fanout_size( ntk.get_node( new_f ) ) );
      }
      else if ( ps.use_mffc )
      {
        recursive_ref( n );
      }

      return true;
    } );

    st.estimated_gain = _estimated_gain;
    st.candidates = _candidates;
  }

private:
  int32_t measure_mffc_ref( node<Ntk> const& n,  cut_t const* cut )
  {
    /* reference cut leaves */
    for ( auto leaf : *cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    int32_t mffc_size = static_cast<int32_t>( recursive_ref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }

    return mffc_size;
  }

  int32_t measure_mffc_deref( node<Ntk> const& n,  cut_t const* cut )
  {
    /* reference cut leaves */
    for ( auto leaf : *cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    int32_t mffc_size = static_cast<int32_t>( recursive_deref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }

    return mffc_size;
  }

  uint32_t recursive_deref( node<Ntk> const& n )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      return 0;

    /* recursively collect nodes */
    uint32_t value{cost_fn( ntk, n )};
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.decr_value( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_deref( ntk.get_node( s ) );
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
    uint32_t value{cost_fn( ntk, n )};
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      if ( ntk.incr_value( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  inline int32_t evaluate_entry( auto& db, node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    db.incr_trav_id();
    return evaluate_entry_rec( db, n, leaves );
  }

  int32_t evaluate_entry_rec( auto& db, node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    if ( db.is_pi( n ) || db.is_constant( n ) )
      return 0;
    if ( db.visited( n ) == db.trav_id() )
      return 0;

    db.set_visited( n, db.trav_id() );

    int32_t area = 0;
    bool hashed = true;

    std::array<signal<Ntk>, Ntk::max_fanin_size> node_data;
    db.foreach_fanin( n, [&]( auto const& f, auto i ) {
      node<Ntk> g = db.get_node( f );
      if ( db.is_constant( g ) )
      {
        node_data[i] = f; /* ntk.get_costant( db.is_complemented( f ) ) */
        return;
      }
      if ( db.is_pi( g ) )
      {
        node_data[i] = leaves[db.node_to_index( g ) - 1] ^ db.is_complemented( f );
        return;
      }

      area += evaluate_entry_rec( db, g, leaves );

      /* check value */
      if ( db.value( g ) < ntk.size() )
      {
        node_data[i] = ntk.make_signal( ntk.index_to_node( db.value( g ) ) ) ^ db.is_complemented( f );
      }
      else
      {
        hashed = false;
      }
    } );

    if ( hashed )
    {
      /* try hash */
      /* only AIG is supported now */
      auto val = ntk.has_and( node_data[0], node_data[1] );

      if ( val.has_value() )
      {
        db.set_value( n, *val );
        return area + ( ntk.value( *val ) > 0 ? 0 : cost_fn( ntk, n ) );
      }
    }

    db.set_value( n, ntk.size() );
    return area + cost_fn( ntk, n );
  }

  void clear_cuts_fanout_rec( network_cuts_t& cuts, cut_manager_t& cut_manager, node<Ntk> const& n )
  {
    ntk.foreach_fanout( n, [&]( auto const& g ) {
      auto const index = ntk.node_to_index( g );
      if ( cuts.cuts( index ).size() > 0 )
      {
        cut_manager.clear_cuts( g );
        clear_cuts_fanout_rec( cuts, cut_manager, g );
      }
    } );
  }

private:
  Ntk& ntk;
  Library&& library;
  rewrite_params const& ps;
  rewrite_stats& st;
  NodeCostFn cost_fn;

  uint32_t _candidates{0};
  uint32_t _estimated_gain{0};
};

} /* namespace detail */

/*! \brief Boolean rewrite.
 *
 * This algorithm rewrites maximal fanout-free cones (MFFCs) or enumerated cuts
 * using new network structure from a database.
 * The algorithm performs changes directly in the input network and keeps the
 * substituted structures dangling in the network. They can be cleaned up using
 * the `cleanup_dangling` algorithm.
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
 * \param library Exact library containing pre-computed structures
 * \param ps Rewrite params
 * \param pst Rewrite statistics
 * \param cost_fn Node cost function (a functor with signature `uint32_t(Ntk const&, node<Ntk> const&)`)
 */
template<class Ntk, class Library, class NodeCostFn = unit_cost<Ntk>>
void rewrite( Ntk& ntk, Library&& library, rewrite_params const& ps = {}, rewrite_stats* pst = nullptr, NodeCostFn const& cost_fn = {} )
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

  rewrite_stats st;

  using fanout_view_t = fanout_view<Ntk>;
  fanout_view_t fanout_view{ntk};

  detail::rewrite_impl<fanout_view_t, Library, NodeCostFn> p( fanout_view, library, ps, st, cost_fn );
  p.run();
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