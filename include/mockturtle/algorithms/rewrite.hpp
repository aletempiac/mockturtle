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
  \brief Rewrite factored form literals version. TODO: merge with standard rewrite

  \author Alessandro Tempia Calvino
*/

#pragma once

#include "../traits.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/node_map.hpp"
#include "../utils/progress_bar.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "cleanup.hpp"
#include "cut_enumeration/rewrite_cut.hpp"
#include "cut_enumeration.hpp"

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

  /*! \brief If true, candidates are only accepted if they do not increase logic depth. */
  bool preserve_depth{false};

  /*! \brief Allow rewrite with multiple structures */
  bool allow_multiple_structures{true};

  /*! \brief Allow zero-gain substitutions */
  bool allow_zero_gain{false};

  /*! \brief Allow zero-gain substitutions */
  bool aggressive_zero_gain{false};

  /*! \brief Optimize literal cost instead of number of nodes */
  bool optimize_literal_cost{false};

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
      : ntk( ntk ), library( library ), ps( ps ), st( st ), cost_fn( cost_fn ), required( ntk, UINT32_MAX )
  {
    register_events();
  }
  
  ~rewrite_impl()
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
    progress_bar pbar{ntk.size(), "rewrite |{0}| node = {1:>4}   cand = {2:>4}   est. reduction = {3:>5}", ps.progress};

    stopwatch t( st.time_total );

    ntk.incr_trav_id();

    if ( ps.optimize_literal_cost )
    {
      /* set node value for POs */
      ntk.clear_values();
      ntk.foreach_po( [&]( auto const& f ) {
        ntk.incr_value( ntk.get_node( f ) );
      } );
    }

    if ( ps.preserve_depth )
    {
      compute_required();
    }

    /* initialize cut manager */
    cut_enumeration_stats cst;
    network_cuts_t cuts( ps.use_mffc ? 0 : ntk.size() + ( ntk.size() >> 1 ) );
    cut_manager_t cut_manager( ntk, ps.cut_enumeration_ps, cst, cuts );

    /* initialize cuts for constant nodes and PIs */
    if ( !ps.use_mffc )
      cut_manager.init_cuts();

    auto& db = library.get_database();

    const auto size = ntk.size();
    ntk.foreach_gate( [&]( auto const& n, auto i ) {
      if ( ntk.fanout_size( n ) == 0u )
      {
        return true;
      }

      pbar( i, i, _candidates, _estimated_gain );

      int32_t best_gain = -1;
      int32_t best_gain2 = -1;
      uint32_t best_level = UINT32_MAX;
      uint32_t best_cut = 0;
      signal<Ntk> best_signal;
      std::vector<signal<Ntk>> best_leaves;
      bool best_phase = false;
      std::vector<signal<Ntk>> leaves( num_vars, ntk.get_constant( false ) );

      /* update level for node */
      if constexpr ( has_level_v<Ntk> )
      {
        if ( ps.preserve_depth )
        {
          uint32_t level = 0;
          ntk.foreach_fanin( n, [&]( auto const& f ) {
            level = std::max( level, ntk.level( ntk.get_node( f ) ) );
          } );
          ntk.set_level( n, level + cost_fn( ntk, n ) );
          best_level = level + cost_fn( ntk, n );
        }
      }

      {
        /* use cuts */
        cut_manager.clear_cuts( n );
        cut_manager.compute_cuts( n );

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
              return library.get_supergates( tt_npn );
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
            int32_t mffc_size, num_lits;
            if ( ps.optimize_literal_cost )
            {
              std::tie( mffc_size, num_lits ) = measure_literals_dereference( n, cut );
            }
            else
            {
              mffc_size = measure_mffc_deref( n, cut );
            }

            for ( auto const& dag : *structures )
            {
              int32_t gain, gain2;

              if ( ps.optimize_literal_cost )
              {
                auto [nodes_added, lits_added, level] = evaluate_entry_literals( db.get_node( dag.root ), leaves, ntk.fanout_size( n ) == 1u );
                gain = num_lits - lits_added;
                gain2 = mffc_size - nodes_added;

                  /* discard if dag.root and n are the same */
                if ( ntk.node_to_index( n ) == db.value( db.get_node( dag.root ) ) )
                  continue;

                /* discard if level increases */
                if constexpr ( has_level_v<Ntk> )
                {
                  if ( ps.preserve_depth && level > required[n] )
                    continue;
                }

                /* discard if no gain or gain2 */
                if ( gain < 0 || ( gain == 0 && gain2 < 0 ) || ( !ps.allow_zero_gain && gain == 0 && gain2 >= 0 ) )
                  continue;

                if ( ( gain > best_gain ) || ( gain == best_gain && gain2 > best_gain2 ) || ( gain == best_gain && gain2 == best_gain2 && level < best_level ) )
                {
                  ++_candidates;
                  best_gain = gain;
                  best_gain2 = gain2;
                  best_signal = dag.root;
                  best_leaves = leaves;
                  best_phase = phase;
                  best_cut = cut_index;
                  best_level = level;
                }
                else if ( ps.allow_zero_gain && ps.aggressive_zero_gain && ( gain > 0 || ( gain2 > 0 && gain == 0 ) || ( ps.allow_zero_gain && gain == 0 ) ) && ( ( gain > best_gain ) || ( gain == best_gain && gain2 > best_gain2 ) ) )
                {
                  ++_candidates;
                  best_gain = gain;
                  best_gain2 = gain2;
                  best_signal = dag.root;
                  best_leaves = leaves;
                  best_phase = phase;
                  best_cut = cut_index;
                  best_level = level;
                }
              }
              else
              {
                auto [nodes_added, level] = evaluate_entry( db.get_node( dag.root ), leaves );
                gain = mffc_size - nodes_added;

                /* discard if dag.root and n are the same */
                if ( ntk.node_to_index( n ) == db.value( db.get_node( dag.root ) ) )
                  continue;

                /* discard if no gain */
                if ( gain < 0 || ( !ps.allow_zero_gain && gain == 0 ) )
                  continue;

                /* discard if level increases */
                if constexpr ( has_level_v<Ntk> )
                {
                  if ( ps.preserve_depth && level > required[n] )
                    continue;
                }

                if ( ( gain > best_gain ) || ( gain == best_gain && level < best_level ) )
                {
                  ++_candidates;
                  best_gain = gain;
                  best_signal = dag.root;
                  best_leaves = leaves;
                  best_phase = phase;
                  best_cut = cut_index;
                  best_level = level;
                }
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

      if ( ( best_gain > 0 || ( best_gain == 0 && best_gain2 > 0 ) ) || ( ps.allow_zero_gain && best_gain == 0 ) )
      {
        /* replace node wth the new structure */
        topo_view topo{db, best_signal};
        auto new_f = cleanup_dangling( topo, ntk, best_leaves.begin(), best_leaves.end() ).front();

        assert( n != ntk.get_node( new_f ) );

        if ( ps.optimize_literal_cost && ntk.value( n ) )
        {
          /* inherit the PO info */
          ntk.set_value( ntk.get_node( new_f ), ntk.value( n ) );
        }

        _estimated_gain += best_gain;
        ntk.substitute_node( n, new_f ^ best_phase );

        if constexpr ( has_level_v<Ntk> )
        {
          /* TODO: propagate new required to leaves */
          if ( ps.preserve_depth )
          {
            propagate_required_rec( ntk.node_to_index( n ), ntk.get_node( new_f ), size, required[n] );
            assert( ntk.level( ntk.get_node( new_f ) ) <= required[n] );
          }
        }

        clear_cuts_fanout_rec( cuts, cut_manager, ntk.get_node( new_f ) );
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
      ntk.incr_fanout_size( ntk.index_to_node( leaf ) );
    }

    int32_t mffc_size = static_cast<int32_t>( recursive_ref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_fanout_size( ntk.index_to_node( leaf ) );
    }

    return mffc_size;
  }

  int32_t measure_mffc_deref( node<Ntk> const& n,  cut_t const* cut )
  {
    /* reference cut leaves */
    for ( auto leaf : *cut )
    {
      ntk.incr_fanout_size( ntk.index_to_node( leaf ) );
    }

    int32_t mffc_size = static_cast<int32_t>( recursive_deref( n ) );

    /* dereference leaves */
    for ( auto leaf : *cut )
    {
      ntk.decr_fanout_size( ntk.index_to_node( leaf ) );
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
      if ( ntk.decr_fanout_size( ntk.get_node( s ) ) == 0 )
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
      if ( ntk.incr_fanout_size( ntk.get_node( s ) ) == 0 )
      {
        value += recursive_ref( ntk.get_node( s ) );
      }
    } );
    return value;
  }

  std::pair<int32_t, uint32_t> measure_literals_dereference( node<Ntk> const& n,  cut_t const* cut )
  {
    ntk.incr_trav_id();
    uint32_t ref_leaves = 0;

    /* reference cut leaves */
    uint32_t index = 0;
    for ( auto leaf : *cut )
    {
      if ( !ntk.is_pi( ntk.index_to_node( leaf ) ) && ntk.fanout_size( ntk.index_to_node( leaf ) ) + ntk.value( ntk.index_to_node( leaf ) ) == 1 )
        ref_leaves |= 1 << index;

      ntk.incr_fanout_size( ntk.index_to_node( leaf ) );
      ++index;
    }

    uint32_t mffc_size = 0;
    int32_t lits = static_cast<int32_t>( measure_literals_dereference_rec( n, mffc_size ) );

    /* dereference leaves */
    index = 0;
    for ( auto leaf : *cut )
    {
      ntk.decr_fanout_size( ntk.index_to_node( leaf ) );

      if ( ntk.value( leaf ) == 0 )
      {
        if ( ntk.fanout_size( leaf ) == 0 && ( ( ref_leaves >> index ) & 1 ) == 1 )
          lits -= 2;
        if ( ntk.fanout_size( leaf ) == 0 && ( ( ref_leaves >> index ) & 1 ) != 1 )
          --lits;
        else if ( !ntk.is_pi( ntk.index_to_node( leaf ) ) && ntk.fanout_size( leaf ) == 1 && ( ( ref_leaves >> index ) & 1 ) != 1 )
          ++lits;
      }
      ++index;
    }

    return {mffc_size, lits};
  }

  uint32_t measure_literals_dereference_rec( node<Ntk> const& n, uint32_t& mffc_size )
  {
    /* terminate? */
    if ( ntk.is_constant( n ) )
      return 0;
    
    if ( ntk.is_pi( n ) )
      return 1u;
    
    mffc_size += cost_fn( ntk, n );

    /* recursively dereference and count literals */
    uint32_t lits = 0;
    ntk.foreach_fanin( n, [&]( auto const& s ) {
      auto g = ntk.get_node( s );
      if ( ntk.is_constant( g ) )
      {
        ntk.decr_fanout_size( g );
        return;
      }
      if ( ntk.is_pi( g ) )
      {
        ntk.decr_fanout_size( g );
        ++lits;
        return;
      }

      auto ref = ntk.decr_fanout_size( g );
      if ( ref == 0 )
      {
        /* visited before, this counts as a literal */
        // if ( ntk.visited( g ) == ntk.trav_id() )
        //   ++lits;
        lits += measure_literals_dereference_rec( g, mffc_size );
      }
      else
      {
        /* add literal */
        ++lits;
        if ( ref + ntk.value( g ) == 1 )
          ++lits;
        ntk.set_visited( g, ntk.trav_id() );
      }
    } );
    return lits;
  }

  inline std::pair<int32_t, uint32_t> evaluate_entry( node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    auto& db = library.get_database();
    db.incr_trav_id();

    return evaluate_entry_rec( n, leaves );
  }

  inline std::tuple<int32_t, int32_t, uint32_t> evaluate_entry_literals( node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves, bool single_fanout_root )
  {
    auto& db = library.get_database();
    db.incr_trav_id();
    auto [gates_added, level] = evaluate_entry_rec( n, leaves );

    /* is const */
    if ( db.is_constant( n ) )
    {
      return {0, 0, 0};
    }

    /* is leaf */
    // if ( db.is_pi( n ) )
    // {
    //   auto leaf = ntk.get_node( leaves[n - 1] );
    //   cost = 0;
    //   /* add cost if hashed leaf becomes a literal */
    //   if ( ntk.fanout_size( leaf ) == 1 && !ntk.is_pi( leaf ) )
    //     ++cost;
    //   if ( single_fanout_root )
    //     ++cost;
    //   return cost;
    // }

    int32_t cost = 0;
    uint32_t ref_leaves = 0;
    uint32_t index = 0;
    for ( auto leaf : leaves )
    {
      if ( !ntk.is_pi( ntk.get_node( leaf ) ) && ntk.fanout_size( ntk.get_node( leaf ) ) + ntk.value( ntk.get_node( leaf ) ) == 1 )
        ref_leaves |= 1 << index;

      ++index;
    }

    db.incr_trav_id();
    entry_reference_rec( n, leaves );

    index = 0;
    for ( auto leaf : leaves )
    {
      if ( ntk.fanout_size( ntk.get_node( leaf ) ) > 1 && ( ( ref_leaves >> index ) & 1 ) == 1 )
        ++cost;

      ++index;
    }

    ntk.incr_trav_id();
    db.incr_trav_id();

    /* hashed, do evaluate */
    if ( db.visited( n ) < ( db.trav_id() - 1 ) && db.value( n ) < ntk.size() )
    {
      cost = 0;
      /* add cost if hashed node becomes a literal */
      if ( ntk.fanout_size( db.value( n ) ) + ntk.value( db.value( n ) ) == 1 && !ntk.is_pi( db.value( n ) ) )
        ++cost;
      /* add cost if root node becomes a literal */
      if ( single_fanout_root )
        ++cost;
      return {gates_added, cost, level};
    }

    cost += evaluate_entry_literals_dereference_rec( n, leaves );

    return {gates_added, cost, level};
  }

  std::pair<int32_t, uint32_t> evaluate_entry_rec( node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    auto& db = library.get_database();
    if ( db.is_pi( n ) || db.is_constant( n ) )
      return {0, 0};
    if ( db.visited( n ) == db.trav_id() )
      return {0, 0};

    db.set_visited( n, db.trav_id() );

    int32_t area = 0;
    uint32_t level = 0;
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
        if constexpr ( has_level_v<Ntk> )
        {
          level = std::max( level, ntk.level( ntk.get_node( leaves[db.node_to_index( g ) - 1] ) ) );
        }
        return;
      }

      auto [area_rec, level_rec] = evaluate_entry_rec( g, leaves );
      area += area_rec;
      level = std::max( level, level_rec );

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
        return {area + ( ntk.fanout_size( *val ) > 0 ? 0 : cost_fn( ntk, n ) ), level + cost_fn( ntk, n )}; /* TODO: have two different cost functions for size and depth */
      }
    }

    db.set_value( n, ntk.size() );
    return {area + cost_fn( ntk, n ), level + cost_fn( ntk, n )};
  }

  void entry_reference_rec( node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    auto& db = library.get_database();
    /* terminate? */
    if ( db.is_constant( n ) )
      return;

    if ( db.is_pi( n ) )
    {
      /* register associated leaf */
      db.set_value( n, ntk.node_to_index( ntk.get_node( leaves[db.node_to_index( n ) - 1] ) ) );
      return;
    }

    /* hashed, do not recur */
    if ( db.visited( n ) != db.trav_id() && db.value( n ) < ntk.size() )
    {
      if ( ntk.fanout_size( ntk.index_to_node( db.value( n ) ) ) == 0 )
      {
        /* remove hash info */
        db.set_value( n, 0 );
        db.set_visited( n, db.trav_id() );
      }
      else
      {
        return;
      }
    }
    // else if ( db.visited( n ) != db.trav_id() )
    // {
    //   /* reset and use as a reference counter */
    //   db.set_value( n, 0u );
    //   db.set_visited( n, db.trav_id() );
    // }

    /* recursively reference */
    db.foreach_fanin( n, [&]( auto const& f ) {
      auto g = db.get_node( f );
      if ( db.is_constant( g ) )
      {
        ntk.incr_fanout_size( g );
        return;
      }
      if ( db.is_pi( g ) )
      {
        auto leaf = leaves[db.node_to_index( g ) - 1];
        ntk.incr_fanout_size( ntk.get_node( leaf ) );
        return;
      }

      /* hashed */
      if ( db.visited( g ) != db.trav_id() && db.value( g ) < ntk.size() )
      {
        /* remove hash info if not referenced */
        if ( ntk.fanout_size( ntk.index_to_node( db.value( g ) ) ) == 0 )
        {
          db.set_value( g, 0 );
          db.set_visited( g, db.trav_id() );
        }
        else
        {
          ntk.incr_fanout_size( ntk.index_to_node( db.value( g ) ) );
          return;
        }
      }
      else if ( db.visited( g ) != db.trav_id() )
      {
        /* reset and use as a reference counter */
        db.set_value( g, 0u );
        db.set_visited( g, db.trav_id() );
      }

      if ( db.incr_value( g ) == 0 )
      {
        entry_reference_rec( g, leaves );
      }
    } );
  }

  int32_t evaluate_entry_literals_dereference_rec( node<Ntk> const& n, std::vector<signal<Ntk>> const& leaves )
  {
    auto& db = library.get_database();

    /* terminate? */
    if ( db.is_constant( n ) )
      return 0;

    /* TODO: export it out of the recusion */
    if ( db.is_pi( n ) )
    {
      auto leaf = ntk.get_node( leaves[db.node_to_index( n ) - 1] );
      if ( ntk.fanout_size( leaf ) > 0 || ntk.is_pi( leaf ) )
        return 1;
      else
        return 0;
    }

    /* hashed, do not recur */
    // if ( db.visited( n ) < ( db.trav_id() - 1 ) && db.value( n ) < ntk.size() )
    // {
    //   return 0;
    // }

    /* recursively dereference and count literals */
    int32_t lits = 0;
    db.foreach_fanin( n, [&]( auto const& f ) {
      auto g = db.get_node( f );
      if ( db.is_constant( g ) )
      {
        ntk.decr_fanout_size( g );
        return;
      }
      if ( db.is_pi( g ) )
      {
        auto leaf = leaves[db.node_to_index( g ) - 1];
        auto fanout_leaf = ntk.decr_fanout_size( ntk.get_node( leaf ) );
        if ( fanout_leaf != 0 || ntk.visited( ntk.get_node( leaf ) ) == ntk.trav_id() || ntk.is_pi( ntk.get_node( leaf ) ) )
        {
          ntk.set_visited( ntk.get_node( leaf ), ntk.trav_id() );
          ++lits;
        }
        return;
      }

      /* hashed */
      if ( db.visited( g ) < ( db.trav_id() - 1 ) )
      {
        assert( db.value( g ) < ntk.size() );
        auto fanout_hashed = ntk.decr_fanout_size( ntk.index_to_node( db.value( g ) ) );
        ++lits;
        /* hashed node had fanout of one -> is a new literal */
        if ( fanout_hashed == 1 && ntk.value( ntk.index_to_node( db.value( g ) ) ) == 0 )
          ++lits;
        return;
      }

      auto ref = db.decr_value( g );
      if ( ref != 0 || db.visited( g ) == db.trav_id() )
      {
        db.set_visited( g, db.trav_id() );
        ++lits;
      }

      if ( ref == 0 )
      {
        lits += evaluate_entry_literals_dereference_rec( g, leaves );
      }
    } );

    return lits;
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
  Library&& library;
  rewrite_params const& ps;
  rewrite_stats& st;
  NodeCostFn cost_fn;

  node_map<uint32_t, Ntk> required;

  uint32_t _candidates{0};
  uint32_t _estimated_gain{0};

  /* events */
  std::shared_ptr<typename network_events<Ntk>::add_event_type> add_event;
  std::shared_ptr<typename network_events<Ntk>::modified_event_type> modified_event;
  std::shared_ptr<typename network_events<Ntk>::delete_event_type> delete_event;
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

  if ( ps.preserve_depth )
  {
    using depth_view_t = depth_view<Ntk, NodeCostFn>;
    depth_view_t depth_ntk{ ntk };
    using fanout_view_t = fanout_view<depth_view_t>;
    fanout_view_t fanout_view{ depth_ntk };

    detail::rewrite_impl<fanout_view_t, Library, NodeCostFn> p( fanout_view, library, ps, st, cost_fn );
    p.run();
  }
  else
  {
    using fanout_view_t = fanout_view<Ntk>;
    fanout_view_t fanout_view{ ntk };

    detail::rewrite_impl<fanout_view_t, Library, NodeCostFn> p( fanout_view, library, ps, st, cost_fn );
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