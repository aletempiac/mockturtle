/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
  \file lut_mapping.hpp
  \brief LUT mapping

  \author Mathias Soeken
*/

#pragma once

#include <cstdint>
#include <limits>

#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../utils/node_map.hpp"
#include "../views/topo_view.hpp"
#include "../views/depth_view.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"

namespace mockturtle
{

/*! \brief Parameters for lut_mapping.
 *
 * The data structure `lut_mapping_params` holds configurable parameters
 * with default arguments for `lut_mapping`.
 */
struct mapping_params
{
  mapping_params()
  {
    cut_enumeration_ps.cut_size = 4;
    cut_enumeration_ps.cut_limit = 12;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut size is 6, the default cut limit is 8.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Number of rounds for area flow optimization.
   *
   * The first round is used for delay optimization.
   */
  uint32_t rounds{2u};


  /*! \brief Number of rounds for exact area optimization. */
  uint32_t rounds_ela{1u};

  /*! \brief Be verbose. */
  bool verbose{false};
};

// template<typename Ntk>
// struct mapping_storage<Ntk>
// {
//   node_map<std::vector<signal<Ntk>>, Ntk> best_gates( ntk );
//   uint32_t mapping_size{0};
// };

/*! \brief Statistics for mapper.
 *
 * The data structure `mapper_stats` provides data collected by running
 * `mapper`.
 */
struct mapping_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  void report() const
  {
    std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

/* function to update all cuts after cut enumeration */
template<typename CutData>
struct mapping_update_cuts
{
  template<typename NetworkCuts, typename Ntk>
  static void apply( NetworkCuts const& cuts, Ntk const& ntk )
  {
    (void)cuts;
    (void)ntk;
  }
};

namespace detail
{

template<class NtkDest, class Ntk, class RewritingFn, typename CutData>
class mapping_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, true, CutData>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  mapping_impl( Ntk& ntk, RewritingFn const& rewriting_fn, mapping_params const& ps, mapping_stats& st )
      : ntk( ntk ),
        rewriting_fn( rewriting_fn ),
        ps( ps ),
        st( st ),
        flow_refs( ntk.size() ),
        map_refs( ntk.size(), 0 ),
        flows( ntk.size() ),
        arrivals( ntk.size() ),
        delays( ntk.size() ),
        cuts( cut_enumeration<Ntk, true, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    mapping_update_cuts<CutData>().apply( cuts, ntk );
  }

  NtkDest run()
  {
    stopwatch t( st.time_total );

    node_map<signal<NtkDest>, Ntk> old2new( ntk );
    NtkDest res;
    old2new[ntk.get_constant( false )] = res.get_constant( false );
    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      old2new[ntk.get_constant( true )] = res.get_constant( true );
    }
    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[n] = res.create_pi();
    } );

    /* compute and save topological order */
    top_order.reserve( ntk.size() );
    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      top_order.push_back( n );
    } );

    // node_map<std::vector<signal<Ntk>>, Ntk> best_gates( ntk );

    init_nodes();

    depth_view<NtkDest> res_depth{res};
    compute_mapping_delay( res_depth, old2new );

    // set_mapping_refs<true>();


    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? res.create_not( old2new[f] ) : old2new[f] );
    } );

    NtkDest ret = cleanup_dangling<NtkDest>( res );
    // set_mapping_refs<false>();
    // //print_state();

    // while ( iteration < ps.rounds )
    // {
    //   compute_mapping<false>();
    // }

    // while ( iteration < ps.rounds + ps.rounds_ela )
    // {
    //   compute_mapping<true>();
    // }

    // derive_mapping();

    return ret;
  }

private:
  uint32_t cut_area( cut_t const& cut ) const
  {
    return static_cast<uint32_t>( cut->data.cost );
  }

  void init_nodes()
  {
    ntk.foreach_node( [this]( auto n, auto ) {
      const auto index = ntk.node_to_index( n );

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
        arrivals[index] = 0.0f;
      }
      else
      {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }

      flows[index] = cuts.cuts( index )[0]->data.flow;
      delays[index] = cuts.cuts( index )[0]->data.delay;
    } );
  }

  void compute_mapping_delay( depth_view<NtkDest>& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      /* foreach cut */
      signal<NtkDest> best_signal;
      uint32_t best_arrival = UINT32_MAX;
      uint32_t best_size = UINT32_MAX;
      uint8_t best_cut = 0u;
      float best_area_flow = std::numeric_limits<float>::max();
      uint8_t cut_index = 0u;
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        /* trivial cuts */
        if ( cut->size() == 1 )
        {
          continue;
        }

        auto cut_flow = cut_leaves_flow( *cut );

        const auto tt = cuts.truth_table( *cut );
        assert( cut->size() == static_cast<unsigned>( tt.num_vars() ) );

        std::vector<signal<NtkDest>> children( cut->size() );

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          children[ctr++] = old2new[ntk.index_to_node( l )];
        }

        const auto on_signal = [&]( auto const& f_new ) {
          auto area_flow = cut_flow + recursive_ref<NtkDest>( res, res.get_node( f_new ) );
          recursive_deref<NtkDest>( res, res.get_node( f_new ) );

          if ( compare_map_delay( res.level( res.get_node( f_new ) ), best_arrival, area_flow, best_area_flow, cut->size(), best_size ) )
          {
            best_signal = f_new;
            best_arrival = res.level( res.get_node( f_new ) );
            best_area_flow = area_flow;
            best_size = cut->size();
            best_cut = cut_index;
          }
          return true;
        };

        rewriting_fn( res, cuts.truth_table( *cut ), children.begin(), children.end(), on_signal );
        cut_index++;
      }
      old2new[n] = best_signal;
      recursive_ref<NtkDest>( res, res.get_node( best_signal ) );
      flows[ntk.index_to_node( n )] = best_area_flow / flow_refs[ntk.index_to_node( n )];
      arrivals[ntk.index_to_node( n )] = best_arrival;
      if ( best_cut != 0 )
      {
        cuts.cuts( ntk.index_to_node( n ) ).update_best( best_cut );
      }
    }
    // set_mapping_refs<ELA>();
    //print_state();
  }

  template<bool ELA>
  void compute_mapping()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      compute_best_cut<ELA>( ntk.node_to_index( n ) );
    }
    set_mapping_refs<ELA>();
    //print_state();
  }

  template<bool ELA>
  void set_mapping_refs()
  {
    const auto coef = 1.0f / ( 1.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    /* compute current delay and update mapping refs */
    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, delays[index] );

      if constexpr ( !ELA )
      {
        map_refs[index]++;
      }
    } );

    /* compute current area and update mapping refs */
    area = 0;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      /* skip constants and PIs (TODO: stop earlier) */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( map_refs[index] == 0 )
        continue;

      if constexpr ( !ELA )
      {
        for ( auto leaf : cuts.cuts( index )[0] )
        {
          map_refs[leaf]++;
        }
      }
      area++;
    }

    /* blend flow referenes */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      flow_refs[i] = coef * flow_refs[i] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( map_refs[i] ) );
    }

    ++iteration;
  }

  std::pair<float, uint32_t> cut_flow( cut_t const& cut )
  {
    uint32_t time{0u};
    float flow{0.0f};

    for ( auto leaf : cut )
    {
      time = std::max( time, delays[leaf] );
      flow += flows[leaf];
    }

    return {flow + cut_area( cut ), time + 1u};
  }

  inline float cut_leaves_flow( cut_t const& cut )
  {
    float flow{0.0f};

    for ( auto leaf : cut )
    {
      flow += flows[leaf];
    }

    return flow;
  }

  /* reference cut:
   *   adds cut to current mapping and recursively adds best cuts of leaf
   *   nodes, if they are not part of the current mapping.
   */
  uint32_t cut_ref( cut_t const& cut )
  {
    uint32_t count = cut_area( cut );
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( map_refs[leaf]++ == 0 )
      {
        count += cut_ref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  /* dereference cut:
   *   removes cut from current mapping and recursively removes best cuts of
   *   leaf nodes, if they are part of the current mapping.
   *   (this is the inverse operation to cut_ref)
   */
  uint32_t cut_deref( cut_t const& cut )
  {
    uint32_t count = cut_area( cut );
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( --map_refs[leaf] == 0 )
      {
        count += cut_deref( cuts.cuts( leaf ).best() );
      }
    }
    return count;
  }

  /* reference cut (special version):
   *   this special version of cut_ref does two additional things:
   *   1. it stops recursing if it has found `limit` cuts
   *   2. it remembers all cuts for which the reference count increases in the
   *      vector `tmp_area`.
   */
  uint32_t cut_ref_limit_save( cut_t const& cut, uint32_t limit )
  {
    uint32_t count = cut_area( cut );
    if ( limit == 0 )
      return count;

    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      tmp_area.push_back( leaf );
      if ( map_refs[leaf]++ == 0 )
      {
        count += cut_ref_limit_save( cuts.cuts( leaf ).best(), limit - 1 );
      }
    }
    return count;
  }

  /* estimates the cost of adding this cut to the mapping:
   *   This algorithm references cuts recursively to estimate how many cuts
   *   would be needed to add to the mapping if `cut` were to be added.  It
   *   temporarily modifies the reference counters but reverts them eventually.
   */
  uint32_t cut_area_estimation( cut_t const& cut )
  {
    tmp_area.clear();
    const auto count = cut_ref_limit_save( cut, 8 );
    for ( auto const& n : tmp_area )
    {
      map_refs[n]--;
    }
    return count;
  }

  template<bool ELA>
  void compute_best_cut( uint32_t index )
  {
    constexpr auto mf_eps{0.005f};

    float flow;
    uint32_t time{0};
    int32_t best_cut{-1};
    float best_flow{std::numeric_limits<float>::max()};
    uint32_t best_time{std::numeric_limits<uint32_t>::max()};
    int32_t cut_index{-1};

    if constexpr ( ELA )
    {
      if ( map_refs[index] > 0 )
      {
        cut_deref( cuts.cuts( index )[0] );
      }
    }

    for ( auto* cut : cuts.cuts( index ) )
    {
      ++cut_index;
      if ( cut->size() == 1 )
        continue;

      if constexpr ( ELA )
      {
        flow = static_cast<float>( cut_area_estimation( *cut ) );
      }
      else
      {
        std::tie( flow, time ) = cut_flow( *cut );
      }

      if ( best_cut == -1 || best_flow > flow + mf_eps || ( best_flow > flow - mf_eps && best_time > time ) )
      {
        best_cut = cut_index;
        best_flow = flow;
        best_time = time;
      }
    }

    if ( best_cut == -1 )
      return;

    if constexpr ( ELA )
    {
      if ( map_refs[index] > 0 )
      {
        cut_ref( cuts.cuts( index )[best_cut] );
      }
    }
    else
    {
      map_refs[index] = 0;
    }
    if constexpr ( ELA )
    {
      best_time = cut_flow( cuts.cuts( index )[best_cut] ).second;
    }
    delays[index] = best_time;
    flows[index] = best_flow / flow_refs[index];

    if ( best_cut != 0 )
    {
      cuts.cuts( index ).update_best( best_cut );
    }
  }

  void derive_mapping()
  {
    ntk.clear_mapping();

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      const auto index = ntk.node_to_index( n );
      if ( map_refs[index] == 0 )
        continue;

      std::vector<node<Ntk>> nodes;
      for ( auto const& l : cuts.cuts( index ).best() )
      {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
    }
  }

  inline bool compare_map_delay( uint32_t arrival, uint32_t best_arrival, float area_flow, float best_area_flow, uint32_t size, uint32_t best_size )
  {
    if ( arrival < best_arrival )
    {
      return true;
    }
    else if ( area_flow < best_area_flow - fEpsilon )
    {
      return true;
    }
    else if ( area_flow > best_area_flow + fEpsilon )
    {
      return false;
    }
    else if ( size < best_size )
    {
      return true;
    }
    /* TODO: add compare on fanout size */
    return false;
  }

  void print_state()
  {
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      std::cout << fmt::format( "*** Obj = {:>3} (node = {:>3})  FlowRefs = {:5.2f}  MapRefs = {:>2}  Flow = {:5.2f}  Delay = {:>3}\n", i, ntk.index_to_node( i ), flow_refs[i], map_refs[i], flows[i], delays[i] );
      //std::cout << cuts.cuts( i );
    }
    std::cout << fmt::format( "Level = {}  Area = {}\n", delay, area );
  }

private:
  Ntk& ntk;
  RewritingFn const& rewriting_fn;
  mapping_params const& ps;
  mapping_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  uint32_t delay{0};     /* current delay of the mapping */
  uint32_t area{0};      /* current area of the mapping */
  const float fEpsilon{0.005f};
  //bool ela{false};       /* compute exact area */

  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<float> arrivals;
  std::vector<uint32_t> delays;
  network_cuts_t cuts;

  std::vector<uint32_t> tmp_area; /* temporary vector to compute exact area */
};

}; /* namespace detail */

/*! \brief LUT mapping.
 *
 * This function implements a LUT mapping algorithm.  It is controlled by two
 * template arguments `StoreFunction` (defaulted to `true`) and `CutData`
 * (defaulted to `cut_enumeration_mf_cut`).  The first argument `StoreFunction`
 * controls whether the LUT function is stored in the mapping.  In that case
 * truth tables are computed during cut enumeration, which requires more
 * runtime.  The second argument is simuilar to the `CutData` argument in
 * `cut_enumeration`, which can specialize the cost function to select priority
 * cuts and store additional data.  For LUT mapping using this function the
 * type passed as `CutData` must implement the following three fields:
 *
 * - `uint32_t delay`
 * - `float flow`
 * - `float costs`
 *
 * See `include/mockturtle/algorithms/cut_enumeration/mf_cut.hpp` for one
 * example of a CutData type that implements the cost function that is used in
 * the LUT mapper `&mf` in ABC.
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
 * - `clear_mapping`
 * - `add_to_mapping`
 * - `set_lut_funtion` (if `StoreFunction` is true)
 *
   \verbatim embed:rst

   .. note::

      The implementation of this algorithm was heavily inspired but the LUT
      mapping command ``&mf`` in ABC.
   \endverbatim
 */
template<class Ntk, class NtkDest = Ntk, class RewritingFn, typename CutData = cut_enumeration_mf_cut>
NtkDest mapping( Ntk& ntk, RewritingFn const& rewriting_fn = {}, mapping_params const& ps = {}, mapping_stats* pst = nullptr )
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
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  // static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  // static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );

  mapping_stats st;
  detail::mapping_impl<NtkDest, Ntk, RewritingFn, CutData> p( ntk, rewriting_fn, ps, st );
  auto res = p.run();
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }

  return res;
}

} /* namespace mockturtle */
