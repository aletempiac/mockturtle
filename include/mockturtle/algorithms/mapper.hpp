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
  \file mapper.hpp
  \brief Mapper

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <limits>

#include <fmt/format.h>

#include "../utils/stopwatch.hpp"
#include "../utils/node_map.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/tech_library.hpp"
#include "../views/topo_view.hpp"
#include "../views/depth_view.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"
#include "cut_enumeration/map_cut.hpp"

namespace mockturtle
{

/*! \brief Parameters for lut_mapping.
 *
 * The data structure `lut_mapping_params` holds configurable parameters
 * with default arguments for `lut_mapping`.
 */
struct map_params
{
  map_params()
  {
    cut_enumeration_ps.cut_size = 4;
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut size is 4, the default cut limit is 8.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Required time for delay optimization. */
  float required_time{0.0f};

  /*! \brief Do area optimization. */
  bool skip_delay_round{false};

  /*! \brief Number of rounds for area flow optimization. */
  uint32_t area_flow_rounds{1u};

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t ela_rounds{1u};

  /*! \brief Use structural choices. */
  bool choices{false};

  /*! \brief Be verbose. */
  bool verbose{false};
};

/*! \brief Statistics for mapper.
 *
 * The data structure `mapper_stats` provides data collected by running
 * `mapper`.
 */
struct map_stats
{
  /* \brief Area and delay */
  float area{0};
  float delay{0};
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{0};

  void report() const
  {
    std::cout << fmt::format( "[i] area = {:>5.2f}; delay = {:>5.2f}\n", area, delay );
    std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

/* function to update all cuts after cut enumeration */
template<typename CutData>
struct map_update_cuts
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

template<class NtkDest, class Ntk, class RewritingFn, class DelayCostFn, class AreaCostFn, typename CutData>
class map_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, true, CutData>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  map_impl( Ntk& ntk, RewritingFn const& rewriting_fn, map_params const& ps, map_stats& st )
      : ntk( ntk ),
        rewriting_fn( rewriting_fn ),
        ps( ps ),
        st( st ),
        flow_refs( ntk.size() ),
        map_refs( ntk.size(), 0u ),
        flows( ntk.size() ),
        arrival( ntk.size() ),
        required( ntk.size() ),
        cut_area( ntk.size() ),
        cuts( cut_enumeration<Ntk, true, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    map_update_cuts<CutData>().apply( cuts, ntk );
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

    init_nodes();

    depth_view<NtkDest, DelayCostFn> res_depth{res};
    /* compute mapping delay */
    compute_mapping<false>( res_depth, old2new );

    /* compute mapping global area flow */
    while ( iteration < ps.area_flow_rounds + 1 )
    {
      compute_required_time( res_depth, old2new );
      compute_mapping<true>( res_depth, old2new );
    }

    /* compute mapping local exact area */
    while ( iteration < ps.ela_rounds + ps.area_flow_rounds + 1 )
    {
      compute_required_time( res_depth, old2new );
      compute_exact_area( res_depth, old2new );
    }

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? res.create_not( old2new[f] ) : old2new[f] );
    } );

    std::cout << "LUTs: " << int( area ) << std::endl;

    return cleanup_dangling<NtkDest>( res );
  }

private:
  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n, auto ) {
      const auto index = ntk.node_to_index( n );

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
        arrival[index] = 0.0f;
      }
      else
      {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }

      flows[index] = cuts.cuts( index )[0]->data.flow;
    } );
  }

  template<bool DO_AREA>
  void compute_mapping( depth_view<NtkDest>& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    if constexpr ( DO_AREA )
    {
      for ( auto i = 0u; i < res.size(); i++ )
      {
        res.set_value( res.index_to_node( i ), 0u );
      }
    }

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      signal<NtkDest> best_signal;
      uint32_t best_arrival = UINT32_MAX;
      uint32_t best_size = UINT32_MAX;
      uint8_t best_cut = 0u;
      float best_area = 0.0f;
      float best_area_flow = std::numeric_limits<float>::max();
      uint8_t cut_index = 0u;

      /* foreach cut */
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        /* trivial cuts */
        if ( cut->size() == 1 )
          continue;

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
          auto area_local = static_cast<float>( recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) ) );
          auto area_flow = cut_flow + area_local;
          recursive_deref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );

          if constexpr ( DO_AREA )
          {
            if ( res.level( res.get_node( f_new ) ) > required[ntk.node_to_index( n )] )
              return true;
          }

          if ( compare_map<DO_AREA>( res.level( res.get_node( f_new ) ), best_arrival, area_flow, best_area_flow, cut->size(), best_size ) )
          {
            best_signal = f_new;
            best_arrival = res.level( res.get_node( f_new ) );
            best_area_flow = area_flow;
            best_size = cut->size();
            best_cut = cut_index;
            best_area = area_local;
          }
          return true;
        };

        rewriting_fn( res, cuts.truth_table( *cut ), children.begin(), children.end(), on_signal );
        cut_index++;
      }
      old2new[n] = best_signal;
      recursive_ref<NtkDest>( res, res.get_node( best_signal ) );
      flows[ntk.index_to_node( n )] = best_area_flow / flow_refs[ntk.index_to_node( n )];
      arrival[ntk.index_to_node( n )] = best_arrival;
      cut_area[ntk.index_to_node( n )] = best_area;
      if ( best_cut != 0 )
      {
        cuts.cuts( ntk.index_to_node( n ) ).update_best( best_cut );
      }
    }
    set_mapping_refs();
  }


  void compute_exact_area( depth_view<NtkDest>& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    for ( auto i = 0u; i < res.size(); i++ )
    {
      res.set_value( res.index_to_node( i ), 0u );
    }

    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      
      signal<NtkDest> best_signal;
      uint32_t best_arrival = UINT32_MAX;
      uint32_t best_size = UINT32_MAX;
      uint8_t best_cut = 0u;
      float best_area = std::numeric_limits<float>::max();
      uint8_t cut_index = 0u;

      if ( map_refs[ntk.node_to_index( n )] )
      {
        cut_deref( cuts.cuts( ntk.node_to_index( n ) ).best(), n );
      }
      /* foreach cut */
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        /* trivial cuts */
        if ( cut->size() == 1 )
          continue;

        const auto tt = cuts.truth_table( *cut );
        assert( cut->size() == static_cast<unsigned>( tt.num_vars() ) );

        std::vector<signal<NtkDest>> children( cut->size() );

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          children[ctr++] = old2new[ntk.index_to_node( l )];
        }

        const auto on_signal = [&]( auto const& f_new ) {
          cut_area[n] = recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );
          recursive_deref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );
          auto area_local = cut_ref( *cut, n );
          cut_deref( *cut, n );

          if ( res.level( res.get_node( f_new ) ) > required[ntk.node_to_index( n )] )
            return true;

          if ( compare_map<true>( res.level( res.get_node( f_new ) ), best_arrival, area_local, best_area, cut->size(), best_size ) )
          {
            best_signal = f_new;
            best_arrival = res.level( res.get_node( f_new ) );
            best_area = area_local;
            best_size = cut->size();
            best_cut = cut_index;
          }
          return true;
        };

        rewriting_fn( res, cuts.truth_table( *cut ), children.begin(), children.end(), on_signal );
        cut_index++;
      }
      old2new[n] = best_signal;
      recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( best_signal ) );
      arrival[ntk.index_to_node( n )] = best_arrival;
      cut_area[ntk.index_to_node( n )] = best_area;
      if ( best_cut != 0 )
      {
        cuts.cuts( ntk.index_to_node( n ) ).update_best( best_cut );
      }
      if ( map_refs[ntk.node_to_index( n )] )
      {
        cut_ref( cuts.cuts( ntk.node_to_index( n ) ).best(), n );
      }
    }
    set_mapping_refs();
  }

  void set_mapping_refs()
  {
    const auto coef = 1.0f / ( 2.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    std::fill( map_refs.begin(), map_refs.end(), 0u );

    /* compute current delay and update mapping refs */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, arrival[index] );
      map_refs[index]++;
    } );

    /* compute current area and update mapping refs */
    area = 0.0f;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( map_refs[index] == 0 )
        continue;

      for ( auto leaf : cuts.cuts( index ).best() )
      {
        map_refs[leaf]++;
      }

      // area += cut_area[index];
      area += 1.0f;
    }

    /* blend flow references */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      flow_refs[i] = coef * flow_refs[i] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( map_refs[i] ) );
    }

    ++iteration;
  }

  void compute_required_time( NtkDest& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    std::vector<float> required_res( res.size(), std::numeric_limits<float>::max() );

    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      const auto index_res = res.node_to_index( res.get_node( old2new[s] ) );
      if ( ps.required_time == 0.0f )
      {
        required[index] = delay;
        required_res[index_res] = delay;
      }
      else
      {
        required[index] = ps.required_time;
        required_res[index_res] = ps.required_time;
      }
    } );

    auto i = res.size();
    while ( i-- > 0u )
    {
      const auto n = res.index_to_node( i );
      if ( res.is_pi( n ) || res.is_constant( n ) )
        break;
      const auto cost = static_cast<float>( DelayCostFn{}( res, n ) );
      res.foreach_fanin( n, [&]( const auto& child ) {
        const auto child_index = res.node_to_index( res.get_node( child ) );
        required_res[child_index] = std::min( required_res[child_index], required_res[i] - cost );
      } );
    }

    ntk.foreach_node( [&]( auto const& n ) {
      const auto index = ntk.node_to_index( n );
      required[index] = required_res[res.node_to_index( res.get_node( old2new[n] ) )];
    } );
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

  float cut_ref( cut_t const& cut, node<Ntk> const& n )
  {
    float count = cut_area[n];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( map_refs[leaf]++ == 0u )
      {
        count += cut_ref( cuts.cuts( leaf ).best(), ntk.index_to_node( leaf ) );
      }
    }
    return count;
  }

  float cut_deref( cut_t const& cut, node<Ntk> const& n )
  {
    float count = cut_area[n];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( --map_refs[leaf] == 0u )
      {
        count += cut_deref( cuts.cuts( leaf ).best(), ntk.index_to_node( leaf ) );
      }
    }
    return count;
  }

  template<bool DO_AREA>
  inline bool compare_map( uint32_t arrival, uint32_t best_arrival, float area_flow, float best_area_flow, uint32_t size, uint32_t best_size )
  {
    if constexpr ( DO_AREA )
    {
      if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
      else if ( arrival < best_arrival )
      {
        return true;
      }
      else if ( arrival > best_arrival )
      {
        return false;
      }
    }
    else
    {
      if ( arrival < best_arrival )
      {
        return true;
      }
      else if ( arrival > best_arrival )
      {
        return false;
      }
      else if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
    }
    if ( size < best_size )
    {
      return true;
    }
    /* TODO: add compare on fanout size */
    return false;
  }


private:
  Ntk& ntk;
  RewritingFn const& rewriting_fn;
  map_params const& ps;
  map_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  float delay{0.0f};     /* current delay of the mapping */
  float area{0.0f};      /* current area of the mapping */
  const float epsilon{0.005f}; /* epsilon */

  std::vector<node<Ntk>> top_order;
  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<float> arrival;
  std::vector<float> required;
  std::vector<float> cut_area;
  network_cuts_t cuts;
};

} /* namespace detail */

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
template<class Ntk, class NtkDest = Ntk, class RewritingFn, class DelayCostFn = unit_cost<NtkDest>, class AreaCostFn = unit_cost<NtkDest>, typename CutData = cut_enumeration_mf_cut>
NtkDest map( Ntk& ntk, RewritingFn const& rewriting_fn = {}, map_params const& ps = {}, map_stats* pst = nullptr )
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

  map_stats st;
  detail::map_impl<NtkDest, Ntk, RewritingFn, DelayCostFn, AreaCostFn, CutData> p( ntk, rewriting_fn, ps, st );
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



namespace detail
{

template<class NtkDest, class Ntk, class RewritingFn, class DelayCostFn, class AreaCostFn, typename CutData>
class map_choices_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, true, CutData>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  map_choices_impl( Ntk& ntk, RewritingFn const& rewriting_fn, map_params const& ps, map_stats& st )
      : ntk( ntk ),
        rewriting_fn( rewriting_fn ),
        ps( ps ),
        st( st ),
        flow_refs( ntk.size() ),
        map_refs( ntk.size(), 0u ),
        flows( ntk.size() ),
        arrival( ntk.size() ),
        required( ntk.size() ),
        cut_area( ntk.size() ),
        cuts( cut_enumeration_choices<Ntk, true, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    map_update_cuts<CutData>().apply( cuts, ntk );
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

    init_nodes();

    depth_view<NtkDest, DelayCostFn> res_depth{res};
    /* compute mapping delay */
    compute_mapping<false>( res_depth, old2new );

    /* compute mapping global area flow */
    while ( iteration < ps.area_flow_rounds + 1 )
    {
      compute_required_time( res_depth, old2new );
      compute_mapping<true>( res_depth, old2new );
    }

    /* compute mapping local exact area */
    while ( iteration < ps.ela_rounds + ps.area_flow_rounds + 1 )
    {
      compute_required_time( res_depth, old2new );
      compute_exact_area( res_depth, old2new );
    }

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? res.create_not( old2new[f] ) : old2new[f] );
    } );

    std::cout << "LUTs: " << int( area ) << std::endl;

    return cleanup_dangling<NtkDest>( res );
  }

private:
  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n, auto ) {
      const auto index = ntk.node_to_index( n );

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        flow_refs[index] = 1.0f;
        arrival[index] = 0.0f;
      }
      else
      {
        flow_refs[index] = static_cast<float>( ntk.fanout_size( n ) );
      }

      flows[index] = cuts.cuts( index )[0]->data.flow;
    } );
  }

  template<bool DO_AREA>
  void compute_mapping( depth_view<NtkDest>& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    if constexpr ( DO_AREA )
    {
      for ( auto i = 0u; i < res.size(); i++ )
      {
        res.set_value( res.index_to_node( i ), 0u );
      }
    }

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_ci( n ) || !ntk.is_choice_representative( n ) )
        return;

      signal<NtkDest> best_signal;
      uint32_t best_arrival = UINT32_MAX;
      uint32_t best_size = UINT32_MAX;
      uint8_t best_cut = 0u;
      float best_area = 0.0f;
      float best_area_flow = std::numeric_limits<float>::max();
      uint8_t cut_index = 0u;

      /* foreach cut */
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        /* trivial cuts */
        if ( cut->size() == 1 )
          continue;

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
          auto area_local = static_cast<float>( recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) ) );
          auto area_flow = cut_flow + area_local;
          recursive_deref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );

          if constexpr ( DO_AREA )
          {
            if ( res.level( res.get_node( f_new ) ) > required[ntk.node_to_index( n )] )
              return true;
          }

          if ( compare_map<DO_AREA>( res.level( res.get_node( f_new ) ), best_arrival, area_flow, best_area_flow, cut->size(), best_size ) )
          {
            best_signal = f_new;
            best_arrival = res.level( res.get_node( f_new ) );
            best_area_flow = area_flow;
            best_size = cut->size();
            best_cut = cut_index;
            best_area = area_local;
          }
          return true;
        };

        rewriting_fn( res, cuts.truth_table( *cut ), children.begin(), children.end(), on_signal );
        cut_index++;
      }
      old2new[n] = best_signal;
      recursive_ref<NtkDest>( res, res.get_node( best_signal ) );
      flows[ntk.index_to_node( n )] = best_area_flow / flow_refs[ntk.index_to_node( n )];
      arrival[ntk.index_to_node( n )] = best_arrival;
      cut_area[ntk.index_to_node( n )] = best_area;
      if ( best_cut != 0 )
      {
        cuts.cuts( ntk.index_to_node( n ) ).update_best( best_cut );
      }
    } );
    set_mapping_refs();
  }


  void compute_exact_area( depth_view<NtkDest>& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    for ( auto i = 0u; i < res.size(); i++ )
    {
      res.set_value( res.index_to_node( i ), 0u );
    }

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_ci( n ) || !ntk.is_choice_representative( n ) )
        return;
      
      signal<NtkDest> best_signal;
      uint32_t best_arrival = UINT32_MAX;
      uint32_t best_size = UINT32_MAX;
      uint8_t best_cut = 0u;
      float best_area = std::numeric_limits<float>::max();
      uint8_t cut_index = 0u;

      if ( map_refs[ntk.node_to_index( n )] )
      {
        cut_deref( cuts.cuts( ntk.node_to_index( n ) ).best(), n );
      }
      /* foreach cut */
      for ( auto& cut : cuts.cuts( ntk.node_to_index( n ) ) )
      {
        /* trivial cuts */
        if ( cut->size() == 1 )
          continue;

        const auto tt = cuts.truth_table( *cut );
        assert( cut->size() == static_cast<unsigned>( tt.num_vars() ) );

        std::vector<signal<NtkDest>> children( cut->size() );

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          children[ctr++] = old2new[ntk.index_to_node( l )];
        }

        const auto on_signal = [&]( auto const& f_new ) {
          cut_area[n] = recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );
          recursive_deref<NtkDest, AreaCostFn>( res, res.get_node( f_new ) );
          auto area_local = cut_ref( *cut, n );
          cut_deref( *cut, n );

          if ( res.level( res.get_node( f_new ) ) > required[ntk.node_to_index( n )] )
            return true;

          if ( compare_map<true>( res.level( res.get_node( f_new ) ), best_arrival, area_local, best_area, cut->size(), best_size ) )
          {
            best_signal = f_new;
            best_arrival = res.level( res.get_node( f_new ) );
            best_area = area_local;
            best_size = cut->size();
            best_cut = cut_index;
          }
          return true;
        };

        rewriting_fn( res, cuts.truth_table( *cut ), children.begin(), children.end(), on_signal );
        cut_index++;
      }
      old2new[n] = best_signal;
      recursive_ref<NtkDest, AreaCostFn>( res, res.get_node( best_signal ) );
      arrival[ntk.index_to_node( n )] = best_arrival;
      cut_area[ntk.index_to_node( n )] = best_area;
      if ( best_cut != 0 )
      {
        cuts.cuts( ntk.index_to_node( n ) ).update_best( best_cut );
      }
      if ( map_refs[ntk.node_to_index( n )] )
      {
        cut_ref( cuts.cuts( ntk.node_to_index( n ) ).best(), n );
      }
    } );
    set_mapping_refs();
  }

  void set_mapping_refs()
  {
    const auto coef = 1.0f / ( 2.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    std::fill( map_refs.begin(), map_refs.end(), 0u );

    /* compute current delay and update mapping refs */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, arrival[index] );
      map_refs[index]++;
    } );

    /* compute current area and update mapping refs */
    area = 0.0f;
    auto index = ntk.size();
    while( index-- > 0u )
    {
      const auto n = ntk.index_to_node( index );
      /* skip constants and PIs */
      if ( ntk.is_ci( n ) || !ntk.is_choice_representative( n ) )
        continue;

      if ( map_refs[index] == 0 )
        continue;

      for ( auto leaf : cuts.cuts( index ).best() )
      {
        map_refs[leaf]++;
      }

      // area += cut_area[index];
      area += 1.0f;
    }

    /* blend flow references */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      flow_refs[i] = coef * flow_refs[i] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( map_refs[i] ) );
    }

    ++iteration;
  }

  void compute_required_time( NtkDest& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    std::vector<float> required_res( res.size(), std::numeric_limits<float>::max() );

    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      const auto index_res = res.node_to_index( res.get_node( old2new[s] ) );
      if ( ps.required_time == 0.0f )
      {
        required[index] = delay;
        required_res[index_res] = delay;
      }
      else
      {
        required[index] = ps.required_time;
        required_res[index_res] = ps.required_time;
      }
    } );

    auto i = res.size();
    while ( i-- > 0u )
    {
      const auto n = res.index_to_node( i );
      if ( res.is_ci( n ) )
        continue;
      const auto cost = static_cast<float>( DelayCostFn{}( res, n ) );
      res.foreach_fanin( n, [&]( const auto& child ) {
        const auto child_index = res.node_to_index( res.get_node( child ) );
        required_res[child_index] = std::min( required_res[child_index], required_res[i] - cost );
      } );
    }

    ntk.foreach_node( [&]( auto const& n ) {
      const auto index = ntk.node_to_index( n );
      required[index] = required_res[res.node_to_index( res.get_node( old2new[n] ) )];
    } );
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

  float cut_ref( cut_t const& cut, node<Ntk> const& n )
  {
    float count = cut_area[n];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( map_refs[leaf]++ == 0u )
      {
        count += cut_ref( cuts.cuts( leaf ).best(), ntk.index_to_node( leaf ) );
      }
    }
    return count;
  }

  float cut_deref( cut_t const& cut, node<Ntk> const& n )
  {
    float count = cut_area[n];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( --map_refs[leaf] == 0u )
      {
        count += cut_deref( cuts.cuts( leaf ).best(), ntk.index_to_node( leaf ) );
      }
    }
    return count;
  }

  template<bool DO_AREA>
  inline bool compare_map( uint32_t arrival, uint32_t best_arrival, float area_flow, float best_area_flow, uint32_t size, uint32_t best_size )
  {
    if constexpr ( DO_AREA )
    {
      if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
      else if ( arrival < best_arrival )
      {
        return true;
      }
      else if ( arrival > best_arrival )
      {
        return false;
      }
    }
    else
    {
      if ( arrival < best_arrival )
      {
        return true;
      }
      else if ( arrival > best_arrival )
      {
        return false;
      }
      else if ( size < best_size )
      {
        return true;
      }
      else if ( size > best_size )
      {
        return false;
      }
      else if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
    }
    if ( size < best_size )
    {
      return true;
    }
    /* TODO: add compare on fanout size */
    return false;
  }


private:
  Ntk& ntk;
  RewritingFn const& rewriting_fn;
  map_params const& ps;
  map_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  float delay{0.0f};     /* current delay of the mapping */
  float area{0.0f};      /* current area of the mapping */
  const float epsilon{0.005f}; /* epsilon */

  std::vector<float> flow_refs;
  std::vector<uint32_t> map_refs;
  std::vector<float> flows;
  std::vector<float> arrival;
  std::vector<float> required;
  std::vector<float> cut_area;
  network_cuts_t cuts;
};

} /* namespace detail */


template<class Ntk, class NtkDest = Ntk, class RewritingFn, class DelayCostFn = unit_cost<NtkDest>, class AreaCostFn = unit_cost<NtkDest>, typename CutData = cut_enumeration_mf_cut>
NtkDest map_choices( Ntk& ntk, RewritingFn const& rewriting_fn = {}, map_params const& ps = {}, map_stats* pst = nullptr )
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
  static_assert( has_foreach_choice_v<Ntk>, "Ntk does not implement the foreach_choice method" );

  map_stats st;
  detail::map_choices_impl<NtkDest, Ntk, RewritingFn, DelayCostFn, AreaCostFn, CutData> p( ntk, rewriting_fn, ps, st );
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

namespace detail
{

template<typename Ntk, unsigned NInputs>
struct cut_match_t
{
  std::vector<exact_supergate<Ntk, NInputs>> const* supergates{NULL};
  std::array<uint8_t, NInputs> permutation;
  uint8_t negation{0};
};

template<typename Ntk, unsigned NInputs>
struct node_match_t
{
  exact_supergate<Ntk, NInputs> const* best_supergate[2] = {NULL, NULL};
  uint32_t best_cut[2] = {0, 0};

  float arrival[2];
  float required[2];
  float area[2];

  uint32_t map_refs[3];
  float flow_refs[3];
  float flows[3];
};

template<class NtkDest, class Ntk, class RewritingFn, typename CutData, unsigned NInputs>
class tech_map_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, true, CutData>;
  using cut_t = typename network_cuts_t::cut_t;
  using ref_t = std::array<uint32_t, 3>;
  using ref_est_t = std::array<float, 3>;

public:
  tech_map_impl( Ntk& ntk, exact_library<NtkDest, RewritingFn, NInputs> const& library, map_params const& ps, map_stats& st )
      : ntk( ntk ),
        library( library ),
        ps( ps ),
        st( st ),
        node_match( ntk.size() ),
        matches(),
        cuts( cut_enumeration<Ntk, true, CutData>( ntk, ps.cut_enumeration_ps ) )
  {
    map_update_cuts<CutData>().apply( cuts, ntk );
    std::tie( lib_inv_area, lib_inv_delay ) = library.get_inverter_info();
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

    compute_matches();
    init_nodes();

    /* compute mapping delay */
    if ( !ps.skip_delay_round )
    {
      compute_mapping<false>();
    }

    /* compute mapping global area flow */
    while ( iteration < ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      compute_mapping<true>();
    }

    /* compute mapping local exact area */
    while ( iteration < ps.ela_rounds + ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      compute_exact_area();
    }

    map_cover( res, old2new );

    st.area = area;
    st.delay = delay;

    return res;
  }

private:
  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n, auto ) {
      const auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        node_data.flow_refs[0] = node_data.flow_refs[1] = node_data.flow_refs[2] = 1.0f;
        node_data.arrival[0] = node_data.arrival[1] = 0.0f;
      }
      else
      {
        node_data.flow_refs[0] = node_data.flow_refs[1] = 0.0f;
        node_data.flow_refs[2] = static_cast<float>( ntk.fanout_size( n ) );
        ntk.foreach_fanin( n, [&]( auto const& s ) {
          if ( !ntk.is_pi( ntk.get_node( s ) ) )
          {
            const auto c_index = ntk.node_to_index( ntk.get_node( s ) );
            if ( ntk.is_complemented( s ) )
              node_match[c_index].flow_refs[1] += 1.0f;
            else
              node_match[c_index].flow_refs[0] += 1.0f;
          }
        } );
      }

      node_match[index].flows[2] = cuts.cuts( index )[0]->data.flow;
    } );
  }


  void compute_matches()
  {
    ntk.foreach_gate( [&]( auto const& n ) {
      const auto index = ntk.node_to_index( n );

      std::vector<cut_match_t<NtkDest, NInputs>> node_matches;

      auto i = 0u;
      for ( auto& cut : cuts.cuts( index ) )
      {
        const auto tt = cuts.truth_table( *cut );
        const auto fe = kitty::extend_to<NInputs>( tt );
        const auto config = kitty::exact_npn_canonization( fe );
        auto const supergates = library.get_supergates( std::get<0>( config ) );
        if ( supergates != NULL )
        {
          auto neg = std::get<1>( config );
          auto perm = std::get<2>( config );
          cut_match_t<NtkDest, NInputs> match;
          match.supergates = supergates;
          for ( auto j = 0u; j < perm.size() && j < NInputs; ++j )
          {
            match.permutation[perm[j]] = j;
          }
          match.negation = neg;
          node_matches.push_back( match );
          ( *cut )->data.match_index = i++;
        }
        else
        {
          std::cout << "Match failed" << std::endl;
          ( *cut )->data.ignore = true;
        }
      }
      
      matches[index] = node_matches;
    } );
  }

  template<bool DO_AREA>
  void compute_mapping()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      match_phase<DO_AREA>( n, 0u );
    }
    set_mapping_refs<false>();
  }


  void compute_exact_area()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      match_phase_exact( n, 0u );
    }
    set_mapping_refs<true>();
  }


  void map_cover( NtkDest& res, node_map<signal<NtkDest>, Ntk>& old2new )
  {
    auto const& db = library.get_database();

    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        return true;
      auto index = ntk.node_to_index( n );
      if ( node_match[index].map_refs[2] == 0u )
        return true;

      auto& best_cut = cuts.cuts( index )[node_match[index].best_cut[0]];

      std::vector<signal<NtkDest>> children( NInputs, res.get_constant( false ) );
      auto const& match = matches[index][best_cut->data.match_index];
      auto const& supergate = node_match[index].best_supergate[0];
      auto ctr = 0u;
      for ( auto l : best_cut )
      {
        children[match.permutation[ctr++]] = old2new[ntk.index_to_node( l )];
      }
      for ( auto i = 0u; i < NInputs; ++i )
      {
        if ( ( match.negation >> i ) & 1 )
        {
          children[match.permutation[i]] = !children[match.permutation[i]];
        }
      }
      topo_view topo{db, supergate->root};
      auto f = cleanup_dangling( topo, res, children.begin(), children.end() ).front();

      if ( ( match.negation >> NInputs ) & 1 )
        f = !f;
      
      old2new[n] = f;
      return true;
    } );

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? res.create_not( old2new[f] ) : old2new[f] );
    } );
  }

  template<bool ELA>
  void set_mapping_refs()
  {
    const auto coef = 1.0f / ( 2.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    if constexpr ( !ELA )
    {
      for ( auto i = 0u; i < node_match.size(); ++i )
      {
        node_match[i].map_refs[0] = node_match[i].map_refs[1] = node_match[i].map_refs[2] = 0u;
      }
    }

    /* compute current delay and update mapping refs */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, node_match[index].arrival[0] );
      if constexpr ( !ELA )
      {
        node_match[index].map_refs[2]++;
      }
    } );

    /* compute current area and update mapping refs */
    area = 0.0f;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      if ( node_match[index].map_refs[2] == 0u )
        continue;

      if constexpr ( !ELA )
      {
        for ( auto leaf : cuts.cuts( index )[node_match[index].best_cut[0]] )
        {
          node_match[leaf].map_refs[2]++;
        }
      }

      area += node_match[index].area[0];
    }

    /* blend flow references */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      node_match[i].flow_refs[2] = coef * node_match[i].flow_refs[2] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( node_match[i].map_refs[2] ) );
    }

    ++iteration;
  }

  void compute_required_time()
  {
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required[0] = node_match[i].required[1] = std::numeric_limits<float>::max();
    }
    
    /* return in case of first round of area optimization */
    if ( iteration == 0 )
      return;

    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      if ( ps.required_time == 0.0f )
      {
        node_match[index].required[0] = delay;
      }
      else
      {
        node_match[index].required[0] = ps.required_time;
      }
    } );

    auto i = ntk.size();
    while ( i-- > 0u )
    {
      const auto n = ntk.index_to_node( i );
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        break;

      if ( node_match[i].map_refs[2] == 0 )
        continue;

      auto ctr = 0u;
      auto best_cut = cuts.cuts( i )[node_match[i].best_cut[0]];
      auto const& match = matches[i][best_cut->data.match_index];
      auto const& supergate = node_match[i].best_supergate[0];
      for ( auto leaf : best_cut )
      {
        node_match[leaf].required[0] = std::min( node_match[leaf].required[0], node_match[i].required[0] + supergate->tdelay[match.permutation[ctr]] );
        ctr++;
      }
    }
  }

  template<bool DO_AREA>
  void match_phase( node<Ntk> const& n, unsigned phase )
  {
    float best_arrival = std::numeric_limits<float>::max();
    float best_area_flow = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];
    auto& cut_matches = matches[index];
    exact_supergate<NtkDest, NInputs> const* best_supergate = node_data.best_supergate[phase];

    /* recompute best match info */
    if ( best_supergate != NULL )
    {
      auto const& cut = cuts.cuts( index )[node_data.best_cut[phase]];
      auto& supergates = cut_matches[( cut )->data.match_index];
      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      best_arrival = 0.0f;
      best_area_flow = best_supergate->area + cut_leaves_flow( cut );
      best_area = best_supergate->area;
      best_cut = node_data.best_cut[phase];
      best_size = cut.size();
      for ( auto pin = 0u; pin < NInputs; pin++ )
      {
        float arrival_pin = node_match[children[pin]].arrival[phase] - best_supergate->tdelay[pin];
        best_arrival = std::max( best_arrival, arrival_pin );
      }
    }

    /* foreach cut */
    for ( auto& cut : cuts.cuts( index ) )
    {
      /* trivial cuts */
      if ( cut->size() == 1 || ( *cut )->data.ignore )
      {
        cut_index++;
        continue;
      }

      auto cut_flow = cut_leaves_flow( *cut );
      auto& supergates = cut_matches[( *cut )->data.match_index];

      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : *cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      for ( auto const& gate : *supergates.supergates )
      {
        // auto phase = gate.polarity ^ supergates.negation;
        float area_local = gate.area + cut_flow ;
        float worst_arrival = 0.0f;
        for ( auto pin = 0u; pin < NInputs; pin++ )
        {
          float arrival_pin = node_match[children[pin]].arrival[phase] - gate.tdelay[pin];
          worst_arrival = std::max( worst_arrival, arrival_pin );
        }

        if constexpr ( DO_AREA )
        {
          if ( worst_arrival > node_data.required[phase] + epsilon )
            continue;
        }

        if ( compare_map<DO_AREA>( worst_arrival, best_arrival, area_local, best_area_flow, cut->size(), best_size ) )
        {
          best_arrival = worst_arrival;
          best_area_flow = area_local;
          best_size = cut->size();
          best_cut = cut_index;
          best_area = gate.area;
          best_supergate = &gate;
        }
      }

      cut_index++;
    }

    node_data.flows[2] = best_area_flow / node_data.flow_refs[2];
    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.best_supergate[phase] = best_supergate;
  }

  void match_phase_exact( node<Ntk> const& n, unsigned phase )
  {
    float best_arrival = std::numeric_limits<float>::max();
    float best_exact_area = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];
    auto& cut_matches = matches[index];
    exact_supergate<NtkDest, NInputs> const* best_supergate = node_data.best_supergate[phase];


    /* recompute best match info */
    if ( best_supergate != NULL )
    {
      best_cut = node_data.best_cut[phase];
      auto const& cut = cuts.cuts( index )[best_cut];
      auto& supergates = cut_matches[( cut )->data.match_index];
      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      best_arrival = 0.0f;
      best_area = best_supergate->area;
      best_size = cut.size();
      for ( auto pin = 0u; pin < NInputs; pin++ )
      {
        float arrival_pin = node_match[children[pin]].arrival[phase] - best_supergate->tdelay[pin];
        best_arrival = std::max( best_arrival, arrival_pin );
      }
    }

    if ( node_data.map_refs[2] )
    {
      best_exact_area = cut_deref( cuts.cuts( index )[best_cut], n, phase );
    }

    /* foreach cut */
    for ( auto& cut : cuts.cuts( index ) )
    {
      /* trivial cuts */
      if ( cut->size() == 1 || ( *cut )->data.ignore )
      {
        cut_index++;
        continue;
      }

      auto& supergates = cut_matches[( *cut )->data.match_index];

      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : *cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      auto area_exact = cut_ref( *cut, n, phase ) - node_data.area[phase];
      cut_deref( *cut, n, phase );

      for ( auto const& gate : *supergates.supergates )
      {
        // auto phase = gate.polarity ^ supergates.negation;
        auto area_exact_gate = area_exact + gate.area;
        float worst_arrival = 0.0f;
        for ( auto pin = 0u; pin < NInputs; pin++ )
        {
          float arrival_pin = node_match[children[pin]].arrival[phase] - gate.tdelay[pin];
          worst_arrival = std::max( worst_arrival, arrival_pin );
        }

        if ( worst_arrival > node_data.required[phase] + epsilon )
          continue;

        if ( compare_map<true>( worst_arrival, best_arrival, area_exact_gate, best_exact_area, cut->size(), best_size ) )
        {
          best_arrival = worst_arrival;
          best_exact_area = area_exact_gate;
          best_area = gate.area;
          best_size = cut->size();
          best_cut = cut_index;
          best_supergate = &gate;
        }
      }

      cut_index++;
    }

    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.best_supergate[phase] = best_supergate;

    if ( node_data.map_refs[2] )
    {
      cut_ref( cuts.cuts( index )[best_cut], n, phase );
    }
  }

  inline float cut_leaves_flow( cut_t const& cut )
  {
    float flow{0.0f};

    for ( auto leaf : cut )
    {
      flow += node_match[leaf].flows[2];
    }

    return flow;
  }

  float cut_ref( cut_t const& cut, node<Ntk> const& n, unsigned phase )
  {
    float count = node_match[ntk.node_to_index( n )].area[phase];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( node_match[leaf].map_refs[2]++ == 0u )
      {
        count += cut_ref( cuts.cuts( leaf )[node_match[leaf].best_cut[phase]], ntk.index_to_node( leaf ), phase );
      }
    }
    return count;
  }

  float cut_deref( cut_t const& cut, node<Ntk> const& n, unsigned phase )
  {
    float count = node_match[ntk.node_to_index( n )].area[phase];
    for ( auto leaf : cut )
    {
      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) || ntk.is_pi( ntk.index_to_node( leaf ) ) )
        continue;

      if ( --node_match[leaf].map_refs[2] == 0u )
      {
        count += cut_deref( cuts.cuts( leaf )[node_match[leaf].best_cut[phase]], ntk.index_to_node( leaf ), phase );
      }
    }
    return count;
  }

  template<bool DO_AREA>
  inline bool compare_map( float arrival, float best_arrival, float area_flow, float best_area_flow, uint32_t size, uint32_t best_size )
  {
    if constexpr ( DO_AREA )
    {
      if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
      else if ( arrival < best_arrival - epsilon )
      {
        return true;
      }
      else if ( arrival > best_arrival + epsilon )
      {
        return false;
      }
    }
    else
    {
      if ( arrival < best_arrival - epsilon )
      {
        return true;
      }
      else if ( arrival > best_arrival + epsilon )
      {
        return false;
      }
      else if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
    }
    if ( size < best_size )
    {
      return true;
    }
    /* TODO: add compare on fanout size */
    return false;
  }


private:
  Ntk& ntk;
  exact_library<NtkDest, RewritingFn, NInputs> const& library;
  map_params const& ps;
  map_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  float delay{0.0f};     /* current delay of the mapping */
  float area{0.0f};      /* current area of the mapping */
  const float epsilon{0.005f}; /* epsilon */

  /* lib inverter info */
  float lib_inv_area;
  float lib_inv_delay;

  std::vector<node<Ntk>> top_order;
  std::vector<node_match_t<NtkDest, NInputs>> node_match;
  std::unordered_map<uint32_t, std::vector<cut_match_t<NtkDest, NInputs>>> matches;
  network_cuts_t cuts;
};

} /* namespace detail */

template<class Ntk, class NtkDest = Ntk, class RewritingFn, unsigned NInputs, typename CutData = cut_enumeration_map_cut>
NtkDest tech_map( Ntk& ntk, exact_library<NtkDest, RewritingFn, NInputs> const& library, map_params const& ps = {}, map_stats* pst = nullptr )
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

  map_stats st;
  detail::tech_map_impl<NtkDest, Ntk, RewritingFn, CutData, NInputs> p( ntk, library, ps, st );
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
