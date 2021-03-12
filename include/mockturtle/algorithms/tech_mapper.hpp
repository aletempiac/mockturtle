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

#include "mapper.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/tech_library.hpp"
#include "../views/topo_view.hpp"
#include "cut_enumeration/tech_map_cut.hpp"

namespace mockturtle
{

/*! \brief Parameters for lut_mapping.
 *
 * The data structure `lut_mapping_params` holds configurable parameters
 * with default arguments for `lut_mapping`.
 */
// struct map_params
// {
//   map_params()
//   {
//     cut_enumeration_ps.cut_size = 4;
//     cut_enumeration_ps.cut_limit = 8;
//     cut_enumeration_ps.minimize_truth_table = true;
//   }

//   /*! \brief Parameters for cut enumeration
//    *
//    * The default cut size is 4, the default cut limit is 8.
//    */
//   cut_enumeration_params cut_enumeration_ps{};

//   /*! \brief Required time for delay optimization. */
//   float required_time{0.0f};

//   /*! \brief Do area optimization. */
//   bool skip_delay_round{false};

//   /*! \brief Number of rounds for area flow optimization. */
//   uint32_t area_flow_rounds{1u};

//   /*! \brief Number of rounds for exact area optimization. */
//   uint32_t ela_rounds{1u};

//   /*! \brief Use structural choices. */
//   bool choices{false};

//   /*! \brief Be verbose. */
//   bool verbose{false};
// };

/*! \brief Statistics for mapper.
 *
 * The data structure `mapper_stats` provides data collected by running
 * `mapper`.
 */
// struct map_stats
// {
//   /* \brief Area and delay */
//   double area{0};
//   double delay{0};
//   /*! \brief Total runtime. */
//   stopwatch<>::duration time_total{0};

//   void report() const
//   {
//     std::cout << fmt::format( "[i] area = {:>5.2f}; delay = {:>5.2f}\n", area, delay );
//     std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
//   }
// };

/* function to update all cuts after cut enumeration */
// template<typename CutData>
// struct map_update_cuts
// {
//   template<typename NetworkCuts, typename Ntk>
//   static void apply( NetworkCuts const& cuts, Ntk const& ntk )
//   {
//     (void)cuts;
//     (void)ntk;
//   }
// };

namespace detail
{

template<unsigned NInputs>
struct cut_match_tech
{
  std::vector<supergate<NInputs>> const* supergates[2] = {NULL, NULL};
  std::array<uint8_t, NInputs> permutation;
  uint8_t negation{0};
};

template<unsigned NInputs>
struct node_match_tech
{
  supergate<NInputs> const* best_supergate[2] = {NULL, NULL};
  uint8_t phase[2];
  uint32_t best_cut[2];
  bool same_match{false};

  float arrival[2];
  float required[2];
  float area[2];

  uint32_t map_refs[3];
  float flow_refs[3];
  float flows[3];
};

template<class Ntk, unsigned NInputs, typename CutData>
class tech_mapping_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, true, CutData>;
  using cut_t = typename network_cuts_t::cut_t;
  using ref_t = std::array<uint32_t, 3>;
  using ref_est_t = std::array<float, 3>;

public:
  tech_mapping_impl( Ntk& ntk, tech_library<NInputs> const& library, map_params const& ps, map_stats& st )
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

  klut_network run()
  {
    stopwatch t( st.time_total );

    auto [res, old2new] = initialize_copy_network<klut_network>( ntk );
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
      if ( !compute_mapping<false>() )
      {
        return res;
      }
    }

    /* compute mapping global area flow */
    while ( iteration < ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      if ( !compute_mapping<true>() )
      {
        return res;
      }
    }

    /* compute mapping local exact area */
    while ( iteration < ps.ela_rounds + ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      if ( !compute_exact_area() )
      {
        return res;
      }
    }

    /* write final results */
    st.area = area;
    st.delay = delay;

    finalize_cover( res, old2new );

    return res;
  }

private:
  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n, auto ) {
      const auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      if ( ntk.is_constant( n ) )
      {
        /* all terminals have flow 1.0 */
        node_data.flow_refs[0] = node_data.flow_refs[1] = node_data.flow_refs[2] = 1.0f;
        node_data.arrival[0] = node_data.arrival[1] = 0.0f;
      }
      else if ( ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        node_data.flow_refs[0] = node_data.flow_refs[1] = node_data.flow_refs[2] = 1.0f;
        node_data.arrival[0] = 0.0f;
        node_data.arrival[1] = lib_inv_delay;
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

      std::vector<cut_match_tech<NInputs>> node_matches;

      auto i = 0u;
      for ( auto& cut : cuts.cuts( index ) )
      {
        if ( cut->size() == 1 )
        {
          ( *cut )->data.ignore = true;
          continue;
        }
        const auto tt = cuts.truth_table( *cut );
        if ( tt.num_vars() > NInputs )
        {
          /* Ignore cuts too big to be mapped */
          ( *cut )->data.ignore = true;
          continue;
        }
        const auto fe = kitty::extend_to<NInputs>( tt );
        const auto config = kitty::exact_npn_canonization( fe );
        auto const supergates_npn = library.get_supergates( std::get<0>( config ) );
        auto const supergates_npn_neg = library.get_supergates( ~std::get<0>( config ) );
        if ( supergates_npn != NULL || supergates_npn_neg != NULL )
        {
          auto neg = std::get<1>( config );
          auto perm = std::get<2>( config );
          uint8_t phase = ( neg >> NInputs ) & 1;
          cut_match_tech<NInputs> match;
          if ( supergates_npn != NULL )
            match.supergates[phase] = supergates_npn;
          if ( supergates_npn_neg != NULL )
            match.supergates[phase ^ 1] = supergates_npn_neg;
          match.negation = 0;
          for ( auto j = 0u; j < perm.size() && j < NInputs; ++j )
          {
            match.permutation[perm[j]] = j;
            match.negation |= ( ( neg >> perm[j] ) & 1 ) << j;
          }
          node_matches.push_back( match );
          ( *cut )->data.match_index = i++;
        }
        else
        {
          ( *cut )->data.ignore = true;
        }
      }
      
      matches[index] = node_matches;
    } );
  }

  template<bool DO_AREA>
  bool compute_mapping()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      /* match positive phase */
      match_phase<DO_AREA>( n, 0u );

      /* match negative phase */
      match_phase<DO_AREA>( n, 1u );

      /* try to drop one phase */
      match_drop_phase<DO_AREA, false>( n, 0 );
    }
    return set_mapping_refs<false>();
  }


  bool compute_exact_area()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      /* recursive deselect the best cut in common if used */
      if ( node_data.same_match && node_data.map_refs[2] != 0 )
      {
        if ( node_data.best_supergate[0] != NULL )
          cut_deref( cuts.cuts( index )[node_data.best_cut[0]], n, 0u );
        else
          cut_deref( cuts.cuts( index )[node_data.best_cut[1]], n, 1u );
      }

      /* match positive phase */
      match_phase_exact( n, 0u );

      /* match negative phase */
      match_phase_exact( n, 1u );

      /* try to drop one phase */
      match_drop_phase<true, true>( n, 0 );
    }
    return set_mapping_refs<true>();
  }

  template<bool ELA>
  bool set_mapping_refs()
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
      if ( ntk.is_complemented( s ) )
        delay = std::max( delay, node_match[index].arrival[1] );
      else
        delay = std::max( delay, node_match[index].arrival[0] );
      if constexpr ( !ELA )
      {
        node_match[index].map_refs[2]++;
        if ( ntk.is_complemented( s ) )
          node_match[index].map_refs[1]++;
        else
          node_match[index].map_refs[0]++;
      }
    } );

    /* compute current area and update mapping refs */
    area = 0.0f;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      const auto index = ntk.node_to_index( *it );
      /* skip constants and PIs */
      if ( ntk.is_pi( *it ) )
      {
        if ( node_match[index].map_refs[1] > 0u )
        {
          /* Add inverter over the negated fanins */
          area += lib_inv_area;
        }
        continue;
      }
      else if ( ntk.is_constant( *it ) )
      {
        /* TODO add map on constants */
        continue;
      }

      if ( node_match[index].map_refs[2] == 0u )
        continue;

      auto& node_data = node_match[index];

      unsigned use_phase = node_data.best_supergate[0] == NULL ? 1u : 0u;

      if ( node_data.best_supergate[use_phase] == NULL )
      {
        /* Library is not complete, mapping is not possible */
        std::cerr << "MAP ERROR: technology library is not complete, impossible to perform mapping" << std::endl;
        return false;
      }

      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts.cuts( index )[node_data.best_cut[use_phase]];
          auto const& match = matches[index][best_cut->data.match_index];
          auto ctr = 0u;
          for ( auto const leaf : best_cut )
          {
            node_match[leaf].map_refs[2]++;
            if ( ( node_data.phase[use_phase] >> match.permutation[ctr++] ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase];
        if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        {
          area += lib_inv_area;
        }
      }
      use_phase = use_phase ^ 1;
        /* if both phases are implemented and used */
      if ( !node_data.same_match && node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts.cuts( index )[node_data.best_cut[use_phase]];
          auto const& match = matches[index][best_cut->data.match_index];
          auto ctr = 0u;
          for ( auto const leaf : best_cut )
          {
            node_match[leaf].map_refs[2]++;
            if ( ( node_data.phase[use_phase] >> match.permutation[ctr++] ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase];
      }
    }

    /* blend flow references */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      node_match[i].flow_refs[2] = coef * node_match[i].flow_refs[2] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( node_match[i].map_refs[2] ) );
      node_match[i].flow_refs[1] = coef * node_match[i].flow_refs[1] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( node_match[i].map_refs[1] ) );
      node_match[i].flow_refs[0] = coef * node_match[i].flow_refs[0] + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( node_match[i].map_refs[0] ) );
    }

    ++iteration;
    return true;
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

    auto required = delay;
    if ( ps.required_time != 0.0f )
    {
      /* Global target time constraint */
      if ( ps.required_time < delay - epsilon )
      {
        if ( !ps.skip_delay_round && iteration == 1 )
          std::cerr << fmt::format( "MAP WARNING: cannot meet the target required time of {:.2f}", ps.required_time ) << std::endl;
      }
      else
      {
        required = ps.required_time;
      }
    }

    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      if ( ntk.is_complemented( s ) )
        node_match[index].required[1] = required;
      else
        node_match[index].required[0] = required;
    } );

    auto i = ntk.size();
    while ( i-- > 0u )
    {
      const auto n = ntk.index_to_node( i );
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        break;

      if ( node_match[i].map_refs[2] == 0 )
        continue;

      auto& node_data = node_match[i];

      unsigned use_phase = node_data.best_supergate[0] == NULL ? 1u : 0u;
      unsigned other_phase = use_phase ^ 1;

      assert( node_data.best_supergate[0] != NULL || node_data.best_supergate[1] != NULL );
      assert( node_data.map_refs[0] || node_data.map_refs[1] );

      /* propagate required time over output inverter if present */
      if ( node_data.same_match && node_data.map_refs[other_phase] > 0 )
      {
        node_data.required[use_phase] = std::min( node_data.required[use_phase], node_data.required[other_phase] - lib_inv_delay );   
      }

      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
      {
        auto ctr = 0u;
        auto best_cut = cuts.cuts( i )[node_data.best_cut[use_phase]];
        auto const& match = matches[i][best_cut->data.match_index];
        auto const& supergate = node_data.best_supergate[use_phase];
        for ( auto leaf : best_cut )
        {
          auto phase = ( node_data.phase[use_phase] >> match.permutation[ctr] ) & 1;
          node_match[leaf].required[phase] = std::min( node_match[leaf].required[phase], node_data.required[use_phase] - supergate->tdelay[match.permutation[ctr]] );
          ctr++;
        }
      }

      if ( !node_data.same_match && node_data.map_refs[other_phase] > 0 )
      {
        auto ctr = 0u;
        auto best_cut = cuts.cuts( i )[node_data.best_cut[other_phase]];
        auto const& match = matches[i][best_cut->data.match_index];
        auto const& supergate = node_data.best_supergate[other_phase];
        for ( auto leaf : best_cut )
        {
          auto phase = ( node_data.phase[other_phase] >> match.permutation[ctr] ) & 1;
          node_match[leaf].required[phase] = std::min( node_match[leaf].required[phase], node_data.required[other_phase] - supergate->tdelay[match.permutation[ctr]] );
          ctr++;
        }
      }
    }
  }

  template<bool DO_AREA>
  void match_phase( node<Ntk> const& n, uint8_t phase )
  {
    float best_arrival = std::numeric_limits<float>::max();
    float best_area_flow = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t best_phase = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];
    auto& cut_matches = matches[index];
    supergate<NInputs> const* best_supergate = node_data.best_supergate[phase];

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

      best_phase = node_data.phase[phase];
      best_arrival = 0.0f;
      best_area_flow = best_supergate->area + cut_leaves_flow( cut, n, phase );
      best_area = best_supergate->area;
      best_cut = node_data.best_cut[phase];
      best_size = cut.size();
      for ( auto pin = 0u; pin < NInputs; pin++ )
      {
        float arrival_pin = node_match[children[pin]].arrival[( best_phase >> pin ) & 1] + best_supergate->tdelay[pin];
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

      auto const& supergates = cut_matches[( *cut )->data.match_index];

      if ( supergates.supergates[phase] == NULL )
      {
        cut_index++;
        continue;
      }

      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : *cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      for ( auto const& gate : *supergates.supergates[phase] )
      {
        uint8_t complement = gate.polarity ^ supergates.negation;
        node_data.phase[phase] = complement;
        float area_local = gate.area + cut_leaves_flow( *cut, n, phase );
        float worst_arrival = 0.0f;
        for ( auto pin = 0u; pin < NInputs; pin++ )
        {
          float arrival_pin = node_match[children[pin]].arrival[( complement >> pin ) & 1] + gate.tdelay[pin];
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
          best_phase = complement;
          best_supergate = &gate;
        }
      }

      cut_index++;
    }

    node_data.flows[phase] = best_area_flow;
    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.phase[phase] = best_phase;
    node_data.best_supergate[phase] = best_supergate;
  }

  void match_phase_exact( node<Ntk> const& n, uint8_t phase )
  {
    float best_arrival = std::numeric_limits<float>::max();
    float best_exact_area = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t best_phase = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];
    auto& cut_matches = matches[index];
    supergate<NInputs> const* best_supergate = node_data.best_supergate[phase];

    /* recompute best match info */
    if ( best_supergate != NULL )
    {
      auto const& cut = cuts.cuts( index )[node_data.best_cut[phase]];
      auto const& supergates = cut_matches[( cut )->data.match_index];
      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      best_phase = best_supergate->polarity ^ supergates.negation;
      best_arrival = 0.0f;
      best_area = best_supergate->area;
      best_cut = node_data.best_cut[phase];
      best_size = cut.size();
      for ( auto pin = 0u; pin < NInputs; pin++ )
      {
        float arrival_pin = node_match[children[pin]].arrival[( best_phase >> pin ) & 1] + best_supergate->tdelay[pin];
        best_arrival = std::max( best_arrival, arrival_pin );
      }

      /* if cut is implemented, remove it from the cover */
      if ( !node_data.same_match && node_data.map_refs[phase] )
      {
        best_exact_area = cut_deref( cuts.cuts( index )[best_cut], n, phase );
      }
      else
      {
        best_exact_area = cut_ref( cuts.cuts( index )[best_cut], n, phase );
        cut_deref( cuts.cuts( index )[best_cut], n, phase );
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

      auto const& supergates = cut_matches[( *cut )->data.match_index];

      if ( supergates.supergates[phase] == NULL )
      {
        cut_index++;
        continue;
      }

      std::vector<uint32_t> children( NInputs, 0u );
      auto ctr = 0u;
      for ( auto l : *cut )
      {
        children[supergates.permutation[ctr++]] = l;
      }

      for ( auto const& gate : *supergates.supergates[phase] )
      {
        uint8_t complement = gate.polarity ^ supergates.negation;
        node_data.phase[phase] = complement;
        node_data.area[phase] = gate.area;
        auto area_exact = cut_ref( *cut, n, phase );
        cut_deref( *cut, n, phase );
        float worst_arrival = 0.0f;
        for ( auto pin = 0u; pin < NInputs; pin++ )
        {
          float arrival_pin = node_match[children[pin]].arrival[( complement >> pin ) & 1] + gate.tdelay[pin];
          worst_arrival = std::max( worst_arrival, arrival_pin );
        }

        if ( worst_arrival > node_data.required[phase] + epsilon )
          continue;

        if ( compare_map<true>( worst_arrival, best_arrival, area_exact, best_exact_area, cut->size(), best_size ) )
        {
          best_arrival = worst_arrival;
          best_exact_area = area_exact;
          best_area = gate.area;
          best_size = cut->size();
          best_cut = cut_index;
          best_phase = complement;
          best_supergate = &gate;
        }
      }

      cut_index++;
    }

    node_data.flows[phase] = best_exact_area;
    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.phase[phase] = best_phase;
    node_data.best_supergate[phase] = best_supergate;

    if ( !node_data.same_match && node_data.map_refs[phase] )
    {
      best_exact_area = cut_ref( cuts.cuts( index )[best_cut], n, phase );
    }
  }

  template<bool DO_AREA, bool ELA>
  void match_drop_phase( node<Ntk> const& n, float required_margin_factor )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    /* compute arrival adding an inverter to the other match phase */
    float worst_arrival_npos = node_data.arrival[1] + lib_inv_delay;
    float worst_arrival_nneg = node_data.arrival[0] + lib_inv_delay;
    bool use_zero = false;
    bool use_one = false;

    /* only one phase is matched */
    if ( node_data.best_supergate[0] == NULL )
    {
      set_match_complemented_phase( index, 1, worst_arrival_npos );
      if constexpr ( ELA )
      {
        if ( node_data.map_refs[2] )
          cut_ref( cuts.cuts( index )[node_data.best_cut[1]], n, 1 );
      }
      return;
    }
    else if ( node_data.best_supergate[1] == NULL )
    {
      set_match_complemented_phase( index, 0, worst_arrival_nneg );
      if constexpr ( ELA )
      {
        if ( node_data.map_refs[2] )
          cut_ref( cuts.cuts( index )[node_data.best_cut[0]], n, 0 );
      }
      return;
    }

    /* try to use only one match to cover both phases */
    if constexpr ( !DO_AREA )
    {
      /* if arrival is less matching the other phase and inserting an inverter */
      if ( worst_arrival_npos < node_data.arrival[0] + epsilon )
      {
        use_one = true;
      }
      if ( worst_arrival_nneg < node_data.arrival[1] + epsilon )
      {
        use_zero = true;
      }
    }
    else
    {
      /* check if both phases + inverter meet the required time */
      use_zero = worst_arrival_nneg < ( node_data.required[1] + epsilon - required_margin_factor * lib_inv_delay );
      use_one = worst_arrival_npos < ( node_data.required[0] + epsilon - required_margin_factor * lib_inv_delay );
    }


    if ( !use_zero && !use_one )
    {
      /* use both phases to improve delay */
      node_data.flows[2] = ( node_data.flows[0] + node_data.flows[1] ) / node_data.flow_refs[2];
      node_data.flows[0] = node_data.flows[0] / node_data.flow_refs[0];
      node_data.flows[1] = node_data.flows[1] / node_data.flow_refs[1];
      node_data.same_match = false;
      return;
    }

    if constexpr ( DO_AREA )
    {
      /* condition on not used phases, evaluate a substitution */
      if ( iteration != 0 )
      {
        if ( node_data.map_refs[0] == 0 || node_data.map_refs[1] == 0 )
        {
          /* select the used match */
          auto phase = 0;
          auto nphase = 0;
          if ( node_data.map_refs[0] == 0 )
          {
            phase = 1;
            use_one = true;
            use_zero = false;
          }
          else
          {
            nphase = 1;
            use_one = false;
            use_zero = true;
          }
          /* select the not used match instead if it leads to improvement and doesn't violate the required time */
          if ( node_data.arrival[nphase] + lib_inv_delay < node_data.required[phase] + epsilon )
          {
            auto size_phase = cuts.cuts( index )[node_data.best_cut[phase]].size();
            auto size_nphase = cuts.cuts( index )[node_data.best_cut[nphase]].size();
            auto inverter_cost = 0;
            if ( ELA )
              inverter_cost = lib_inv_area;
            if ( compare_map<DO_AREA>( node_data.arrival[nphase] + lib_inv_delay, node_data.arrival[phase],  node_data.flows[nphase] + inverter_cost, node_data.flows[phase], size_nphase, size_phase ) )
            {
              /* invert the choice */
              use_zero = !use_zero;
              use_one = !use_one;
            }
          }
        }
      }
    }

    /* use area flow as a tiebreaker */
    if ( use_zero && use_one )
    {
      auto size_zero = cuts.cuts( index )[node_data.best_cut[0]].size();
      auto size_one = cuts.cuts( index )[node_data.best_cut[1]].size();
      if ( compare_map<DO_AREA>( worst_arrival_nneg, worst_arrival_npos,  node_data.flows[0], node_data.flows[1], size_zero, size_one ) )
        use_one = false;
      else
        use_zero = false;
    }

    if ( use_zero )
    {
      if constexpr ( ELA )
      {
        if ( !node_data.same_match ) 
        {
          if ( node_data.map_refs[1] > 0 )
            cut_deref( cuts.cuts( index )[node_data.best_cut[1]], n, 1 );
          if ( node_data.map_refs[0] == 0 && node_data.map_refs[2] )
            cut_ref( cuts.cuts( index )[node_data.best_cut[0]], n, 0 );
        }
        else if ( node_data.map_refs[2] )
          cut_ref( cuts.cuts( index )[node_data.best_cut[0]], n, 0 );
      }
      set_match_complemented_phase( index, 0, worst_arrival_nneg );
    }
    else
    {
      if constexpr ( ELA )
      {
        if ( !node_data.same_match ) 
        {
          if ( node_data.map_refs[0] > 0 )
            cut_deref( cuts.cuts( index )[node_data.best_cut[0]], n, 0 );
          if ( node_data.map_refs[1] == 0 && node_data.map_refs[2] )
            cut_ref( cuts.cuts( index )[node_data.best_cut[1]], n, 1 );
        }
        else if ( node_data.map_refs[2] )
          cut_ref( cuts.cuts( index )[node_data.best_cut[1]], n, 1 );
      }
      set_match_complemented_phase( index, 1, worst_arrival_npos );
    }
  }

  inline void set_match_complemented_phase( uint32_t index, uint8_t phase, float worst_arrival_n )
  {
    auto& node_data = node_match[index];
    auto phase_n = phase ^ 1;
    node_data.same_match = true;
    node_data.best_supergate[phase_n] = NULL;
    node_data.best_cut[phase_n] = node_data.best_cut[phase];
    node_data.phase[phase_n] = node_data.phase[phase] ^ ( 1 << NInputs );
    node_data.arrival[phase_n] = worst_arrival_n;
    node_data.area[phase_n] = node_data.area[phase];
    node_data.flows[phase_n] = node_data.flows[phase] / node_data.flow_refs[2];
    node_data.flows[2] = node_data.flows[phase] / node_data.flow_refs[2];
    node_data.flows[phase] = node_data.flows[phase] / node_data.flow_refs[2];
  }

  inline float cut_leaves_flow( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    float flow{0.0f};
    auto const& node_data = node_match[ntk.node_to_index( n )];
    auto const& match = matches[ntk.node_to_index( n )][cut->data.match_index];

    uint8_t ctr = 0u;
    for ( auto leaf : cut )
    {
      uint8_t leaf_phase = ( node_data.phase[phase] >> match.permutation[ctr++] ) & 1;
      flow += node_match[leaf].flows[leaf_phase];
    }

    return flow;
  }

  float cut_ref( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    auto const& node_data = node_match[ntk.node_to_index( n )];
    auto const& match = matches[ntk.node_to_index( n )][cut->data.match_index];
    float count = node_data.area[phase];
    uint8_t ctr = 0;
    for ( auto leaf : cut )
    {
      /* compute leaf phase using the current gate */
      uint8_t leaf_phase = ( node_data.phase[phase] >> match.permutation[ctr] ) & 1;

      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        ++ctr;
        continue;
      }
      else if ( ntk.is_pi( ntk.index_to_node( leaf ) ) )
      {
        /* inverter cost for PIs */
        if ( leaf_phase == 1u )
        {
          if ( node_match[leaf].map_refs[1]++ == 0u )
            count += lib_inv_area;
        }
        ++ctr;
        continue;
      }

      if ( node_match[leaf].same_match )
      {
        /* Add inverter area if not present yet and leaf node is implemented in the opposite phase */
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u && node_match[leaf].best_supergate[leaf_phase] == NULL )
          count += lib_inv_area;
        /* Recursive referencing if leaf was not referenced */
        if ( node_match[leaf].map_refs[2]++ == 0u )
        {
          count += cut_ref( cuts.cuts( leaf )[node_match[leaf].best_cut[leaf_phase]], ntk.index_to_node( leaf ), leaf_phase );
        }
      }
      else
      {
        ++node_match[leaf].map_refs[2];
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u )
        {
          count += cut_ref( cuts.cuts( leaf )[node_match[leaf].best_cut[leaf_phase]], ntk.index_to_node( leaf ), leaf_phase );
        }
      }
      ++ctr;
    }
    return count;
  }

  float cut_deref( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    auto const& node_data = node_match[ntk.node_to_index( n )];
    auto const& match = matches[ntk.node_to_index( n )][cut->data.match_index];
    float count = node_data.area[phase];
    uint8_t ctr = 0;
    for ( auto leaf : cut )
    {
      /* compute leaf phase using the current gate */
      uint8_t leaf_phase = ( node_data.phase[phase] >> match.permutation[ctr] ) & 1;

      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        ++ctr;
        continue;
      }
      else if ( ntk.is_pi( ntk.index_to_node( leaf ) ) )
      {
        /* inverter cost for PIs */
        if ( leaf_phase == 1u )
        {
          if ( --node_match[leaf].map_refs[1] == 0u )
            count += lib_inv_area;
        }
        ++ctr;
        continue;
      }

      if ( node_match[leaf].same_match )
      {
        /* Add inverter area if it is used only by the current gate and leaf node is implemented in the opposite phase */
        if ( --node_match[leaf].map_refs[leaf_phase] == 0u && node_match[leaf].best_supergate[leaf_phase] == NULL )
          count += lib_inv_area;
        /* Recursive dereferencing */
        if ( --node_match[leaf].map_refs[2] == 0u )
        {
          count += cut_deref( cuts.cuts( leaf )[node_match[leaf].best_cut[leaf_phase]], ntk.index_to_node( leaf ), leaf_phase );
        }
      }
      else
      {
        --node_match[leaf].map_refs[2];
        if ( --node_match[leaf].map_refs[leaf_phase] == 0u )
        {
          count += cut_deref( cuts.cuts( leaf )[node_match[leaf].best_cut[leaf_phase]], ntk.index_to_node( leaf ), leaf_phase );
        }
      }
      ++ctr;
    }
    return count;
  }

  void finalize_cover( klut_network& res, node_map<signal<klut_network>, Ntk>& old2new )
  {
    ntk.foreach_node( [&]( auto const& n ) {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        return true;
      auto index = ntk.node_to_index( n );
      if ( node_match[index].map_refs[2] == 0u )
        return true;

      unsigned phase = ( node_match[index].best_supergate[0] != NULL ) ? 0 : 1;
      auto& best_cut = cuts.cuts( index )[node_match[index].best_cut[phase]];

      std::vector<signal<klut_network>> children( best_cut.size() );

      auto ctr = 0u;
      for ( auto l : best_cut )
      {
        children[ctr++] = old2new[ntk.index_to_node( l )];
      }
      auto f = res.create_node( children, cuts.truth_table( best_cut ) );

      old2new[n] = f;

      /* check correctness */
      // auto& node_data = node_match[index];
      // auto& cut_match = matches[index][best_cut->data.match_index];
      // auto tt = kitty::extend_to<NInputs>( cuts.truth_table( best_cut ) );
      // if ( ( phase ) == 1 )
      //   tt = ~tt;
      // std::array<uint8_t, NInputs> perm;
      // for ( auto i = 0u; i < NInputs; ++i )
      // {
      //   perm[cut_match.permutation[i]] = i;
      // }
      // for ( auto i = 0u; i < NInputs - 1; ++i )
      // {
      //   /* swap */
      //   if ( i != perm[i] )
      //   {
      //     kitty::swap_inplace( tt, i, perm[i] );
      //     for ( auto j = i + 1u; j < NInputs; ++j )
      //     {
      //       if ( perm[j] == i )
      //       {
      //         perm[j] = perm[i];
      //         break;
      //       }
      //     }
      //     perm[i] = i;
      //   }
      //   /* flip */
      //   if ( ( ( node_data.phase[phase] >> i ) & 1 ) == 1 )
      //     kitty::flip_inplace( tt, i );
      // }
      // if ( ( ( node_data.phase[phase] >> ( NInputs - 1 ) ) & 1 ) == 1 )
      //     kitty::flip_inplace( tt, NInputs - 1 );

      // /* obtained NP canonical truth table */
      // /* apply supergate negation on top and permutate again */
      // auto const supergate = node_data.best_supergate[phase];
      // for ( auto i = 0u; i < NInputs; ++i )
      // {
      //   perm[i] = supergate->permutation[i];
      // }
      // for ( auto i = 0u; i < NInputs - 1; ++i )
      // {
      //   if ( i != perm[i] )
      //   {
      //     kitty::swap_inplace( tt, i, perm[i] );
      //     for ( auto j = i + 1u; j < NInputs; ++j )
      //     {
      //       if ( perm[j] == i )
      //       {
      //         perm[j] = perm[i];
      //         break;
      //       }
      //     }
      //     perm[i] = i;
      //   }
      // }
      // /* check truth table is equal */
      // auto gate_tt = kitty::extend_to<NInputs>( supergate->root->function );
      // assert( gate_tt == tt );
      return true;
    } );

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      res.create_po( ntk.is_complemented( f ) ? res.create_not( old2new[f] ) : old2new[f] );
    } );
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
  tech_library<NInputs> const& library;
  map_params const& ps;
  map_stats& st;

  uint32_t iteration{0}; /* current mapping iteration */
  float delay{0.0f};     /* current delay of the mapping */
  double area{0.0f};      /* current area of the mapping */
  const float epsilon{0.005f}; /* epsilon */

  /* lib inverter info */
  float lib_inv_area;
  float lib_inv_delay;

  std::vector<node<Ntk>> top_order;
  std::vector<node_match_tech<NInputs>> node_match;
  std::unordered_map<uint32_t, std::vector<cut_match_tech<NInputs>>> matches;
  network_cuts_t cuts;
};

} /* namespace detail */

template<class Ntk, unsigned NInputs, typename CutData = cut_enumeration_tech_map_cut>
klut_network tech_mapping( Ntk& ntk, tech_library<NInputs> const& library, map_params const& ps = {}, map_stats* pst = nullptr )
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
  detail::tech_mapping_impl<Ntk, NInputs, CutData> p( ntk, library, ps, st );
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
