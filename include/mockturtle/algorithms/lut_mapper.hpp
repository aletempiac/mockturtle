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
  \file lut_mapper.hpp
  \brief Mapper

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <limits>

#include <fmt/format.h>

#include "../networks/klut.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/tech_library.hpp"
#include "../views/topo_view.hpp"
#include "../views/mapping_view.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/lut_delay_cut.hpp"
#include "cut_enumeration/mf_cut.hpp"
#include "detail/mffc_utils.hpp"

namespace mockturtle
{

/*! \brief Parameters for map.
 *
 * The data structure `map_params` holds configurable parameters
 * with default arguments for `map`.
 */
struct lut_map_params
{
  lut_map_params()
  {
    cut_enumeration_ps.cut_size = 6u;
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = false;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 249. By default,
   * truth table minimization is not performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Required depth for depth relaxation. */
  uint32_t required_delay{ 0u };

  /*! \brief Skip depth round for size optimization. */
  bool skip_delay_round{ false };

  /*! \brief Number of rounds for area flow optimization. */
  uint32_t area_flow_rounds{ 1u };

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t ela_rounds{ 2u };

  /*! \brief Use edge count reduction. */
  bool edge_optimization{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for mapper.
 *
 * The data structure `map_stats` provides data collected by running
 * `map`.
 */
struct lut_map_stats
{
  /*! \brief Area result. */
  uint32_t area{ 0 };
  /*! \brief Worst delay result. */
  uint32_t delay{ 0 };
  /*! \brief Edge result. */
  uint32_t edges{ 0 };

  /*! \brief Runtime for covering. */
  stopwatch<>::duration time_mapping{ 0 };
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Cut enumeration stats. */
  cut_enumeration_stats cut_enumeration_st{};

  /*! \brief Depth and size stats for each round. */
  std::vector<std::string> round_stats{};

  void report() const
  {
    for ( auto const& stat : round_stats )
    {
      std::cout << stat;
    }
    std::cout << fmt::format( "[i] Area = {:8d}; Delay = {:8d}; Edge = {:8d};\n", area, delay, edges );
    std::cout << fmt::format( "[i] Mapping runtime = {:>5.2f} secs\n", to_seconds( time_mapping ) );
    std::cout << fmt::format( "[i] Total runtime   = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

#pragma region cut enumeration
/* cut data */
struct cut_enumeration_lut_cut
{
  uint32_t delay{0};
  float flow{0};
  float cost{0};
};

enum lut_cut_sort_type
{
  DELAY,
  AREA_FLOW,
  AREA,
  NONE
};

template<typename CutType, int MaxCuts>
class lut_cut_set
{
public:
  /*! \brief Standard constructor.
   */
  lut_cut_set()
  {
    clear();
  }

  /*! \brief Clears a cut set.
   */
  void clear()
  {
    _pcend = _pend = _pcuts.begin();
    auto pit = _pcuts.begin();
    for ( auto& c : _cuts )
    {
      *pit++ = &c;
    }
  }

  /*! \brief Adds a cut to the end of the set.
   *
   * This function should only be called to create a set of cuts which is known
   * to be sorted and irredundant (i.e., no cut in the set dominates another
   * cut).
   *
   * \param begin Begin iterator to leaf indexes
   * \param end End iterator (exclusive) to leaf indexes
   * \return Reference to the added cut
   */
  template<typename Iterator>
  CutType& add_cut( Iterator begin, Iterator end )
  {
    assert( _pend != _pcuts.end() );

    auto& cut = **_pend++;
    cut.set_leaves( begin, end );

    ++_pcend;
    return cut;
  }

  /*! \brief Checks whether cut is dominates by any cut in the set.
   *
   * \param cut Cut outside of the set
   */
  bool is_dominated( CutType const& cut ) const
  {
    return std::find_if( _pcuts.begin(), _pcend, [&cut]( auto const* other ) { return other->dominates( cut ); } ) != _pcend;
  }

  static bool sort_delay( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{0.005f};
    if ( c1->data.delay < c2->data.delay )
      return true;
    if ( c1->data.delay > c2->data.delay )
      return false;
    if ( c1.size() < c2.size() )
      return true;
    if ( c1.size() > c2.size() )
      return false;
    return c1->data.flow < c2->data.flow - eps;
  }

  static bool sort_area_flow( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{0.005f};
    if ( c1->data.flow < c2->data.flow - eps )
      return true;
    if ( c1->data.flow > c2->data.flow + eps )
      return false;
    if ( c1->data.delay < c2->data.delay )
      return true;
    if ( c1->data.delay > c2->data.delay )
      return false;
    return c1.size() < c2.size();
  }

  /*! \brief Inserts a cut into a set.
   *
   * This method will insert a cut into a set and maintain an order.  Before the
   * cut is inserted into the correct position, it will remove all cuts that are
   * dominated by `cut`.
   *
   * If `cut` is dominated by any of the cuts in the set, it will still be
   * inserted.  The caller is responsible to check whether `cut` is dominated
   * before inserting it into the set.
   *
   * \param cut Cut to insert.
   */
  void insert( CutType const& cut, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    /* remove elements that are dominated by new cut */
    _pcend = _pend = std::stable_partition( _pcuts.begin(), _pend, [&cut]( auto const* other ) { return !cut.dominates( *other ); } );

    /* insert cut in a sorted way */
    typename std::array<CutType*, MaxCuts>::iterator ipos = _pcuts.begin();

    if ( sort == lut_cut_sort_type::DELAY )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_delay( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::AREA_FLOW )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area_flow( *a, *b ); } );
    }
    else /* NONE */
    {
       ipos == _pend;
    }

    /* too many cuts, we need to remove one */
    if ( _pend == _pcuts.end() )
    {
      /* cut to be inserted is worse than all the others, return */
      if ( ipos == _pend )
      {
        return;
      }
      else
      {
        /* remove last cut */
        --_pend;
        --_pcend;
      }
    }

    /* copy cut */
    auto& icut = *_pend;
    icut->set_leaves( cut.begin(), cut.end() );
    icut->data() = cut.data();

    if ( ipos != _pend )
    {
      auto it = _pend;
      while ( it > ipos )
      {
        std::swap( *it, *( it - 1 ) );
        --it;
      }
    }

    /* update iterators */
    _pcend++;
    _pend++;
  }

  /*! \brief Begin iterator (constant).
   *
   * The iterator will point to a cut pointer.
   */
  auto begin() const { return _pcuts.begin(); }

  /*! \brief End iterator (constant). */
  auto end() const { return _pcend; }

  /*! \brief Begin iterator (mutable).
   *
   * The iterator will point to a cut pointer.
   */
  auto begin() { return _pcuts.begin(); }

  /*! \brief End iterator (mutable). */
  auto end() { return _pend; }

  /*! \brief Number of cuts in the set. */
  auto size() const { return _pcend - _pcuts.begin(); }

  /*! \brief Returns reference to cut at index.
   *
   * This function does not return the cut pointer but dereferences it and
   * returns a reference.  The function does not check whether index is in the
   * valid range.
   *
   * \param index Index
   */
  auto const& operator[]( uint32_t index ) const { return *_pcuts[index]; }

  /*! \brief Returns the best cut, i.e., the first cut.
   */
  auto const& best() const { return *_pcuts[0]; }

  /*! \brief Updates the best cut.
   *
   * This method will set the cut at index `index` to be the best cut.  All
   * cuts before `index` will be moved one position higher.
   *
   * \param index Index of new best cut
   */
  void update_best( uint32_t index )
  {
    auto* best = _pcuts[index];
    for ( auto i = index; i > 0; --i )
    {
      _pcuts[i] = _pcuts[i - 1];
    }
    _pcuts[0] = best;
  }

  /*! \brief Resize the cut set, if it is too large.
   *
   * This method will resize the cut set to `size` only if the cut set has more
   * than `size` elements.  Otherwise, the size will remain the same.
   */
  void limit( uint32_t size )
  {
    if ( std::distance( _pcuts.begin(), _pend ) > static_cast<long>( size ) )
    {
      _pcend = _pend = _pcuts.begin() + size;
    }
  }

  /*! \brief Prints a cut set. */
  friend std::ostream& operator<<( std::ostream& os, lut_cut_set const& set )
  {
    for ( auto const& c : set )
    {
      os << *c << "\n";
    }
    return os;
  }

private:
  std::array<CutType, MaxCuts> _cuts;
  std::array<CutType*, MaxCuts> _pcuts;
  typename std::array<CutType*, MaxCuts>::const_iterator _pcend{_pcuts.begin()};
  typename std::array<CutType*, MaxCuts>::iterator _pend{_pcuts.begin()};
};

template<typename Ntk, bool ComputeTruth>
struct lut_network_cuts
{
public:
  static constexpr uint32_t max_cut_num = 250;
  using cut_t = cut_type<ComputeTruth, cut_enumeration_lut_cut>;
  using cut_set_t = lut_cut_set<cut_t, max_cut_num>;
  static constexpr bool compute_truth = ComputeTruth;
  using node = typename Ntk::node;

public:
  explicit lut_network_cuts( Ntk const& ntk, cut_enumeration_params const& ps, cut_enumeration_stats& st )
      : _ntk( ntk ), _ps( ps ), _st( st ), _cuts( _ntk.size() )
  {
    assert( _ps.cut_limit < max_cut_num && "cut_limit exceeds the compile-time limit for the maximum number of cuts" );

    kitty::dynamic_truth_table zero( 0u ), proj( 1u );
    kitty::create_nth_var( proj, 0u );

    _truth_tables.insert( zero );
    _truth_tables.insert( proj );
  }

public:
  /*! \brief Computes the cuts for each node in the network */
  void compute_cuts( lut_cut_sort_type sort = lut_cut_sort_type::DELAY )
  {
    _ntk.foreach_node( [this, &sort]( auto n ) {
      compute_cuts( n, sort );
    } );
  }

  /*! \brief Computes the cuts of one node in the network */
  void compute_cuts( node const& n, lut_cut_sort_type sort )
  {
    const auto index = _ntk.node_to_index( n );

    if ( _ps.very_verbose )
    {
      std::cout << fmt::format( "[i] compute cut for node at index {}\n", index );
    }

    if ( _ntk.is_constant( n ) )
    {
      add_zero_cut( index );
    }
    else if ( _ntk.is_pi( n ) )
    {
      add_unit_cut( index );
    }
    else
    {
      if constexpr ( Ntk::min_fanin_size == 2 && Ntk::max_fanin_size == 2 )
      {
        merge_cuts2( index, sort );
      }
      else
      {
        merge_cuts( index, sort );
      }
    }
  }

  /*! \brief Returns the cut set of a node */
  cut_set_t& cuts( uint32_t node_index ) { return _cuts[node_index]; }

  /*! \brief Returns the cut set of a node */
  cut_set_t const& cuts( uint32_t node_index ) const { return _cuts[node_index]; }

  /*! \brief Returns the truth table of a cut */
  template<bool enabled = ComputeTruth, typename = std::enable_if_t<std::is_same_v<Ntk, Ntk> && enabled>>
  auto truth_table( cut_t const& cut ) const
  {
    return _truth_tables[cut->func_id];
  }

  /*! \brief Returns the total number of tuples that were tried to be merged */
  auto total_tuples() const
  {
    return _total_tuples;
  }

  /*! \brief Returns the total number of cuts in the database. */
  auto total_cuts() const
  {
    return _total_cuts;
  }

  /*! \brief Returns the number of nodes for which cuts are computed */
  auto nodes_size() const
  {
    return _cuts.size();
  }

  /* compute positions of leave indices in cut `sub` (subset) with respect to
   * leaves in cut `sup` (super set).
   *
   * Example:
   *   compute_truth_table_support( {1, 3, 6}, {0, 1, 2, 3, 6, 7} ) = {1, 3, 4}
   */
  std::vector<uint8_t> compute_truth_table_support( cut_t const& sub, cut_t const& sup ) const
  {
    std::vector<uint8_t> support;
    support.reserve( sub.size() );

    auto itp = sup.begin();
    for ( auto i : sub )
    {
      itp = std::find( itp, sup.end(), i );
      support.push_back( static_cast<uint8_t>( std::distance( sup.begin(), itp ) ) );
    }

    return support;
  }

  /*! \brief Inserts a truth table into the truth table cache.
   *
   * This message can be used when manually adding or modifying cuts from the
   * cut sets.
   *
   * \param tt Truth table to add
   * \return Literal id from the truth table store
   */
  uint32_t insert_truth_table( kitty::dynamic_truth_table const& tt )
  {
    return _truth_tables.insert( tt );
  }

private:
  void add_zero_cut( uint32_t index )
  {
    auto& cut = _cuts[index].add_cut( &index, &index ); /* fake iterator for emptyness */

    if constexpr ( ComputeTruth )
    {
      cut->func_id = 0;
    }
  }

  void add_unit_cut( uint32_t index )
  {
    auto& cut = _cuts[index].add_cut( &index, &index + 1 );

    if constexpr ( ComputeTruth )
    {
      cut->func_id = 2;
    }
  }

  void compute_cut_data( cut_t& cut, node const& n )
  {
    uint32_t delay{0};
    float flow = cut->data.cost = cut.size() < 2 ? 0.0f : 1.0f;

    for ( auto leaf : cut )
    {
      const auto& best_leaf_cut = _cuts[leaf][0];
      delay = std::max( delay, best_leaf_cut->data.delay );
      flow += best_leaf_cut->data.flow;
    }

    cut->data.delay = 1 + delay;
    cut->data.flow = flow / _ntk.fanout_size( n );
  }

  uint32_t compute_truth_table( uint32_t index, std::vector<cut_t const*> const& vcuts, cut_t& res )
  {
    stopwatch t( _st.time_truth_table );

    std::vector<kitty::dynamic_truth_table> tt( vcuts.size() );
    auto i = 0;
    for ( auto const& cut : vcuts )
    {
      tt[i] = kitty::extend_to( _truth_tables[( *cut )->func_id], res.size() );
      const auto supp = compute_truth_table_support( *cut, res );
      kitty::expand_inplace( tt[i], supp );
      ++i;
    }

    auto tt_res = _ntk.compute( _ntk.index_to_node( index ), tt.begin(), tt.end() );

    if ( _ps.minimize_truth_table )
    {
      const auto support = kitty::min_base_inplace( tt_res );
      if ( support.size() != res.size() )
      {
        auto tt_res_shrink = shrink_to( tt_res, static_cast<unsigned>( support.size() ) );
        std::vector<uint32_t> leaves_before( res.begin(), res.end() );
        std::vector<uint32_t> leaves_after( support.size() );

        auto it_support = support.begin();
        auto it_leaves = leaves_after.begin();
        while ( it_support != support.end() )
        {
          *it_leaves++ = leaves_before[*it_support++];
        }
        res.set_leaves( leaves_after.begin(), leaves_after.end() );
        return _truth_tables.insert( tt_res_shrink );
      }
    }

    return _truth_tables.insert( tt_res );
  }

  void merge_cuts2( uint32_t index, lut_cut_sort_type sort )
  {
    const auto fanin = 2;

    uint32_t pairs{1};
    _ntk.foreach_fanin( _ntk.index_to_node( index ), [this, &pairs]( auto child, auto i ) {
      lcuts[i] = &_cuts[_ntk.node_to_index( _ntk.get_node( child ) )];
      pairs *= static_cast<uint32_t>( lcuts[i]->size() );
    } );
    lcuts[2] = &_cuts[index];
    auto& rcuts = *lcuts[fanin];
    rcuts.clear();

    cut_t new_cut;

    std::vector<cut_t const*> vcuts( fanin );

    _total_tuples += pairs;
    for ( auto const& c1 : *lcuts[0] )
    {
      for ( auto const& c2 : *lcuts[1] )
      {
        if ( !c1->merge( *c2, new_cut, _ps.cut_size ) )
        {
          continue;
        }

        if ( rcuts.is_dominated( new_cut ) )
        {
          continue;
        }

        if constexpr ( ComputeTruth )
        {
          vcuts[0] = c1;
          vcuts[1] = c2;
          new_cut->func_id = compute_truth_table( index, vcuts, new_cut );
        }

        compute_cut_data( new_cut, _ntk.index_to_node( index ) );

        rcuts.insert( new_cut, sort );
      }
    }

    /* limit the maximum number of cuts */
    rcuts.limit( _ps.cut_limit - 1 );

    _total_cuts += rcuts.size();

    if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
    {
      add_unit_cut( index );
    }
  }

  void merge_cuts( uint32_t index, lut_cut_sort_type sort )
  {
    uint32_t pairs{1};
    std::vector<uint32_t> cut_sizes;
    _ntk.foreach_fanin( _ntk.index_to_node( index ), [this, &pairs, &cut_sizes]( auto child, auto i ) {
      lcuts[i] = &_cuts[_ntk.node_to_index( _ntk.get_node( child ) )];
      cut_sizes.push_back( static_cast<uint32_t>( lcuts[i]->size() ) );
      pairs *= cut_sizes.back();
    } );

    const auto fanin = cut_sizes.size();
    lcuts[fanin] = &_cuts[index];

    auto& rcuts = *lcuts[fanin];

    if ( fanin > 1 && fanin <= _ps.fanin_limit )
    {
      rcuts.clear();

      cut_t new_cut, tmp_cut;

      std::vector<cut_t const*> vcuts( fanin );

      _total_tuples += pairs;
      foreach_mixed_radix_tuple( cut_sizes.begin(), cut_sizes.end(), [&]( auto begin, auto end ) {
        auto it = vcuts.begin();
        auto i = 0u;
        while ( begin != end )
        {
          *it++ = &( ( *lcuts[i++] )[*begin++] );
        }

        if ( !vcuts[0]->merge( *vcuts[1], new_cut, _ps.cut_size ) )
        {
          return true; /* continue */
        }

        for ( i = 2; i < fanin; ++i )
        {
          tmp_cut = new_cut;
          if ( !vcuts[i]->merge( tmp_cut, new_cut, _ps.cut_size ) )
          {
            return true; /* continue */
          }
        }

        if ( rcuts.is_dominated( new_cut ) )
        {
          return true; /* continue */
        }

        if constexpr ( ComputeTruth )
        {
          new_cut->func_id = compute_truth_table( index, vcuts, new_cut );
        }

        comput_cut_data( new_cut, index );

        rcuts.insert( new_cut, sort );

        return true;
      } );

      /* limit the maximum number of cuts */
      rcuts.limit( _ps.cut_limit - 1 );
    } else if ( fanin == 1 ) {
      rcuts.clear();

      for ( auto const& cut : *lcuts[0] ) {
        cut_t new_cut = *cut;

        if constexpr ( ComputeTruth )
        {
          new_cut->func_id = compute_truth_table( index, {cut}, new_cut );
        }

        compute_cut_data( new_cut, _ntk.index_to_node( index ) );

        rcuts.insert( new_cut, sort );
      }

      /* limit the maximum number of cuts */
      rcuts.limit( _ps.cut_limit - 1 );
    }

    _total_cuts += static_cast<uint32_t>( rcuts.size() );

    add_unit_cut( index );
  }

private:
  Ntk const& _ntk;
  cut_enumeration_params const& _ps;
  cut_enumeration_stats& _st;

  /* compressed representation of cuts */
  std::vector<cut_set_t> _cuts;

  /* node cut computation container */
  std::array<cut_set_t*, Ntk::max_fanin_size + 1> lcuts;

  /* cut truth tables */
  truth_table_cache<kitty::dynamic_truth_table> _truth_tables;

  /* statistics */
  uint32_t _total_tuples{};
  std::size_t _total_cuts{};
};
#pragma endregion

#pragma region LUT mapper
struct node_lut
{
  /* arrival time at node output */
  uint32_t arrival;
  /* required time at node output */
  uint32_t required;
  /* area of the best match */
  uint32_t area;
  /* edge count of the best match */
  uint32_t edges;

  /* number of references in the cover 0: pos, 1: neg, 2: pos+neg */
  uint32_t map_refs;
  /* references estimation */
  float est_refs;
  /* area flow */
  float flows;
  /* edge flow */
  float edge_flows;
};

template<class Ntk, bool StoreFunction>
class lut_map_impl
{
public:
  using network_cuts_t = lut_network_cuts<Ntk, StoreFunction>;
  using cut_t = typename network_cuts_t::cut_t;

public:
  explicit lut_map_impl( Ntk& ntk, lut_map_params const& ps, lut_map_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        node_match( ntk.size() ),
        cuts( ntk, ps.cut_enumeration_ps, st.cut_enumeration_st )
  {}

  void run()
  {
    stopwatch t( st.time_mapping );

    /* compute and save topological order */
    top_order.reserve( ntk.size() );
    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      top_order.push_back( n );
    } );

    /* init the data structure */
    init_nodes();

    /* compute cuts */
    cuts.compute_cuts();

    /* compute mapping for depth */
    if ( !ps.skip_delay_round )
    {
      compute_mapping<false>();
    }

    /* compute mapping using global area flow */
    while ( iteration < ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      compute_mapping<true>();
    }

    /* compute mapping using exact area */
    while ( iteration < ps.ela_rounds + ps.area_flow_rounds + 1 )
    {
      compute_required_time();
      compute_mapping_exact();
    }

    /* generate the output network */
    derive_mapping();
  }

private:
  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n, auto ) {
      const auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      node_data.est_refs = static_cast<float>( ntk.fanout_size( n ) );

      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        /* all terminals have flow 1.0 */
        node_data.flows = 0.0f;
        node_data.edge_flows = 0.0f;
        node_data.arrival = 0.0f;
      }
    } );
  }

  template<bool DO_AREA>
  void compute_mapping()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        continue;
      }

      compute_best_cut<DO_AREA>( n );
    }

    uint32_t area_old = area;
    set_mapping_refs<false>();

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats{};

      if constexpr ( DO_AREA )
      {
        stats << fmt::format( "[i] AreaFlow : Delay = {:8d}  Area = {:8d}  Edges = {:8d}\n", delay, area, edges );
      }
      else
      {
        stats << fmt::format( "[i] Delay    : Delay = {:8d}  Area = {:8d}  Edges = {:8d}\n", delay, area, edges );
      }
      st.round_stats.push_back( stats.str() );
    }
  }

  void compute_mapping_exact()
  {
    for ( auto const& n : top_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      compute_best_cut_exact( n );
    }

    uint32_t area_old = area;
    set_mapping_refs<true>();

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats{};
      stats << fmt::format( "[i] Area     : Delay = {:8d}  Area = {:8d}  Edges = {:8d}\n", delay, area, edges );
      st.round_stats.push_back( stats.str() );
    }
  }

  template<bool ELA>
  void set_mapping_refs()
  {
    const auto coef = 1.0f / ( 2.0f + ( iteration + 1 ) * ( iteration + 1 ) );

    if constexpr ( !ELA )
    {
      for ( auto i = 0u; i < node_match.size(); ++i )
      {
        node_match[i].map_refs = 0u;
      }
    }

    /* compute the current worst delay and update the mapping refs */
    delay = 0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );

      delay = std::max( delay, node_match[index].arrival );

      if constexpr ( !ELA )
      {
        node_match[index].map_refs++;
      }
    } );

    /* compute current area and update mapping refs in top-down order */
    area = 0;
    edges = 0;
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      const auto index = ntk.node_to_index( *it );
      auto& node_data = node_match[index];

      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
      {
        continue;
      }

      /* continue if not referenced in the cover */
      if ( node_match[index].map_refs == 0u )
        continue;

      if constexpr ( !ELA )
      {
        for ( auto const leaf : cuts.cuts( index )[0] )
        {
          node_match[leaf].map_refs++;
        }
      }
      ++area;
      edges += cuts.cuts( index )[0].size();
    }

    /* blend estimated references */
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      node_match[i].est_refs = coef * node_match[i].est_refs + ( 1.0f - coef ) * std::max( 1.0f, static_cast<float>( node_match[i].map_refs ) );
    }

    ++iteration;
  }

  void compute_required_time()
  {
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required = UINT32_MAX;
    }

    /* return in case of `skip_delay_round` */
    if ( iteration == 0 )
      return;

    auto required = delay;

    if ( ps.required_delay != 0 )
    {
      /* Global target time constraint */
      if ( ps.required_delay < delay )
      {
        if ( !ps.skip_delay_round && iteration == 1 )
          std::cerr << fmt::format( "[i] MAP WARNING: cannot meet the target required time of {:.2f}", ps.required_delay ) << std::endl;
      }
      else
      {
        required = ps.required_delay;
      }
    }

    /* set the required time at POs */
    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      node_match[index].required = required;
    } );

    /* propagate required time to the PIs */
    for ( auto it = top_order.rbegin(); it != top_order.rend(); ++it )
    {
      if ( ntk.is_pi( *it ) || ntk.is_constant( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );

      if ( node_match[index].map_refs == 0 )
        continue;

      for ( auto leaf : cuts.cuts( index )[0] )
      {
        node_match[leaf].required = std::min( node_match[leaf].required, node_match[index].required - 1 );
      }
    }
  }

  template<bool DO_AREA>
  void compute_best_cut( node<Ntk> const& n )
  {
    uint32_t best_arrival = UINT32_MAX;
    float best_area_flow = std::numeric_limits<float>::max();
    float best_edge_flow = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];

    /* foreach cut */
    for ( auto& cut : cuts.cuts( index ) )
    {
      /* ignore trivial cut */
      if ( cut->size() == 1 && ( *cut->begin() ) == index )
      {
        ++cut_index;
        continue;
      }

      uint32_t worst_arrival = 0;
      float flow = 0;
      float edge_flow = 0;

      for ( auto leaf : *cut )
      {
        worst_arrival = std::max( worst_arrival, node_match[leaf].arrival + 1 );
        flow += node_match[leaf].flows;
        edge_flow += node_match[leaf].edge_flows;
      }

      float area_local = 1 + flow;
      float edge_local = cut->size() + edge_flow;

      if constexpr ( DO_AREA )
      {
        if ( worst_arrival > node_data.required )
        {
          ++cut_index;
          continue;
        }
      }

      bool result = false;
      if ( ps.edge_optimization )
      {
        result = compare_map_edge<DO_AREA>( worst_arrival, best_arrival, area_local, best_area_flow, edge_local, best_edge_flow, cut->size(), best_size );
      }
      else
      {
        result = compare_map<DO_AREA>( worst_arrival, best_arrival, area_local, best_area_flow, cut->size(), best_size );
      }

      if ( result )
      {
        best_arrival = worst_arrival;
        best_area_flow = area_local;
        best_edge_flow = edge_local;
        best_size = cut->size();
        best_cut = cut_index;
      }
      ++cut_index;
    }

    node_data.flows = best_area_flow / node_data.est_refs;
    node_data.edge_flows = best_edge_flow / node_data.est_refs;
    node_data.arrival = best_arrival;

    if ( best_cut != 0 )
    {
      cuts.cuts( index ).update_best( best_cut );
    }
  }

  void compute_best_cut_exact( node<Ntk> const& n )
  {
    uint32_t best_arrival = UINT32_MAX;
    uint32_t best_exact_area = UINT32_MAX;
    uint32_t best_exact_edge = UINT32_MAX;
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint8_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];

    /* recursively deselect the best cut shared between
     * the two phases if in use in the cover */
    if ( node_data.map_refs )
    {
      cut_deref( cuts.cuts( index )[0] );
    }

    /* foreach cut */
    for ( auto& cut : cuts.cuts( index ) )
    {
      /* ignore trivial cut */
      if ( cut->size() == 1 && ( *cut->begin() ) == index )
      {
        ++cut_index;
        continue;
      }

      uint32_t area_exact = cut_ref( *cut );
      uint32_t edge_exact = cut_edge_deref( *cut );
      uint32_t worst_arrival = 0;

      for ( auto l : *cut )
      {
        worst_arrival = std::max( worst_arrival, node_match[l].arrival + 1 );
      }

      if ( worst_arrival > node_data.required )
      {
        ++cut_index;
        continue;
      }

      bool result = false;
      if ( ps.edge_optimization )
      {
        result = compare_map_edge<true>( worst_arrival, best_arrival, area_exact, best_exact_area, edge_exact, best_exact_edge, cut->size(), best_size );
      }
      else
      {
        result = compare_map<true>( worst_arrival, best_arrival, area_exact, best_exact_area, cut->size(), best_size );
      }

      if ( result )
      {
        best_arrival = worst_arrival;
        best_exact_area = area_exact;
        best_exact_edge = edge_exact;
        best_size = cut->size();
        best_cut = cut_index;
      }

      ++cut_index;
    }

    node_data.flows = best_exact_area;
    node_data.arrival = best_arrival;

    if ( best_cut != 0 )
    {
      cuts.cuts( index ).update_best( best_cut );
    }

    if ( node_data.map_refs )
    {
      cut_ref( cuts.cuts( index )[0] );
    }
  }

  uint32_t cut_ref( cut_t const& cut )
  {
    uint32_t count = 1;

    for ( auto leaf : cut )
    {
      if ( ntk.is_pi( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( node_match[leaf].map_refs++ == 0u )
      {
        count += cut_ref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  uint32_t cut_deref( cut_t const& cut )
  {
    uint32_t count = 1;

    for ( auto leaf : cut )
    {
      if ( ntk.is_pi( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( --node_match[leaf].map_refs == 0u )
      {
        count += cut_deref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  uint32_t cut_edge_ref( cut_t const& cut )
  {
    uint32_t count = cut.size();

    for ( auto leaf : cut )
    {
      if ( ntk.is_pi( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( node_match[leaf].map_refs++ == 0u )
      {
        count += cut_ref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  uint32_t cut_edge_deref( cut_t const& cut )
  {
    uint32_t count = cut.size();

    for ( auto leaf : cut )
    {
      if ( ntk.is_pi( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( --node_match[leaf].map_refs == 0u )
      {
        count += cut_deref( cuts.cuts( leaf )[0] );
      }
    }
    return count;
  }

  void derive_mapping()
  {
    ntk.clear_mapping();

    for ( auto const& n : top_order )
    {
      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
        continue;

      const auto index = ntk.node_to_index( n );
      if ( node_match[index].map_refs == 0 )
        continue;

      std::vector<node<Ntk>> nodes;
      auto const& best_cut = cuts.cuts( index ).best();

      for ( auto const& l : best_cut )
      {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      if constexpr ( StoreFunction )
      {
        ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
      }
    }

    st.area = area;
    st.delay = delay;
    st.edges = edges;
  }

  template<bool DO_AREA>
  inline bool compare_map_edge( uint32_t arrival, uint32_t best_arrival, float area_flow, float best_area_flow, float edge_flow, float best_edge_flow, uint32_t size, uint32_t best_size )
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
      else if ( edge_flow < best_edge_flow - epsilon )
      {
        return true;
      }
      else if ( edge_flow > best_edge_flow + epsilon )
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
      else if ( edge_flow < best_edge_flow - epsilon )
      {
        return true;
      }
      else if ( edge_flow > best_edge_flow + epsilon )
      {
        return false;
      }
    }
    if ( size < best_size )
    {
      return true;
    }
    return false;
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
    return false;
  }

private:
  Ntk& ntk;
  lut_map_params const& ps;
  lut_map_stats& st;

  uint32_t iteration{ 0 };       /* current mapping iteration */
  uint32_t delay{ 0 };          /* current delay of the mapping */
  uint32_t area{ 0 };           /* current area of the mapping */
  uint32_t edges{ 0 };           /* current edges of the mapping */
  const float epsilon{ 0.005f }; /* epsilon */

  std::vector<node<Ntk>> top_order;
  std::vector<node_lut> node_match;
  network_cuts_t cuts;
};
#pragma endregion

} /* namespace detail */

/*! \brief LUT mapper.
 *
 * This function implements a LUT mapping algorithm. It is controlled by a
 * template argument `CutData` (defaulted to `cut_enumeration_tech_map_cut`).
 * The argument is similar to the `CutData` argument in `cut_enumeration`, which can
 * specialize the cost function to select priority cuts and store additional data.
 * The default argument gives priority firstly to the cut size, then delay, and lastly
 * to area flow. Thus, it is more suited for delay-oriented mapping.
 * The type passed as `CutData` must implement the following four fields:
 *
 * - `uint32_t delay`
 * - `float flow`
 * - `uint8_t match_index`
 * - `bool ignore`
 *
 * See `include/mockturtle/algorithms/cut_enumeration/cut_enumeration_tech_map_cut.hpp`
 * for one example of a CutData type that implements the cost function that is used in
 * the technology mapper.
 * 
 * The function takes the size of the cuts in the template parameter `CutSize`.
 *
 * The function returns a k-LUT network. Each LUT abstacts a gate of the technology library.
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
 *
 * \param ntk Network
 * \param library Technology library
 * \param ps Mapping params
 * \param pst Mapping statistics
 * 
 * The implementation of this algorithm was inspired by the
 * mapping command ``map`` in ABC.
 */
template<class Ntk, bool StoreFunction = false>
void lut_map( Ntk& ntk, lut_map_params const& ps = {}, lut_map_stats* pst = nullptr )
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
  static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );

  lut_map_stats st;
  detail::lut_map_impl<Ntk, StoreFunction> p( ntk, ps, st );
  p.run();

  st.time_total = st.time_mapping + st.cut_enumeration_st.time_total;
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
