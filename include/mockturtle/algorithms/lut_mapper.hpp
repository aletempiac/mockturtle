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
  \file lut_mapper.hpp
  \brief LUT mapper

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <fmt/format.h>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "cleanup.hpp"
#include "../networks/klut.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/cuts.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/truth_table_cache.hpp"
#include "../views/choice_view.hpp"
#include "../views/mapping_view.hpp"
#include "../views/mffc_view.hpp"
#include "../views/topo_view.hpp"
#include "cut_enumeration.hpp"
#include "simulation.hpp"

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
    cut_enumeration_ps.cut_limit = 8u;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 8. The maximum value
   * is 249. By default, truth table minimization
   * is performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Do area-oriented mapping. */
  bool area_oriented_mapping{ false };

  /*! \brief Required depth for depth relaxation. */
  uint32_t required_delay{ 0u };

  /*! \brief Required depth relaxation ratio (%). */
  uint32_t relax_required{ 0u };

  /*! \brief Recompute cuts at each step. */
  bool recompute_cuts{ true };

  /*! \brief Number of rounds for area sharing optimization. */
  uint32_t area_share_rounds{ 2u };

  /*! \brief Number of rounds for area flow optimization. */
  uint32_t area_flow_rounds{ 1u };

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t ela_rounds{ 2u };

  /*! \brief Use edge count reduction. */
  bool edge_optimization{ false };

  /*! \brief Try to expand the cuts. */
  bool cut_expansion{ true };

  /*! \brief Remove the cuts that are contained in others */
  bool remove_dominated_cuts{ true };

  /*! \brief Decomposes multi-input ANDs/XORs */
  bool multi_decomposition{ false };

  /*! \brief Maximum large cut size for multi-input ANDs */
  uint32_t large_cut_size{ 32 };

  /*! \brief Maps by collapsing MFFCs */
  bool collapse_mffcs{ false };

  /*! \brief Maximum number variables for cost function caching */
  uint32_t cost_cache_vars{ 3u };

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
    std::cout << fmt::format( "[i] Total runtime           = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

#pragma region cut set
/* cut data */
struct cut_enumeration_lut_cut
{
  uint32_t delay{ 0 };
  uint32_t lut_area{ 0 };
  uint32_t lut_delay{ 0 };
  float area_flow{ 0 };
  float edge_flow{ 0 };
};

struct cut_enumeration_adaptive_cut
{
  uint32_t delay{ 0 };
  uint32_t lut_area{ 0 };
  uint32_t lut_delay{ 0 };
  uint32_t lut_edges{ 0 };
  float area_flow{ 0 };
  float edge_flow{ 0 };
  uint64_t connections{ 0 };
};

/* wide cut data */
template<uint32_t CutSize, uint32_t LutSize>
struct cut_enumeration_wide_cut
{
  uint32_t delay{ 0 };
  uint32_t lut_area{ 0 };
  uint32_t lut_delay{ 0 };
  uint32_t lut_edges{ 0 };
  float area_flow{ 0 };
  float edge_flow{ 0 };
  
  std::vector<uint64_t> lut_roots;
  std::vector<uint32_t> delays;
};

enum class lut_cut_sort_type
{
  DELAY,
  DELAY2,
  AREA,
  AREA2,
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

  /*! \brief Assignment operator.
   */
  lut_cut_set& operator=( lut_cut_set const& other )
  {
    if ( this != &other )
    {
      _pcend = _pend = _pcuts.begin();

      auto it = other.begin();
      while( it != other.end() )
      {
        **_pend++ = **it++;
        ++_pcend;
      }
    }

    return *this;
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
    constexpr auto eps{ 0.005f };
    if ( c1->data.delay < c2->data.delay )
      return true;
    if ( c1->data.delay > c2->data.delay )
      return false;
    if ( c1.size() < c2.size() )
      return true;
    if ( c1.size() > c2.size() )
      return false;
    if ( c1->data.area_flow < c2->data.area_flow - eps )
      return true;
    if ( c1->data.area_flow > c2->data.area_flow + eps )
      return false;
    return c1->data.edge_flow < c2->data.edge_flow - eps;
  }

  static bool sort_delay2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1->data.delay < c2->data.delay )
      return true;
    if ( c1->data.delay > c2->data.delay )
      return false;
    if ( c1->data.area_flow < c2->data.area_flow - eps )
      return true;
    if ( c1->data.area_flow > c2->data.area_flow + eps )
      return false;
    if ( c1->data.edge_flow < c2->data.edge_flow - eps )
      return true;
    if ( c1->data.edge_flow > c2->data.edge_flow + eps )
      return false;
    return c1.size() < c2.size();
  }

  static bool sort_area( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1->data.area_flow < c2->data.area_flow - eps )
      return true;
    if ( c1->data.area_flow > c2->data.area_flow + eps )
      return false;
    if ( c1->data.delay < c2->data.delay )
      return true;
    if ( c1->data.delay > c2->data.delay )
      return false;
    return c1.size() < c2.size();
  }

  static bool sort_area2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1->data.area_flow < c2->data.area_flow - eps )
      return true;
    if ( c1->data.area_flow > c2->data.area_flow + eps )
      return false;
    if ( c1->data.edge_flow < c2->data.edge_flow - eps )
      return true;
    if ( c1->data.edge_flow > c2->data.edge_flow + eps )
      return false;
    if ( c1.size() < c2.size() )
      return true;
    if ( c1.size() > c2.size() )
      return false;
    return c1->data.delay < c2->data.delay;
  }

  /*! \brief Compare two cuts using sorting functions.
   *
   * This method compares two cuts using a sorting function.
   *
   * \param cut1 first cut.
   * \param cut2 second cut.
   * \param sort sorting function.
   */
  static bool compare( CutType const& cut1, CutType const& cut2, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    if ( sort == lut_cut_sort_type::DELAY )
    {
      return sort_delay( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::DELAY2 )
    {
      return sort_delay2( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::AREA )
    {
      return sort_area( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::AREA2 )
    {
      return sort_area2( cut1, cut2 );
    }
    else
    {
      return false;
    }
  }

  /*! \brief Inserts a cut into a set without checking dominance.
   *
   * This method will insert a cut into a set and maintain an order.  This
   * method doesn't remove the cuts that are dominated by `cut`.
   *
   * If `cut` is dominated by any of the cuts in the set, it will still be
   * inserted.  The caller is responsible to check whether `cut` is dominated
   * before inserting it into the set.
   *
   * \param cut Cut to insert.
   * \param sort Cut prioritization function.
   */
  void simple_insert( CutType const& cut, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    /* insert cut in a sorted way */
    typename std::array<CutType*, MaxCuts>::iterator ipos = _pcuts.begin();

    if ( sort == lut_cut_sort_type::DELAY )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_delay( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::DELAY2 )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_delay2( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::AREA )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::AREA2 )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area2( *a, *b ); } );
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

  /*! \brief Inserts a cut into a set.
   *
   * This method will insert a cut into a set and maintain an order.  Before the
   * cut is inserted into the correct position, it will remove all cuts that are
   * dominated by `cut`. Variable `skip0` tell to skip the dominance check on
   * cut zero.
   *
   * If `cut` is dominated by any of the cuts in the set, it will still be
   * inserted.  The caller is responsible to check whether `cut` is dominated
   * before inserting it into the set.
   *
   * \param cut Cut to insert.
   * \param skip0 Skip dominance check on cut zero.
   * \param sort Cut prioritization function.
   */
  void insert( CutType const& cut, bool skip0 = false, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    auto begin = _pcuts.begin();

    if ( skip0 && _pend != _pcuts.begin() )
      ++begin;

    /* remove elements that are dominated by new cut */
    _pcend = _pend = std::stable_partition( begin, _pend, [&cut]( auto const* other ) { return !cut.dominates( *other ); } );

    /* insert cut in a sorted way */
    simple_insert( cut, sort );
  }

  /*! \brief Replaces a cut of the set.
   *
   * This method replaces the cut at position `index` in the set by `cut`
   * and maintains the cuts order. The function does not check whether
   * index is in the valid range.
   *
   * \param index Index of the cut to replace.
   * \param cut Cut to insert.
   */
  void replace( uint32_t index, CutType const& cut )
  {
    *_pcuts[index] = cut;
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
  typename std::array<CutType*, MaxCuts>::const_iterator _pcend{ _pcuts.begin() };
  typename std::array<CutType*, MaxCuts>::iterator _pend{ _pcuts.begin() };
};
#pragma endregion

#pragma region Adaptive cuts bucket
template<typename CutType, unsigned MaxCuts>
class adaptive_cut_set
{
public:
  /*! \brief Standard constructor.
   */
  adaptive_cut_set()
  {
    clear();
  }

  /*! \brief Clears the buckets.
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

  template<typename Iterator>
  CutType& add_cut( CutType const& cut_copy )
  {
    assert( _pend != _pcuts.end() );

    auto& cut = **_pend++;
    cut = cut_copy;

    ++_pcend;
    return cut;
  }

  static bool sort_delay( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1.data.delay < c2.data.delay )
      return true;
    if ( c1.data.delay > c2.data.delay )
      return false;
    if ( c1.size < c2.size )
      return true;
    if ( c1.size > c2.size )
      return false;
    if ( c1.data.area_flow < c2.data.area_flow - eps )
      return true;
    if ( c1.data.area_flow > c2.data.area_flow + eps )
      return false;
    return c1.data.edge_flow < c2.data.edge_flow - eps;
  }

  static bool sort_delay2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1.data.delay < c2.data.delay )
      return true;
    if ( c1.data.delay > c2.data.delay )
      return false;
    if ( c1.data.area_flow < c2.data.area_flow - eps )
      return true;
    if ( c1.data.area_flow > c2.data.area_flow + eps )
      return false;
    if ( c1.data.edge_flow < c2.data.edge_flow - eps )
      return true;
    if ( c1.data.edge_flow > c2.data.edge_flow + eps )
      return false;
    return c1.size < c2.size;
  }

  static bool sort_area( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1.data.area_flow < c2.data.area_flow - eps )
      return true;
    if ( c1.data.area_flow > c2.data.area_flow + eps )
      return false;
    if ( c1.data.delay < c2.data.delay )
      return true;
    if ( c1.data.delay > c2.data.delay )
      return false;
    return c1.size < c2.size;
  }

  static bool sort_area2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( c1.data.area_flow < c2.data.area_flow - eps )
      return true;
    if ( c1.data.area_flow > c2.data.area_flow + eps )
      return false;
    if ( c1.data.edge_flow < c2.data.edge_flow - eps )
      return true;
    if ( c1.data.edge_flow > c2.data.edge_flow + eps )
      return false;
    if ( c1.size < c2.size )
      return true;
    if ( c1.size > c2.size )
      return false;
    return c1.data.delay < c2.data.delay;
  }

  static bool compare( CutType const& cut1, CutType const& cut2, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    if ( sort == lut_cut_sort_type::DELAY )
    {
      return sort_delay( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::DELAY2 )
    {
      return sort_delay2( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::AREA )
    {
      return sort_area( cut1, cut2 );
    }
    else if ( sort == lut_cut_sort_type::AREA2 )
    {
      return sort_area2( cut1, cut2 );
    }
    else
    {
      return false;
    }
  }

  void insert( CutType const& cut, lut_cut_sort_type sort = lut_cut_sort_type::NONE )
  {
    /* insert cut in a sorted way */
    typename std::array<CutType*, MaxCuts>::iterator ipos = _pcuts.begin();

    if ( sort == lut_cut_sort_type::DELAY )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_delay( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::DELAY2 )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_delay2( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::AREA )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area( *a, *b ); } );
    }
    else if ( sort == lut_cut_sort_type::AREA2 )
    {
      ipos = std::lower_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area2( *a, *b ); } );
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
    **_pend = cut;

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

  auto begin() const { return _pcuts.begin(); }

  auto end() const { return _pcend; }

  auto begin() { return _pcuts.begin(); }

  auto end() { return _pend; }

  auto size() const { return _pcend - _pcuts.begin(); }

  auto const& operator[]( uint32_t index ) const { return *_pcuts[index]; }

  auto const& best() const { return *_pcuts[0]; }

private:
  std::array<CutType, MaxCuts> _cuts;
  std::array<CutType*, MaxCuts> _pcuts;
  typename std::array<CutType*, MaxCuts>::const_iterator _pcend{ _pcuts.begin() };
  typename std::array<CutType*, MaxCuts>::iterator _pend{ _pcuts.begin() };
};
#pragma endregion

#pragma region Cut merger data
template<uint32_t MaxFanin>
struct adaptive_cut_merger
{
  std::array<uint32_t, MaxFanin> cut_pointers;
  cut_enumeration_lut_cut data;
  uint32_t size;
};

template<typename CutType, uint32_t MaxFanin>
struct bin_merger
{
  std::array<CutType, MaxFanin> cuts;
  std::array<CutType*, MaxFanin> pcuts;
  uint32_t size{ 0 };

  CutType& operator[]( uint32_t index ) const { return *pcuts[index]; }

  void clear()
  {
    size = 0;
    uint32_t i;
    for ( i = 0; i < MaxFanin; ++i )
    {
      pcuts[i] = &cuts[i];
    }
  }
};

struct bin_data
{
  uint32_t size{ 0 };
  uint32_t delay{ 0 };
  bool connected{ false };
};
#pragma endregion

#pragma region LUT mapper
struct node_lut
{
  /* required time at node output */
  uint32_t required;
  /* number of references in the cover */
  uint32_t map_refs;
  /* references estimation */
  float est_refs;
};

template<class Ntk, bool StoreFunction, class LUTCostFn>
class lut_map_impl
{
public:
  static constexpr uint32_t max_cut_num = 16;
  static constexpr uint32_t max_cut_size = 6;
  static constexpr uint32_t max_wide_cut_num = 1;
  static constexpr uint32_t max_gate_fanin = 32;
  static constexpr uint32_t max_wide_cut_size = 48;
  static constexpr uint32_t max_cut_data_decomp_bucket = 8;
  using cut_t = cut<max_cut_size, cut_data<StoreFunction, cut_enumeration_lut_cut>>;
  using cut_set_t = lut_cut_set<cut_t, max_cut_num>;
  using wide_cut_t = cut<max_wide_cut_size, cut_data<false, cut_enumeration_wide_cut<max_wide_cut_size, max_cut_size>>>;
  using wide_cut_set_t = lut_cut_set<wide_cut_t, max_wide_cut_num>;
  using adaptive_cut_t = adaptive_cut_merger<max_gate_fanin>;
  using adaptive_cut_set_t = adaptive_cut_set<adaptive_cut_t, max_cut_data_decomp_bucket>;
  using adaptive_buckets_t = std::array<adaptive_cut_set_t, max_wide_cut_size>;
  using bin_merger_t = bin_merger<cut<max_cut_size, cut_data<false, cut_enumeration_adaptive_cut>>, max_gate_fanin>;
  using node = typename Ntk::node;
  using cut_merge_t = typename std::array<cut_set_t*, Ntk::max_fanin_size + 1>;
  using wide_cut_merge_t = typename std::array<wide_cut_set_t*, Ntk::max_fanin_size + 1>;
  using TT = kitty::dynamic_truth_table;
  using tt_cache = truth_table_cache<TT>;
  using cost_cache = std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>;
  using lut_info = std::pair<kitty::dynamic_truth_table, std::vector<signal<klut_network>>>;
  using pack_index_list = std::array<uint32_t, max_gate_fanin>;

public:
  explicit lut_map_impl( Ntk& ntk, lut_map_params const& ps, lut_map_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        node_match( ntk.size() ),
        cuts( ntk.size() )
  {
    assert( ps.cut_enumeration_ps.cut_limit < max_cut_num && "cut_limit exceeds the compile-time limit for the maximum number of cuts" );

    if constexpr ( StoreFunction )
    {
      TT zero( 0u ), proj( 1u );
      kitty::create_nth_var( proj, 0u );

      tmp_visited.reserve( 100 );
      truth_tables.resize( 20000 );

      truth_tables.insert( zero );
      truth_tables.insert( proj );

      if constexpr ( !std::is_same<LUTCostFn, lut_unitary_cost>::value )
      {
        truth_tables_cost.reserve( 1000 );
      }
    }
  }

  klut_network run()
  {
    stopwatch t( st.time_total );

    /* compute and save topological order */
    topo_order.reserve( ntk.size() );
    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      topo_order.push_back( n );
    } );

    perform_mapping();
    return create_lut_network();
  }

  void run_inplace()
  {
    stopwatch t( st.time_total );

    if ( ps.multi_decomposition )
    {
      std::cerr << "[i] MAP ERROR: multi-input decompositions are not supported in \"inplace\" mapping\n" << std::endl;
      return;
    }

    /* compute and save topological order */
    topo_order.reserve( ntk.size() );
    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      topo_order.push_back( n );
    } );

    if ( ps.collapse_mffcs )
    {
      compute_mffcs_mapping();
      return;
    }

    perform_mapping();
    derive_mapping();
  }

private:
  void perform_mapping()
  {
    /* define area sorting function */
    lut_cut_sort_type area_sort = ( ps.area_oriented_mapping && !ps.edge_optimization ) ? lut_cut_sort_type::AREA : lut_cut_sort_type::AREA2;

    /* init the data structure */
    init_nodes();
    init_cuts();

    /* compute mapping for depth or area */
    if ( !ps.area_oriented_mapping )
    {
      compute_required_time();

      if ( ps.recompute_cuts )
      {
        compute_mapping<false, false>( lut_cut_sort_type::DELAY, true, true, true );
        compute_required_time();
        compute_mapping<false, false>( lut_cut_sort_type::DELAY2, true, true, false );
        compute_required_time();
        compute_mapping<true, false>( area_sort, true, true, false );
      }
      else
      {
        compute_mapping<false, false>( lut_cut_sort_type::DELAY2, true, true, true );
      }
    }
    else
    {
      compute_required_time();
      compute_mapping<true, false>( area_sort, false, true, true );
    }

    if ( ps.cut_expansion )
    {
      compute_required_time();
      expand_cuts<false>();
    }

    /* try backward area iterations */
    uint32_t i = 0;
    while ( i < ps.area_share_rounds )
    {
      compute_share_mapping( area_sort, i == 0 );

      if ( ps.cut_expansion )
      {
        expand_cuts<false>();
      }
      ++i;
    }

    /* compute mapping using global area flow */
    i = 0;
    while ( i < ps.area_flow_rounds )
    {
      compute_required_time();
      compute_mapping<true, false>( area_sort, false, ps.recompute_cuts, i == 0 );

      if ( ps.cut_expansion )
      {
        compute_required_time();
        expand_cuts<false>();
      }
      ++i;
    }

    /* compute mapping using exact area/edge */
    i = 0;
    while ( i < ps.ela_rounds )
    {
      compute_required_time();
      compute_mapping<true, true>( area_sort, false, ps.recompute_cuts, false );

      if ( ps.cut_expansion )
      {
        compute_required_time();
        expand_cuts<true>();
      }
      ++i;
    }
  }

  void init_nodes()
  {
    ntk.foreach_node( [this]( auto const& n ) {
      const auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      node_data.map_refs = ntk.fanout_size( n );
      node_data.est_refs = static_cast<float>( ntk.fanout_size( n ) );
    } );
  }

  void init_cuts()
  {
    /* init constant cut */
    add_zero_cut( ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) ), false );
    if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
      add_zero_cut( ntk.node_to_index( ntk.get_node( ntk.get_constant( true ) ) ), true );

    /* init PIs cuts */
    ntk.foreach_ci( [&]( auto const& n ) {
      add_unit_cut( ntk.node_to_index( n ) );
    } );

    if ( ps.multi_decomposition )
    {
      bins.reserve( max_wide_cut_size );
      bins_ext.reserve( max_wide_cut_size );
      packs.reserve( max_wide_cut_size );
      topo_order2.reserve( 2 * ntk.size() );
    }
  }

  template<bool DO_AREA, bool ELA>
  void compute_mapping( lut_cut_sort_type const sort, bool preprocess, bool recompute_cuts, bool decompose = false )
  {
    cuts_total = 0;
    for ( auto const& n : topo_order )
    {
      if constexpr ( !ELA )
      {
        auto const index = ntk.node_to_index( n );
        if ( !preprocess && iteration != 0 )
        {
          node_match[index].est_refs = ( 2.0 * node_match[index].est_refs + 1.0 * node_match[index].map_refs ) / 3.0;
        }
        else
        {
          node_match[index].est_refs = static_cast<float>( node_match[index].map_refs );
        }
      }

      if ( ntk.is_constant( n ) || ntk.is_ci( n ) )
      {
        if ( ps.multi_decomposition )
          topo_order2.push_back( n );
        continue;
      }

      if ( recompute_cuts )
      {
        if constexpr ( Ntk::min_fanin_size == 2 && Ntk::max_fanin_size == 2 )
        {
          compute_best_cut2<DO_AREA, ELA>( n, sort, preprocess );
        }
        else
        {
          if constexpr ( has_foreach_choice_v<Ntk> )
          {
            if ( ps.multi_decomposition && ntk.fanin_size( n ) > Ntk::min_fanin_size )
            {
              compute_best_cut_decompose<DO_AREA, ELA>( n, sort, preprocess, decompose );
            }
            else
            {
              compute_best_cut<DO_AREA, ELA>( n, sort, preprocess );
            }
          }
          else
          {
            compute_best_cut<DO_AREA, ELA>( n, sort, preprocess );
          }
        }
      }
      else
      {
        /* update cost the function and move the best one first */
        update_cut_data<DO_AREA, ELA>( n, sort );
      }

      if ( ps.multi_decomposition )
          topo_order2.push_back( n );
    }

    /* commit the new order with the additional nodes */
    if ( ps.multi_decomposition )
    {
      topo_order = topo_order2;
      topo_order2.clear();
    }

    set_mapping_refs<ELA>();

    if constexpr ( DO_AREA )
    {
      ++area_iteration;
    }

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats;

      if ( ( sort == lut_cut_sort_type::AREA || sort == lut_cut_sort_type::AREA2 ) && ELA )
      {
        stats << fmt::format( "[i] Area     : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
      }
      else if ( sort == lut_cut_sort_type::AREA || sort == lut_cut_sort_type::AREA2 )
      {
        stats << fmt::format( "[i] AreaFlow : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
      }
      else if ( sort == lut_cut_sort_type::DELAY2 )
      {
        stats << fmt::format( "[i] Delay2   : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
      }
      else
      {
        stats << fmt::format( "[i] Delay    : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
      }
      st.round_stats.push_back( stats.str() );
    }
  }

  void compute_share_mapping( lut_cut_sort_type const sort, bool first )
  {
    /* reset required times and references except for POs */
    compute_share_mapping_init( first );

    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      auto const index = ntk.node_to_index( *it );

      /* skip not used nodes */
      if ( node_match[index].map_refs == 0 )
        continue;

      /* update cost the function and move the best one first */
      update_cut_data_share( *it, sort );
    }

    /* propagate correct arrival times and compute stats */
    propagate_arrival_times();

    /* round stats */
    if ( ps.verbose )
    {
      st.round_stats.push_back( fmt::format( "[i] AreaSh   : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total ) );
    }
  }

  template<bool ELA>
  void expand_cuts()
  {
    /* cut expansion is not yet compatible with truth table computation */
    if constexpr ( StoreFunction )
      return;

    /* don't expand if cut recomputed cuts is off */
    if ( !ps.recompute_cuts )
      return;

    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_ci( n ) )
      {
        continue;
      }

      expand_cuts_node( n );
    }

    set_mapping_refs<ELA>();

    std::string stats = fmt::format( "[i] Reduce   : Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
    st.round_stats.push_back( stats );
  }

  template<bool ELA>
  void set_mapping_refs()
  {
    if constexpr ( !ELA )
    {
      for ( auto i = 0u; i < node_match.size(); ++i )
      {
        node_match[i].map_refs = 0u;
      }
    }

    /* compute the current worst delay and update the mapping refs */
    delay = 0;
    ntk.foreach_co( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );

      delay = std::max( delay, cuts[index][0]->data.delay );

      if constexpr ( !ELA )
      {
        ++node_match[index].map_refs;
      }
    } );

    /* compute current area and update mapping refs in top-down order */
    area = 0;
    edges = 0;
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) || ntk.is_ci( *it ) )
      {
        continue;
      }

      const auto index = ntk.node_to_index( *it );
      auto& node_data = node_match[index];

      /* continue if not referenced in the cover */
      if ( node_match[index].map_refs == 0u )
        continue;

      auto& best_cut = cuts[index][0];

      if constexpr ( !ELA )
      {
        for ( auto const leaf : best_cut )
        {
          node_match[leaf].map_refs++;
        }
      }
      area += best_cut->data.lut_area;
      edges += best_cut.size();
    }

    ++iteration;
  }

  void compute_required_time()
  {
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required = UINT32_MAX >> 1;
    }

    /* return in case of area_oriented_mapping */
    if ( iteration == 0 || ps.area_oriented_mapping )
      return;

    uint32_t required = delay;

    /* relax delay constraints */
    if ( ps.required_delay == 0.0f && ps.relax_required > 0.0f )
    {
      required *= ( 100.0 + ps.relax_required ) / 100.0;
    }

    if ( ps.required_delay != 0 )
    {
      /* Global target time constraint */
      if ( ps.required_delay < delay )
      {
        if ( !ps.area_oriented_mapping && iteration == 1 )
          std::cerr << fmt::format( "[i] MAP WARNING: cannot meet the target required time of {:.2f}", ps.required_delay ) << std::endl;
      }
      else
      {
        required = ps.required_delay;
      }
    }

    /* set the required time at POs */
    ntk.foreach_co( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      node_match[index].required = required;
    } );

    /* propagate required time to the PIs */
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if ( ntk.is_ci( *it ) || ntk.is_constant( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );

      if ( node_match[index].map_refs == 0 )
        continue;

      /* propagate required to choice nodes */
      if constexpr ( has_foreach_choice_v<Ntk> )
      {
        if ( ntk.has_choice( *it ) )
        {
          assert( ntk.is_choice_representative( *it ) );
          ntk.foreach_choice( *it, [&]( auto const& g ) {
            node_match[ntk.node_to_index( g )].required = node_match[index].required;
          } );
        }
      }

      for ( auto leaf : cuts[index][0] )
      {
        node_match[leaf].required = std::min( node_match[leaf].required, node_match[index].required - cuts[index][0]->data.lut_delay );
      }
    }
  }

  void propagate_arrival_times()
  {
    area = 0;
    edges = 0;
    for ( auto const& n : topo_order )
    {
      auto index = ntk.node_to_index( n );

      if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      {
        continue;
      }

      /* propagate arrival time */
      uint32_t node_delay = 0;
      auto& best_cut = **( cuts[index].begin() );

      for ( auto leaf : best_cut )
      {
        const auto& best_leaf_cut = cuts[leaf][0];
        node_delay = std::max( node_delay, best_leaf_cut->data.delay );
      }

      best_cut->data.delay = node_delay + best_cut->data.lut_delay;

      /* continue if not referenced in the cover */
      if ( node_match[index].map_refs == 0u )
        continue;

      /* update stats */
      area += best_cut->data.lut_area;
      edges += best_cut.size();
    }

    /* update worst delay */
    delay = 0;
    ntk.foreach_co( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, cuts[index][0]->data.delay );
    } );
  }

  void compute_share_mapping_init( bool first )
  {
    /* reset the mapping references and the required time */
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required = UINT32_MAX >> 1;
      if ( !first )
        node_match[i].est_refs = ( 2.0 * node_match[i].est_refs + 1.0 * node_match[i].map_refs ) / 3.0;
      else
        node_match[i].est_refs = std::max( 1.0, ( 1.0 * node_match[i].est_refs + 2.0 * node_match[i].map_refs ) / 3.0 );
      node_match[i].map_refs = 0;

      /* update flows if in area-oriented mapping */
      if ( ps.area_oriented_mapping )
        compute_cut_data<false>( **( cuts[i].begin() ), ntk.index_to_node( i ), false );
    }

    uint32_t required = delay;
    if ( ps.required_delay == 0.0f && ps.relax_required > 0.0f )
    {
      required *= ( 100.0 + ps.relax_required ) / 100.0;
    }

    if ( ps.required_delay != 0 )
    {
      /* Global target time constraint */
      if ( ps.required_delay < delay )
      {
        if ( !ps.area_oriented_mapping && iteration == 1 )
          std::cerr << fmt::format( "[i] MAP WARNING: cannot meet the target required time of {:.2f}", ps.required_delay ) << std::endl;
      }
      else
      {
        required = ps.required_delay;
      }
    }

    /* set the required time at POs */
    ntk.foreach_co( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      node_match[index].required = required;
      node_match[index].map_refs++;
    } );
  }

  template<bool DO_AREA, bool ELA>
  void compute_best_cut2( node const& n, lut_cut_sort_type const sort, bool preprocess, bool is_new = false )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    cut_t best_cut;

    /* compute cuts */
    const auto fanin = 2;
    uint32_t pairs{ 1 };
    ntk.foreach_fanin( ntk.index_to_node( index ), [this, &pairs]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
      pairs *= static_cast<uint32_t>( lcuts[i]->size() );
    } );
    lcuts[2] = &cuts[index];
    auto& rcuts = *lcuts[fanin];

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_deref( rcuts[0] );
      }
    }

    /* recompute the data of the best cut */
    if ( iteration != 0 && !is_new )
    {
      best_cut = rcuts[0];
      compute_cut_data<ELA>( best_cut, n, true );
    }

    /* clear cuts */
    rcuts.clear();

    /* insert the previous best cut */
    if ( iteration != 0 && !preprocess && !is_new )
    {
      rcuts.simple_insert( best_cut, sort );
    }

    cut_t new_cut;
    std::vector<cut_t const*> vcuts( fanin );

    for ( auto const& c1 : *lcuts[0] )
    {
      for ( auto const& c2 : *lcuts[1] )
      {
        if ( !c1->merge( *c2, new_cut, ps.cut_enumeration_ps.cut_size ) )
        {
          continue;
        }

        if ( ps.remove_dominated_cuts && rcuts.is_dominated( new_cut ) )
        {
          continue;
        }

        if constexpr ( StoreFunction )
        {
          vcuts[0] = c1;
          vcuts[1] = c2;
          new_cut->func_id = compute_truth_table( index, vcuts, new_cut );
        }

        compute_cut_data<ELA>( new_cut, ntk.index_to_node( index ), true );

        /* check required time */
        if constexpr ( DO_AREA )
        {
          if ( preprocess || new_cut->data.delay <= node_data.required )
          {
            if ( ps.remove_dominated_cuts )
              rcuts.insert( new_cut, false, sort );
            else
              rcuts.simple_insert( new_cut, sort );
          }
        }
        else
        {
          if ( ps.remove_dominated_cuts )
            rcuts.insert( new_cut, false, sort );
          else
            rcuts.simple_insert( new_cut, sort );
        }
      }
    }

    cuts_total += rcuts.size();

    /* limit the maximum number of cuts */
    rcuts.limit( ps.cut_enumeration_ps.cut_limit );

    /* replace the new best cut with previous one */
    if ( preprocess && !is_new && rcuts[0]->data.delay > node_data.required )
      rcuts.replace( 0, best_cut );

    /* add trivial cut */
    if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
    {
      add_unit_cut( index );
    }

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_ref( rcuts[0] );
      }
    }
  }

  template<bool DO_AREA, bool ELA>
  void compute_best_cut( node const& n, lut_cut_sort_type const sort, bool preprocess, bool is_new = false )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    cut_t best_cut;

    /* compute cuts */
    uint32_t pairs{ 1 };
    std::vector<uint32_t> cut_sizes;
    ntk.foreach_fanin( ntk.index_to_node( index ), [this, &pairs, &cut_sizes]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
      cut_sizes.push_back( static_cast<uint32_t>( lcuts[i]->size() ) );
      pairs *= cut_sizes.back();
    } );
    const auto fanin = cut_sizes.size();
    lcuts[fanin] = &cuts[index];
    auto& rcuts = *lcuts[fanin];

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_deref( rcuts[0] );
      }
    }

    /* recompute the data of the best cut */
    if ( iteration != 0 && !is_new )
    {
      best_cut = rcuts[0];
      compute_cut_data<ELA>( best_cut, n, true );
    }

    /* clear cuts */
    rcuts.clear();

    /* insert the previous best cut */
    if ( iteration != 0 && !preprocess && !is_new )
    {
      rcuts.simple_insert( best_cut, sort );
    }

    if ( fanin > 1 && fanin <= ps.cut_enumeration_ps.fanin_limit )
    {
      cut_t new_cut, tmp_cut;

      std::vector<cut_t const*> vcuts( fanin );

      foreach_mixed_radix_tuple( cut_sizes.begin(), cut_sizes.end(), [&]( auto begin, auto end ) {
        auto it = vcuts.begin();
        auto i = 0u;
        while ( begin != end )
        {
          *it++ = &( ( *lcuts[i++] )[*begin++] );
        }

        if ( !vcuts[0]->merge( *vcuts[1], new_cut, ps.cut_enumeration_ps.cut_size ) )
        {
          return true; /* continue */
        }

        for ( i = 2; i < fanin; ++i )
        {
          tmp_cut = new_cut;
          if ( !vcuts[i]->merge( tmp_cut, new_cut, ps.cut_enumeration_ps.cut_size ) )
          {
            return true; /* continue */
          }
        }

        if ( ps.remove_dominated_cuts && rcuts.is_dominated( new_cut ) )
        {
          return true; /* continue */
        }

        if constexpr ( StoreFunction )
        {
          new_cut->func_id = compute_truth_table( index, vcuts, new_cut );
        }

        compute_cut_data<ELA>( new_cut, index, true );

        /* check required time */
        if constexpr ( DO_AREA )
        {
          if ( preprocess || new_cut->data.delay <= node_data.required )
          {
            if ( ps.remove_dominated_cuts )
              rcuts.insert( new_cut, false, sort );
            else
              rcuts.simple_insert( new_cut, sort );
          }
        }
        else
        {
          if ( ps.remove_dominated_cuts )
            rcuts.insert( new_cut, false, sort );
          else
            rcuts.simple_insert( new_cut, sort );
        }

        return true;
      } );

      /* limit the maximum number of cuts */
      rcuts.limit( ps.cut_enumeration_ps.cut_limit );
    }
    else if ( fanin == 1 )
    {
      for ( auto const& cut : *lcuts[0] )
      {
        cut_t new_cut = *cut;

        if constexpr ( StoreFunction )
        {
          new_cut->func_id = compute_truth_table( index, { cut }, new_cut );
        }

        compute_cut_data<ELA>( new_cut, n, true );

        if constexpr ( DO_AREA )
        {
          if ( preprocess || new_cut->data.delay <= node_data.required )
          {
            if ( ps.remove_dominated_cuts )
              rcuts.insert( new_cut, false, sort );
            else
              rcuts.simple_insert( new_cut, sort );
          }
        }
        else
        {
          if ( ps.remove_dominated_cuts )
            rcuts.insert( new_cut, false, sort );
          else
            rcuts.simple_insert( new_cut, sort );
        }
      }

      /* limit the maximum number of cuts */
      rcuts.limit( ps.cut_enumeration_ps.cut_limit );
    }

    cuts_total += rcuts.size();

    /* replace the new best cut with previous one */
    if ( preprocess && !is_new && rcuts[0]->data.delay > node_data.required )
      rcuts.replace( 0, best_cut );

    add_unit_cut( index );

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_ref( rcuts[0] );
      }
    }
  }

  template<bool DO_AREA, bool ELA>
  void compute_best_cut_decompose( node const& n, lut_cut_sort_type const sort, bool preprocess, bool decompose )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto& rcuts = cuts[index];
    cut_t best_cut;

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_deref( rcuts[0] );
      }
    }

    /* recompute the data of the best cut */
    if ( iteration != 0 )
    {
      best_cut = rcuts[0];
      compute_cut_data<ELA>( best_cut, n, true );
    }

    /* clear cuts */
    rcuts.clear();

    /* insert the previous best wide cut */
    if ( iteration != 0 && !preprocess )
    {
      rcuts.simple_insert( best_cut, sort );
    }

    /* if is not a decomposition iteration, copy cuts from previous decompositions, and commit the best one */
    if ( !decompose )
    {
      commit_decomposition( n, rcuts, sort );

      if ( preprocess && rcuts[0]->data.delay > node_data.required )
      {
        rcuts.replace( 0, best_cut );
      }

      /* add trivial cut */
      if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
      {
        add_unit_cut( index );
      }

      if constexpr ( DO_AREA )
      {
        if ( iteration != 0 && node_data.map_refs > 0 )
        {
          cut_ref( rcuts[0] );
        }
      }
      return;
    }

    const auto fanin = ntk.fanin_size( n );
    uint32_t pairs{ 1 };
    ntk.foreach_fanin( ntk.index_to_node( index ), [this, &pairs]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
    } );
    lcuts[fanin] = &rcuts;

    /* load cut set of child 0 into rcuts */
    /* TODO: implement a better data structure with cut prioritization */
    wide_set.clear();

    adaptive_buckets_t *source = &dcuts1;
    adaptive_buckets_t *dest = &dcuts2;
    init_dcuts();

    for ( uint32_t i = 1; i < fanin; ++i )
    {
      merge_cuts2_subset<DO_AREA, ELA>( n, source, dest, i, sort, preprocess );
      std::swap( source, dest );
    }

    /* create official cuts */
    wide_cut_t new_cut;
    uint32_t limit = 0;
    for ( uint32_t i = fanin - 1; i < max_wide_cut_size; ++i )
    {
      for ( auto const& dcut : source->at( i ) )
      {
        dcut_merge( n, *dcut, new_cut, fanin );
        compute_wide_cut_data<ELA>( new_cut, ntk.index_to_node( index ) );

        // assert( dcut->data.lut_area >= new_cut->data.lut_area );
        // assert( dcut->data.delay >= new_cut->data.delay );
        // assert( dcut->data.area_flow >= new_cut->data.area_flow - epsilon );
        // assert( dcut->data.edge_flow >= new_cut->data.edge_flow - epsilon);

        /* check required time */
        if constexpr ( DO_AREA )
        {
          if ( preprocess || new_cut->data.delay <= node_data.required )
          {
              wide_set.simple_insert( new_cut, sort );
          }
        }
        else
        {
            wide_set.simple_insert( new_cut, sort );
        }
      }
    }

    /* decompose the new wide cut */
    if ( wide_set.size() && ( !preprocess || wide_set[0]->data.delay <= node_data.required ) )
    {
      /* add AND nodes to represent the decomposition */
      decompose_wide_cut<DO_AREA, ELA>( n, wide_set[0], sort, preprocess );
    }

    commit_decomposition( n, rcuts, sort );

    if ( preprocess && rcuts[0]->data.delay > node_data.required )
    {
      rcuts.replace( 0, best_cut );
    }

    /* add trivial cut */
    if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
    {
      add_unit_cut( index );
    }

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_ref( rcuts[0] );
      }
    }
  }


  void commit_decomposition( node const& n, cut_set_t& rcuts, lut_cut_sort_type const sort )
  {
    /* merge the choice cuts */
    ntk.foreach_choice( n, [&]( auto const& g ) {
      if ( g == n )
        return;

      for ( auto const& cut : cuts[ntk.node_to_index( g )] )
      {
        if ( cut->size() == 1 )
          continue;
        if ( ps.remove_dominated_cuts )
          rcuts.insert( *cut, false, sort );
        else
          rcuts.simple_insert( *cut, sort );
      }
    } );

    /* limit the maximum number of cuts */
    rcuts.limit( ps.cut_enumeration_ps.cut_limit );
  }

  void init_dcuts()
  {
    uint32_t i;

    /* clear the buckets from cuts */
    for ( i = 0; i < max_wide_cut_size; ++i )
      dcuts1[i].clear();

    /* insert cuts */
    i = 0;
    for ( auto const& c : *lcuts[0] )
    {
      adaptive_cut_t acut;
      acut.cut_pointers[0] = i++;
      acut.data = ( *c )->data;
      acut.size = c->size();
      dcuts1[c->size() - 1].insert( acut );
    }
  }

  template<bool DO_AREA, bool ELA>
  void merge_cuts2_subset( node const& n, adaptive_buckets_t *source, adaptive_buckets_t *dest, uint32_t child, lut_cut_sort_type const sort, bool preprocess )
  {
    uint32_t index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    adaptive_cut_t new_cut;

    /* clear the dest buckets from cuts */
    for ( uint32_t i = 0; i < max_wide_cut_size; ++i )
      dest->at( i ).clear();

    
    for ( uint32_t i = 0; i < max_wide_cut_size; ++i )
    {
      for ( auto const& dc1 : source->at( i ) )
      {
        uint32_t j = 0;
        for ( auto const& c2 : *lcuts[child] )
        {
          /* skip large cuts exceeding max_wide_cut_size */
          if ( dc1->size + c2->size() > max_wide_cut_size )
          {
            ++j;
            continue;
          }

          /* create the cut */
          new_cut = *dc1;
          new_cut.cut_pointers[child] = j++;
          compute_packing_data<ELA>( n, new_cut, child + 1 );

          dest->at( new_cut.size - 1 ).insert( new_cut );
        }
      }
    }
  }

  template<bool ELA>
  void compute_packing_data( node const& n, adaptive_cut_t& new_cut, uint32_t leaves )
  {
    uint32_t size = 0;
    float area_flow = 0;
    float edge_flow = 0;
    bins_ext.clear();

    /* perform 0-level bin-packing */
    for ( uint32_t i = 0; i < leaves; ++i )
    {
      cut_t const& cut = ( *lcuts[i] )[new_cut.cut_pointers[i]];
      uint32_t cut_size = cut.size();

      uint32_t j = 0;
      while ( j < bins_ext.size() )
      {
        if ( bins_ext[j].size + cut_size <= ps.cut_enumeration_ps.cut_size )
          break;
        ++j;
      }

      if ( j == bins_ext.size() )
        bins_ext.emplace_back();

      uint32_t cut_delay = cut->data.delay;
      float cut_area_flow = cut->data.area_flow - cut->data.lut_area;
      float cut_edge_flow = cut->data.edge_flow;

      /* compute cost if trivial cut */
      if ( cut.size() == 1 )
      {
        cut_delay = ( *lcuts[i] )[0]->data.delay + 1;
        
        if constexpr ( ELA )
        {
          cut_area_flow = ( *lcuts[i] )[0]->data.area_flow;
          cut_edge_flow = ( *lcuts[i] )[0]->data.edge_flow + 1;
        }
        else
        {
          uint32_t leaf = *( cut.begin() );
          if ( node_match[leaf].map_refs > 0 )
          {
            cut_area_flow = ( *lcuts[i] )[0]->data.area_flow / node_match[leaf].est_refs;
            cut_edge_flow = ( *lcuts[i] )[0]->data.edge_flow / node_match[leaf].est_refs + 1;
          }
          else
          {
            cut_area_flow = ( *lcuts[i] )[0]->data.area_flow;
            cut_edge_flow = ( *lcuts[i] )[0]->data.edge_flow + 1;
          }
        }
      }

      size += cut_size;
      bins_ext[j].size += cut_size;
      bins_ext[j].delay = std::max( bins_ext[j].delay, cut_delay );
      area_flow += cut_area_flow;
      edge_flow += cut_edge_flow;
    }

    uint32_t area = bins_ext.size();
    uint32_t roots = bins_ext.size();

    // /* perform multi-level bin-packing */
    // std::sort( bins_ext.begin(), bins_ext.end(), []( auto const& a, auto const& b ) {
    //   return a.size > b.size;
    // } );

    // /* check fully-packable */
    // uint32_t max_available = ps.cut_enumeration_ps.cut_size - bins_ext.back().size;
    // if ( roots != 1 && roots - 1 <= max_available )
    // {
    //   roots = 1;
    //   /* update delay */
    //   for ( uint32_t i = 0; i < bins_ext.size() - 1; ++i )
    //     bins_ext.back().delay = std::max( bins_ext.back().delay, bins_ext[i].delay + 1 );
    // }

    // /* check if it cannot be packed but a top LUT can map the network */
    // bool skip_multi_level_packing = false;
    // if ( roots != 1 && roots <= ps.cut_enumeration_ps.cut_size )
    // {
    //   /* check availability */
    //   uint32_t availability = 0;
    //   for ( uint32_t i = 1; i < bins_ext.size(); ++i )
    //   {
    //     availability += ps.cut_enumeration_ps.cut_size - bins_ext[i].size;
    //   }
    //   if ( availability < roots - 1 )
    //     skip_multi_level_packing = true;
    // }

    // uint32_t i;
    // uint32_t unconnected = 1u;
    // for ( i = 1; i < bins_ext.size() && roots != 1 && !skip_multi_level_packing; ++i )
    // {
    //   if ( bins_ext[i].size < ps.cut_enumeration_ps.cut_size )
    //   {
    //     uint32_t connect = std::min( ps.cut_enumeration_ps.cut_size - bins_ext[i].size, unconnected );
    //     bins_ext[i].size += connect;
    //     unconnected -= connect;
    //     roots -= connect;

    //     /* update delay */
    //     for ( uint32_t j = 0; j < i && connect > 0; ++j )
    //     {
    //       if ( bins_ext[j].connected )
    //         continue;
          
    //       bins_ext[j].connected = true;
    //       bins_ext[i].delay = std::max( bins_ext[i].delay, bins_ext[j].delay + 1 );
    //       --connect;
    //     }

    //     /* check terminal condition minimizing delay */
    //     max_available = ps.cut_enumeration_ps.cut_size - bins.back();
    //     if ( roots != 1 && roots - 1 <= max_available )
    //     {
    //       roots = 1;
    //       /* update delay */
    //       for ( uint32_t j = 0; j < bins_ext.size() - 1; ++j )
    //       {
    //         if ( bins_ext[j].connected )
    //           continue;
    //         bins_ext.back().delay = std::max( bins_ext.back().delay, bins_ext[j].delay + 1 );
    //       }
    //       break;
    //     }
    //   }
    //   else
    //   {
    //     ++unconnected;
    //   }
    // }

    // /* add top LUTs until roots are 1 */
    // i = 0;
    // while( roots != 1 )
    // {
    //   bins_ext.emplace_back();
    //   ++area;
    //   ++roots;

    //   /* find unconnected LUTs */
    //   uint64_t lut_fanin = 0;
    //   for ( ; i < bins_ext.size() - 1 && bins_ext.back().size <= ps.cut_enumeration_ps.cut_size; ++i )
    //   {
    //     if ( bins_ext[i].connected )
    //       continue;
    //     bins_ext[i].connected = true;
    //     bins_ext.back().size++;
    //     bins_ext.back().delay = std::max( bins_ext.back().delay, bins_ext[i].delay + 1 );
    //     --roots;
    //   }
    // }

    /* write  implementation cost */
    new_cut.data.lut_delay = 1;
    new_cut.data.lut_area = area;
    new_cut.data.delay = bins_ext.back().delay;
    new_cut.data.area_flow = area_flow + area;
    new_cut.data.edge_flow = edge_flow + area - 1;
    new_cut.size = size;
  }

  void dcut_merge( node const& n, adaptive_cut_t const& cut, wide_cut_t& new_cut, uint32_t num_leaves )
  {
    /* clear new_cut fields */
    new_cut->data.lut_roots.clear();
    new_cut->data.lut_roots.reserve( 2 * max_gate_fanin );
    new_cut->data.delays.clear();
    new_cut->data.delays.reserve( 2 * max_gate_fanin );

    /* perform bin-packing */
    bins.clear();
    packs.clear();
    lbins.clear();

    /* save indexes for variable swapping */
    std::iota( pack_indexes.begin(), pack_indexes.begin() + num_leaves, 0 );

    /* sort cut sizes in decreasing order */
    for ( uint32_t i = 0; i < num_leaves; ++i )
    {
      packs.push_back( ( *lcuts[i] )[cut.cut_pointers[i]].size() );
    }
    std::sort( pack_indexes.begin(), pack_indexes.begin() + num_leaves, [this]( auto const& a, auto const& b ) {
      return packs[a] > packs[b];
    } );

    /* create cuts for bins */
    // std::vector<cut_t> bins_cuts;

    /* TODO: verify correctness */
    /* perform 0-level bin-packing: use the LUT delay field to annotate which nodes are merged */
    cut_t merge_cut;
    cut_t bin_cut;
    for ( uint32_t i = 0; i < num_leaves; ++i )
    {
      cut_t leaf_cut = ( *lcuts[pack_indexes[i]] )[cut.cut_pointers[pack_indexes[i]]];

      uint32_t j = 0;
      while ( j < lbins.size )
      {
        bin_cut.set_leaves( lbins[j] );
        if ( bin_cut.merge( leaf_cut, merge_cut, ps.cut_enumeration_ps.cut_size ) )
        {
          lbins[j].set_leaves( merge_cut.begin(), merge_cut.end() );
          /* TODO: remove */
          lbins[j]->data.connections |= static_cast<uint64_t>( 1 ) << pack_indexes[i];
          lbins[j]->data.delay = std::max( lbins[j]->data.delay, get_cut_delay_lcuts( leaf_cut, pack_indexes[i] ) );
          break;
        }
        ++j;
      }

      if ( j == lbins.size )
      {
        assert( leaf_cut.size() <= ps.cut_enumeration_ps.cut_size );
        lbins[j].set_leaves( leaf_cut );
        lbins[j]->data.delay = get_cut_delay_lcuts( leaf_cut, pack_indexes[i] );
        lbins[j]->data.connections = static_cast<uint64_t>( 1 ) << pack_indexes[i];
        lbins[j]->data.lut_area = 0;
        lbins.size++;
      }
    }

    uint32_t area = lbins.size;
    uint32_t roots = lbins.size;
    uint32_t edges = 0;

    /* sort cuts based on size */
    std::sort( lbins.pcuts.begin(), lbins.pcuts.begin() + lbins.size, []( auto const& a, auto const& b ) {
      return a->size() > b->size();
    } );

    /* load bins and compute edges */
    for ( uint32_t i = 0; i < lbins.size; ++i )
    {
      bins.push_back( lbins[i].size() );
      edges += lbins[i].size();
    }

    /* check fully-packable using the last LUT */
    uint32_t max_available = ps.cut_enumeration_ps.cut_size - bins.back();
    if ( roots != 1 && roots - 1 <= max_available )
    {
      bins.back() += roots - 1;
      roots = 1;
      /* connect LUTs */
      for ( uint32_t i = 0; i < lbins.size - 1; ++i )
      {
        lbins[lbins.size - 1]->data.connections |= static_cast<uint64_t>( 1 ) << ( 32 + i );
        lbins[lbins.size - 1]->data.delay = std::max( lbins[lbins.size - 1]->data.delay, lbins[i]->data.delay + 1 );
      }
    }

    /* check if it cannot be packed but a top LUT can map the network */
    bool skip_multi_level_packing = false;
    if ( roots != 1 && roots <= ps.cut_enumeration_ps.cut_size )
    {
      /* check availability */
      uint32_t availability = 0;
      for ( uint32_t i = 1; i < bins.size(); ++i )
      {
        availability += ps.cut_enumeration_ps.cut_size - bins[i];
      }
      if ( availability < roots - 1 )
        skip_multi_level_packing = true;
    }

    uint32_t i;
    uint32_t unconnected = 1;
    std::array<bool, 2 * max_gate_fanin> connected = { false };
    for ( i = 1; i < lbins.size && roots != 1 && !skip_multi_level_packing; ++i )
    {
      uint32_t cut_size = bins[i];
      if ( cut_size >= ps.cut_enumeration_ps.cut_size )
      {
        ++unconnected;
        continue;
      }

      uint32_t connect = std::min( ps.cut_enumeration_ps.cut_size - cut_size, unconnected );
      bins[i] += connect;
      unconnected -= connect;
      roots -= connect;

      for ( uint32_t j = 0; j < i && connect != 0; ++j )
      {
        if ( connected[j] )
          continue;
        connected[j] = true;
        lbins[i]->data.connections |= static_cast<uint64_t>( 1 ) << ( 32 + j );
        lbins[i]->data.delay = std::max( lbins[i]->data.delay, lbins[j]->data.delay + 1 );
        --connect;
      }

      /* check terminal condition minimizing delay */
      max_available = ps.cut_enumeration_ps.cut_size - bins.back();
      if ( roots != 1 && roots - 1 <= max_available )
      {
        bins.back() += roots - 1;
        roots = 1;

        for ( uint32_t j = 0; j < lbins.size - 1; ++j )
        {
          if ( connected[j] )
            continue;
          lbins[lbins.size - 1]->data.connections |= static_cast<uint64_t>( 1 ) << ( 32 + j );
          lbins[lbins.size - 1]->data.delay = std::max( lbins[lbins.size - 1]->data.delay, lbins[j]->data.delay + 1 );
        }
        break;
      }
    }

    /* create cut */
    uint32_t size = 0;
    std::array<uint32_t, max_wide_cut_size> leaves;
    for ( i = 0; i < lbins.size; ++i )
    {
      auto const& leaf_cut = lbins[i];
      for ( auto leaf : leaf_cut )
      {
        leaves[size++] = leaf;
        ntk.set_value( ntk.node_to_index( leaf ), 0 );
      }

      /* add LUT connection info */
      new_cut->data.lut_roots.push_back( leaf_cut->data.connections );
      new_cut->data.delays.push_back( leaf_cut->data.delay );
    }

    new_cut.set_leaves( leaves.begin(), leaves.begin() + size );

    /* add top LUTs until roots are 1 */
    uint32_t k = 0;
    while( roots != 1 )
    {
      bins.push_back( 0 );
      ++area;
      ++roots;
      /* find unconnected LUTs */
      uint64_t lut_fanin = 0;
      uint32_t worst_delay = 0;
      for ( ; k < bins.size() - 1 && bins.back() <= ps.cut_enumeration_ps.cut_size; ++k )
      {
        if ( connected[k] )
          continue;
        connected[k] = true;
        lut_fanin |= static_cast<uint64_t>( 1 ) << k;
        worst_delay = std::max( worst_delay, new_cut->data.delays[k] + 1 );
        ++bins.back();
        --roots;
      }
      new_cut->data.lut_roots.push_back( lut_fanin << 32 );
      new_cut->data.delays.push_back( worst_delay );
    }

    /* write  implementation cost */
    new_cut->data.delay = new_cut->data.delays.back();
    new_cut->data.lut_delay = 1;
    new_cut->data.lut_area = area;
    new_cut->data.lut_edges = edges + area - 1;
  }

  template<bool DO_AREA, bool ELA>
  void decompose_wide_cut( node const& n, wide_cut_t const& cut, lut_cut_sort_type const sort, bool preprocess )
  {
    std::array<signal<Ntk>, max_wide_cut_size> lut_signals;
    std::vector<signal<Ntk>> children;
    children.reserve( ps.cut_enumeration_ps.cut_size );
    uint32_t ntk_size = ntk.size();

    /* for each LUT create a set of ANDs */
    for ( uint32_t i = 0; i < cut->data.lut_roots.size(); ++i )
    {
      children.clear();
      uint64_t lut_connections = cut->data.lut_roots[i];

      /* TODO: try creating a multi-input AND here */
      /* get signal pairs to connect */
      if ( lut_connections & UINT32_MAX )
      {
        for ( uint32_t j = 0; j < ntk.fanin_size( n ); ++j )
        {
          if ( ( lut_connections >> j ) & 1 )
            children.push_back( ntk._storage->nodes[n].children[j] );
        }

        create_and2_decomposition( children );
      }

      /* check LUTs connections */
      uint32_t in_luts = static_cast<uint32_t>( lut_connections >> 32 );
      if ( in_luts )
      {
        for ( uint32_t j = 0; j < i; ++j )
        {
          if ( ( in_luts >> j ) & 1 )
            children.push_back( lut_signals[j] );
        }

        create_and2_decomposition( children );
      }

      assert( children.size() == 1 );
      lut_signals[i] = children[0];
    }

    /* add new node as a choice */
    ntk.add_choice( n, children[0] );
    /* compute cuts recursively, add best cuts of the root to n, expand the data structures */
    /* wide cuts should only contain an history of selected wide cuts to avoid to reallocate nodes */
    /* alternatively, wide cuts representation in the data structure can be suppressed */
    /* add new nodes to the topo order */
    for ( uint32_t i = ntk_size; i < ntk.size(); ++i )
    {
      /* create cut set, add it to the data structure, compute cuts */
      cuts.emplace_back();
      node_match.emplace_back();
      auto& node_data = node_match[i];
      node const g = ntk.index_to_node( i );

      node_data.map_refs = 0;
      if ( i == ntk.size() - 1 )
        node_data.est_refs = node_match[ntk.node_to_index( n )].est_refs;
      else
        node_data.est_refs = 1.0;

      /* TODO: correct the required time */
      node_data.required = node_match[ntk.node_to_index( n )].required;

      compute_best_cut2<DO_AREA, ELA>( g, sort, preprocess, true );
      topo_order2.push_back( g );
    }
  }

  void create_and2_decomposition( std::vector<signal<Ntk>>& children )
  {
    while ( children.size() > 1 )
    {
      /* pick the last 2 signals and create an AND */
      auto const f0 = children.back();
      children.pop_back();
      auto const f1 = children.back();
      children.pop_back();
      children.push_back( ntk.create_and( f0, f1 ) );
    }
  }

  template<bool DO_AREA, bool ELA>
  void update_cut_data( node const& n, lut_cut_sort_type const sort )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto& node_cut_set = cuts[index];
    uint32_t best_cut_index = 0;
    uint32_t cut_index = 0;

    cut_t const* best_cut = &node_cut_set.best();

    if constexpr ( DO_AREA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_deref( *best_cut );
      }
    }

    /* recompute the data for all the cuts and pick the best */
    for ( cut_t* cut : node_cut_set )
    {
      /* skip trivial cut */
      if ( cut->size() == 1 && *cut->begin() == index )
      {
        ++cut_index;
        continue;
      }

      compute_cut_data<ELA>( *cut, n, false );

      /* update best */
      if constexpr ( DO_AREA )
      {
        if ( ( *cut )->data.delay <= node_data.required )
        {
          if ( node_cut_set.compare( *cut, *best_cut, sort ) )
          {
            best_cut = cut;
            best_cut_index = cut_index;
          }
        }
      }
      else
      {
        if ( node_cut_set.compare( *cut, *best_cut, sort ) )
        {
          best_cut = cut;
          best_cut_index = cut_index;
        }
      }

      ++cut_index;
    }

    if constexpr ( DO_AREA || ELA )
    {
      if ( iteration != 0 && node_data.map_refs > 0 )
      {
        cut_ref( *best_cut );
      }
    }

    /* update the best cut */
    node_cut_set.update_best( best_cut_index );
  }

  void update_cut_data_share( node const& n, lut_cut_sort_type const sort )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto& node_cut_set = cuts[index];
    uint32_t best_cut_index = 0;
    uint32_t cut_index = 0;

    cut_t const* best_cut = &node_cut_set.best();

    /* recompute the data for all the cuts and pick the best */
    for ( cut_t* cut : node_cut_set )
    {
      /* skip trivial cut */
      if ( cut->size() == 1 && *cut->begin() == index )
      {
        ++cut_index;
        continue;
      }

      compute_cut_data_share( *cut );

      /* update best */
      if ( ( *cut )->data.delay <= node_data.required )
      {
        if ( node_cut_set.compare( *cut, *best_cut, sort ) )
        {
          best_cut = cut;
          best_cut_index = cut_index;
        }
      }

      ++cut_index;
    }

    /* propagate required times backward and reference the leaves */
    for ( auto leaf : *best_cut )
    {
      node_match[leaf].required = std::min( node_match[leaf].required, node_data.required - ( *best_cut )->data.lut_delay );
      node_match[leaf].map_refs++;
    }

    /* update the best cut */
    node_cut_set.update_best( best_cut_index );
  }

  // template<bool DO_AREA, bool ELA>
  // void update_wide_cut_data( node const& n, lut_cut_sort_type const sort )
  // {
  //   auto index = ntk.node_to_index( n );
  //   auto& node_data = node_match[index];
  //   auto& node_cut_set = wide_set_cuts;
  //   uint32_t best_cut_index = 0;
  //   uint32_t cut_index = 0;

  //   wide_cut_t const* best_cut = &node_cut_set.best();

  //   if constexpr ( DO_AREA )
  //   {
  //     if ( iteration != 0 && node_data.map_refs > 0 )
  //     {
  //       cut_deref<wide_cut_t>( *best_cut );
  //     }
  //   }

  //   /* recompute the data for all the cuts and pick the best */
  //   for ( wide_cut_t* cut : node_cut_set )
  //   {
  //     /* skip trivial cut */
  //     if ( cut->size() == 1 && *cut->begin() == index )
  //     {
  //       ++cut_index;
  //       continue;
  //     }

  //     compute_wide_cut_data<ELA>( *cut, n );

  //     /* update best */
  //     if constexpr ( DO_AREA )
  //     {
  //       if ( ( *cut )->data.delay <= node_data.required )
  //       {
  //         if ( node_cut_set.compare( *cut, *best_cut, sort ) )
  //         {
  //           best_cut = cut;
  //           best_cut_index = cut_index;
  //         }
  //       }
  //     }
  //     else
  //     {
  //       if ( node_cut_set.compare( *cut, *best_cut, sort ) )
  //       {
  //         best_cut = cut;
  //         best_cut_index = cut_index;
  //       }
  //     }

  //     ++cut_index;
  //   }

  //   if constexpr ( DO_AREA || ELA )
  //   {
  //     if ( iteration != 0 && node_data.map_refs > 0 )
  //     {
  //       cut_ref<wide_cut_t>( *best_cut );
  //     }
  //   }

  //   /* update the best cut */
  //   node_cut_set.update_best( best_cut_index );
  //   cuts[index].update_best( best_cut_index );
  // }

  void expand_cuts_node( node const& n )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    cut_t best_cut = cuts[index][0];

    if ( node_data.map_refs == 0 )
      return;

    /* update delay */
    uint32_t delay_update = 0;
    for ( auto const leaf : best_cut )
    {
      delay_update = std::max( delay_update, cuts[leaf][0]->data.delay + best_cut->data.lut_delay );
    }
    best_cut->data.delay = delay_update;

    auto const area_before = cut_deref( best_cut );

    uint32_t cost_before = 0;

    std::vector<uint32_t> leaves;

    /* mark volume */
    ntk.incr_trav_id();
    for ( auto const leaf : best_cut )
    {
      ntk.set_visited( ntk.index_to_node( leaf ), ntk.trav_id() );
      leaves.push_back( leaf );

      /* MFFC leaves */
      if ( node_match[leaf].map_refs == 0 )
        ++cost_before;
    }
    mark_cut_volume_rec( n );

    /* improve cut */
    while ( improve_cut( leaves ) )
      ;

    /* measure improvement */
    uint32_t cost_after = 0;
    for ( auto const leaf : leaves )
    {
      /* MFFC leaves */
      if ( node_match[leaf].map_refs == 0 )
        ++cost_after;
    }

    assert( cost_after <= cost_before );

    /* create the new cut */
    cut_t new_cut;
    new_cut.set_leaves( leaves.begin(), leaves.end() );
    new_cut->data = best_cut->data;

    uint32_t delay_after = 0;
    for ( auto const leaf : leaves )
    {
      delay_after = std::max( delay_after, cuts[leaf][0]->data.delay + new_cut->data.lut_delay );
    }
    new_cut->data.delay = delay_after;

    auto const area_after = cut_ref( new_cut );

    /* new cut is better */
    if ( area_after <= area_before && new_cut->data.delay <= node_data.required )
    {
      cuts[index].replace( 0, new_cut );
    }
    else
    {
      /* restore */
      cut_deref( new_cut );
      cut_ref( best_cut );
    }
  }

  bool improve_cut( std::vector<uint32_t>& leaves )
  {
    if ( improve_cut_expand0( leaves ) )
      return true;

    if ( leaves.size() < ps.cut_enumeration_ps.cut_size && improve_cut_expand1( leaves ) )
      return true;

    assert( leaves.size() <= ps.cut_enumeration_ps.cut_size );
    return false;
  }

  bool improve_cut_expand0( std::vector<uint32_t>& leaves )
  {
    for ( auto it = leaves.begin(); it != leaves.end(); ++it )
    {
      if ( ntk.is_ci( *it ) )
        continue;

      /* test if expansion would increase the number of leaves */
      int marked = 0;
      ntk.foreach_fanin( ntk.index_to_node( *it ), [&]( auto const& f ) {
        if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() )
          ++marked;
      } );

      if ( marked > 1 )
        continue;

      /* check that the cost does not increase */
      marked = 0;
      if ( node_match[*it].map_refs == 0 )
        --marked;

      ntk.foreach_fanin( ntk.index_to_node( *it ), [&]( auto const& f ) {
        if ( ntk.is_constant( ntk.get_node( f ) ) )
          return;
        auto const index = ntk.node_to_index( ntk.get_node( f ) );
        if ( ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() && node_match[index].map_refs == 0 )
          ++marked;
      } );

      /* not referenced leaves don't increase from the transformation */
      if ( marked <= 0 )
      {
        /* update leaves */
        uint32_t n = *it;
        leaves.erase( it );
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          auto const index = ntk.node_to_index( ntk.get_node( f ) );
          if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() )
          {
            leaves.push_back( index );
            ntk.set_visited( ntk.get_node( f ), ntk.trav_id() );
          }
        } );
        return true;
      }
    }

    return false;
  }

  bool improve_cut_expand1( std::vector<uint32_t>& leaves )
  {
    for ( auto it = leaves.begin(); it != leaves.end(); ++it )
    {
      if ( ntk.is_ci( *it ) )
        continue;

      /* test if expansion would increase the number of leaves by more than 1*/
      int marked = 0;
      ntk.foreach_fanin( ntk.index_to_node( *it ), [&]( auto const& f ) {
        if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() )
          ++marked;
      } );

      if ( marked > 2 )
        continue;

      /* check that the cost reduces */
      marked = 0;
      if ( node_match[*it].map_refs == 0 )
        --marked;

      ntk.foreach_fanin( ntk.index_to_node( *it ), [&]( auto const& f ) {
        if ( ntk.is_constant( ntk.get_node( f ) ) )
          return;
        auto const index = ntk.node_to_index( ntk.get_node( f ) );
        if ( ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() && node_match[index].map_refs == 0 )
          ++marked;
      } );

      /* not referenced leaves should be reduced by the transformation */
      if ( marked < 0 )
      {
        /* update leaves */
        uint32_t n = *it;
        leaves.erase( it );
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          auto const index = ntk.node_to_index( ntk.get_node( f ) );
          if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() )
          {
            leaves.push_back( index );
            ntk.set_visited( ntk.get_node( f ), ntk.trav_id() );
          }
        } );
        return true;
      }
    }

    return false;
  }

  template<typename CutType = cut_t>
  uint32_t cut_ref( CutType const& cut )
  {
    uint32_t count = cut->data.lut_area;

    for ( auto leaf : cut )
    {
      if ( ntk.is_ci( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( node_match[leaf].map_refs++ == 0u )
      {
        count += cut_ref<cut_t>( cuts[leaf][0] );
      }
    }

    return count;
  }

  template<typename CutType = cut_t>
  uint32_t cut_deref( CutType const& cut )
  {
    uint32_t count = cut->data.lut_area;

    for ( auto leaf : cut )
    {
      if ( ntk.is_ci( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( --node_match[leaf].map_refs == 0u )
      {
        count += cut_deref<cut_t>( cuts[leaf][0] );
      }
    }

    return count;
  }

  template<typename CutType = cut_t>
  uint32_t cut_measure_mffc( CutType const& cut )
  {
    tmp_visited.clear();

    uint32_t count = cut_ref_visit<CutType>( cut );

    /* dereference visited */
    for ( auto const& s : tmp_visited )
    {
      --node_match[s].map_refs;
    }

    return count;
  }

  template<typename CutType = cut_t>
  uint32_t cut_ref_visit( CutType const& cut )
  {
    uint32_t count = cut->data.lut_area;

    for ( auto leaf : cut )
    {
      if ( ntk.is_ci( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* add to visited */
      tmp_visited.push_back( leaf );

      /* Recursive referencing if leaf was not referenced */
      if ( node_match[leaf].map_refs++ == 0u )
      {
        count += cut_ref_visit<cut_t>( cuts[leaf][0] );
      }
    }

    return count;
  }

  uint32_t cut_edge_ref( cut_t const& cut )
  {
    uint32_t count = cut.size();

    for ( auto leaf : cut )
    {
      if ( ntk.is_ci( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( node_match[leaf].map_refs++ == 0u )
      {
        count += cut_edge_ref( cuts[leaf][0] );
      }
    }
    return count;
  }

  template<typename CutType = cut_t>
  uint32_t cut_edge_deref( CutType const& cut )
  {
    uint32_t count;
    if constexpr ( std::is_same<CutType, wide_cut_t>::value )
      count = cut->data.lut_edges;
    else
      count = cut.size();

    for ( auto leaf : cut )
    {
      if ( ntk.is_ci( ntk.index_to_node( leaf ) ) || ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* Recursive referencing if leaf was not referenced */
      if ( --node_match[leaf].map_refs == 0u )
      {
        count += cut_edge_deref<cut_t>( cuts[leaf][0] );
      }
    }
    return count;
  }

  void mark_cut_volume_rec( node const& n )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      mark_cut_volume_rec( ntk.get_node( f ) );
    } );
  }

  /* compute positions of leave indices in cut `sub` (subset) with respect to
   * leaves in cut `sup` (super set).
   *
   * Example:
   *   compute_truth_table_support( {1, 3, 6}, {0, 1, 2, 3, 6, 7} ) = {1, 3, 4}
   */
  void compute_truth_table_support( cut_t const& sub, cut_t const& sup, TT& tt )
  {
    std::vector<uint8_t> support( sub.size() );

    size_t j = 0;
    auto itp = sup.begin();
    for ( auto i : sub )
    {
      itp = std::find( itp, sup.end(), i );
      support[j++] = static_cast<uint8_t>( std::distance( sup.begin(), itp ) );
    }

    /* swap variables in the truth table */
    for ( int i = j - 1; i >= 0; --i )
    {
      assert( i <= support[i] );
      kitty::swap_inplace( tt, i, support[i] );
    }
  }

  template<bool ELA>
  void compute_cut_data( cut_t& cut, node const& n, bool recompute_cut_cost )
  {
    uint32_t lut_area;
    uint32_t lut_delay;

    if ( recompute_cut_cost )
    {
      if constexpr ( StoreFunction )
      {
        if constexpr ( !std::is_same<LUTCostFn, lut_unitary_cost>::value )
        {
          if ( auto it = truth_tables_cost.find( cut->func_id ); it != truth_tables_cost.end() )
          {
            std::tie( lut_area, lut_delay ) = it->second;
          }
          else
          {
            auto cost = lut_cost( truth_tables[cut->func_id] );
            if ( truth_tables[cut->func_id].num_vars() <= ps.cost_cache_vars )
            {
              /* cache it */
              truth_tables_cost[cut->func_id] = cost;
            }
            lut_area = cost.first;
            lut_delay = cost.second;
          }
        }
        else
        {
          std::tie( lut_area, lut_delay ) = lut_cost( truth_tables[cut->func_id] );
        }
      }
      else
      {
        std::tie( lut_area, lut_delay ) = lut_cost( cut.size() );
      }
    }
    else
    {
      lut_area = cut->data.lut_area;
      lut_delay = cut->data.lut_delay;
    }

    if constexpr ( ELA )
    {
      uint32_t delay{ 0 };
      for ( auto leaf : cut )
      {
        const auto& best_leaf_cut = cuts[leaf][0];
        delay = std::max( delay, best_leaf_cut->data.delay );
      }

      cut->data.delay = lut_delay + delay;
      cut->data.lut_area = lut_area;
      cut->data.lut_delay = lut_delay;
      if ( ps.edge_optimization )
      {
        cut->data.area_flow = static_cast<float>( cut_ref( cut ) );
        cut->data.edge_flow = static_cast<float>( cut_edge_deref( cut ) );
      }
      else
      {
        cut->data.area_flow = static_cast<float>( cut_measure_mffc( cut ) );
        cut->data.edge_flow = 0;
      }
    }
    else
    {
      uint32_t delay{ 0 };

      float area_flow = static_cast<float>( lut_area );
      float edge_flow = cut.size();

      for ( auto leaf : cut )
      {
        const auto& best_leaf_cut = cuts[leaf][0];
        delay = std::max( delay, best_leaf_cut->data.delay );
        if ( node_match[leaf].map_refs > 0 && leaf != 0 )
        {
          area_flow += best_leaf_cut->data.area_flow / node_match[leaf].est_refs;
          edge_flow += best_leaf_cut->data.edge_flow / node_match[leaf].est_refs;
        }
        else
        {
          area_flow += best_leaf_cut->data.area_flow;
          edge_flow += best_leaf_cut->data.edge_flow;
        }
      }

      cut->data.delay = lut_delay + delay;
      cut->data.lut_area = lut_area;
      cut->data.lut_delay = lut_delay;
      cut->data.area_flow = area_flow;
      cut->data.edge_flow = edge_flow;
    }
  }

  void compute_cut_data_share( cut_t& cut )
  {
    uint32_t delay{ 0 };
    float area_flow = static_cast<float>( cut->data.lut_area );
    float edge_flow = cut.size();

    for ( auto leaf : cut )
    {
      const auto& best_leaf_cut = cuts[leaf][0];
      delay = std::max( delay, best_leaf_cut->data.delay );
      /* flow contribution is added only for not shared leaves */
      if ( node_match[leaf].map_refs == 0 && leaf != 0 )
      {
        area_flow += best_leaf_cut->data.area_flow / node_match[leaf].est_refs;
        edge_flow += best_leaf_cut->data.edge_flow / node_match[leaf].est_refs;
      }
    }

    cut->data.delay = cut->data.lut_delay + delay;
    cut->data.area_flow = area_flow;
    cut->data.edge_flow = edge_flow;
  }

  template<bool ELA>
  void compute_wide_cut_data( wide_cut_t& cut, node const& n )
  {
    if constexpr ( ELA )
    {
      if ( ps.edge_optimization )
      {
        cut->data.area_flow = static_cast<float>( cut_ref<wide_cut_t>( cut ) );
        cut->data.edge_flow = static_cast<float>( cut_edge_deref<wide_cut_t>( cut ) );
      }
      else
      {
        cut->data.area_flow = static_cast<float>( cut_measure_mffc<wide_cut_t>( cut ) );
        cut->data.edge_flow = 0;
      }
    }
    else
    {
      // uint32_t delay{ 0 };
      float area_flow = static_cast<float>( cut->data.lut_area );
      float edge_flow = cut->data.lut_edges;

      for ( auto leaf : cut )
      {
        const auto& best_leaf_cut = cuts[leaf][0];
        if ( node_match[leaf].map_refs > 0 && leaf != 0 )
        {
          area_flow += best_leaf_cut->data.area_flow / node_match[leaf].est_refs;
          edge_flow += best_leaf_cut->data.edge_flow / node_match[leaf].est_refs;
        }
        else
        {
          area_flow += best_leaf_cut->data.area_flow;
          edge_flow += best_leaf_cut->data.edge_flow;
        }
      }

      cut->data.area_flow = area_flow;
      cut->data.edge_flow = edge_flow;
    }
  }

  void add_zero_cut( uint32_t index, bool phase )
  {
    auto& cut = cuts[index].add_cut( &index, &index ); /* fake iterator for emptyness */

    if constexpr ( StoreFunction )
    {
      if ( phase )
        cut->func_id = 1;
      else
        cut->func_id = 0;
    }
  }

  void add_unit_cut( uint32_t index )
  {
    auto& cut = cuts[index].add_cut( &index, &index + 1 );

    if constexpr ( StoreFunction )
    {
      cut->func_id = 2;
    }
  }

  uint32_t get_cut_delay_lcuts( cut_t const& cut, uint32_t lcut_i )
  {
    uint32_t cut_delay = cut->data.delay;

    if ( cut.size() == 1 )
    {
      cut_delay = ( *lcuts[lcut_i] )[0]->data.delay + 1;
    }

    return cut_delay;
  }

  inline bool fast_support_minimization( TT& tt, cut_t& res )
  {
    uint32_t support = 0u;
    uint32_t support_size = 0u;
    for ( uint32_t i = 0u; i < tt.num_vars(); ++i )
    {
      if ( kitty::has_var( tt, i ) )
      {
        support |= 1u << i;
        ++support_size;
      }
    }

    /* has not minimized support? */
    if ( ( support & ( support + 1u ) ) != 0u )
    {
      return false;
    }

    /* variables not in the support are the most significative */
    if ( support_size != res.size() )
    {
      std::vector<uint32_t> leaves( res.begin(), res.begin() + support_size );
      res.set_leaves( leaves.begin(), leaves.end() );
      tt = kitty::shrink_to( tt, support_size );
    }

    return true;
  }

  uint32_t compute_truth_table( uint32_t index, std::vector<cut_t const*> const& vcuts, cut_t& res )
  {
    // stopwatch t( st.cut_enumeration_st.time_truth_table ); /* runtime optimized */

    std::vector<TT> tt( vcuts.size() );
    auto i = 0;
    for ( auto const& cut : vcuts )
    {
      tt[i] = kitty::extend_to( truth_tables[( *cut )->func_id], res.size() );
      compute_truth_table_support( *cut, res, tt[i] );
      ++i;
    }

    auto tt_res = ntk.compute( ntk.index_to_node( index ), tt.begin(), tt.end() );

    if ( ps.cut_enumeration_ps.minimize_truth_table && !fast_support_minimization( tt_res, res ) )
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
        return truth_tables.insert( tt_res_shrink );
      }
    }

    return truth_tables.insert( tt_res );
  }

  void compute_mffcs_mapping()
  {
    ntk.clear_mapping();

    /* map POs */
    ntk.foreach_co( [&]( auto const& f ) {
      node const& n = ntk.get_node( f );
      if ( ntk.is_ci( n ) || ntk.is_constant( n ) )
        return;

      compute_mffc_mapping_node( n );
    } );

    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      node const& n = *it;
      if ( ntk.is_ci( n ) || ntk.is_constant( n ) )
        continue;
      if ( ntk.fanout_size( n ) <= 1 ) /* it should be unnecessary */
        continue;

      /* create MFFC cut */
      compute_mffc_mapping_node( n );
    }

    st.area = area;
    st.delay = delay;
    st.edges = edges;

    if ( ps.verbose )
    {
      std::stringstream stats;
      stats << fmt::format( "[i] Area MFFC: Delay = {:8d}  Area = {:8d}  Edges = {:8d}  Cuts = {:8d}\n", delay, area, edges, cuts_total );
      st.round_stats.push_back( stats.str() );
    }
  }

  void compute_mffc_mapping_node( node const& n )
  {
    uint32_t lut_area, lut_delay;

    /* create FC cut */
    std::vector<node> inner, leaves;
    ntk.incr_trav_id();
    get_fc_nodes_rec( n, inner );

    /* extract leaves */
    for ( auto const& g : inner )
    {
      ntk.foreach_fanin( g, [&]( auto const& f ) {
        if ( ntk.visited( ntk.get_node( f ) ) != ntk.trav_id() && !ntk.is_constant( ntk.get_node( f ) ) )
        {
          leaves.push_back( ntk.get_node( f ) );
          ntk.set_visited( ntk.get_node( f ), ntk.trav_id() );
        }
      } );
    }

    /* sort leaves in topo order */
    std::sort( leaves.begin(), leaves.end() );

    ntk.add_to_mapping( n, leaves.begin(), leaves.end() );

    delay = std::max( delay, static_cast<uint32_t>( leaves.size() ) );

    if constexpr ( StoreFunction )
    {
      default_simulator<TT> sim( leaves.size() );
      unordered_node_map<TT, Ntk> node_to_value( ntk );

      /* populate simulation values for constants */
      node_to_value[ntk.get_node( ntk.get_constant( false ) )] = sim.compute_constant( ntk.constant_value( ntk.get_node( ntk.get_constant( false ) ) ) );
      if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
      {
        node_to_value[ntk.get_node( ntk.get_constant( true ) )] = sim.compute_constant( ntk.constant_value( ntk.get_node( ntk.get_constant( true ) ) ) );
      }

      /* populate simulation values for leaves */
      uint32_t i = 0u;
      for ( auto const& g : leaves )
      {
        node_to_value[g] = sim.compute_pi( i++ );
      }

      /* simulate recursively */
      simulate_fc_rec( n, node_to_value );

      ntk.set_cell_function( n, node_to_value[n] );

      std::tie( lut_area, lut_delay ) = lut_cost( node_to_value[n] );
    }
    else
    {
      std::tie( lut_area, lut_delay ) = lut_cost( leaves.size() );
    }

    area += lut_area;
  }

  void get_fc_nodes_rec( node const& n, std::vector<node>& nodes )
  {
    if ( ntk.is_ci( n ) || ntk.is_constant( n ) )
      return;

    nodes.push_back( n );
    ntk.set_visited( n, ntk.trav_id() );

    /* expand cut for single fanout nodes */
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      auto g = ntk.get_node( f );
      if ( ntk.fanout_size( g ) == 1 )
      {
        get_fc_nodes_rec( g, nodes );
      }
    } );
  }

  void simulate_fc_rec( node const& n, unordered_node_map<TT, Ntk>& node_to_value )
  {
    std::vector<TT> fanin_values( ntk.fanin_size( n ) );

    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
      if ( !node_to_value.has( ntk.get_node( f ) ) )
      {
        simulate_fc_rec( ntk.get_node( f ), node_to_value );
      }

      fanin_values[i] = node_to_value[ntk.get_node( f )];
    } );

    node_to_value[n] = ntk.compute( n, fanin_values.begin(), fanin_values.end() );
  }

#pragma region Dump network
  klut_network create_lut_network()
  {
    klut_network res;
    node_map<signal<klut_network>, Ntk> node_to_signal( ntk );

    /* special map for output drivers to perform some optimizations */
    enum class driver_type : uint32_t
    {
      none = 0,
      pos = 1,
      neg = 2,
      mixed = 3
    };
    node_map<driver_type, Ntk> node_driver_type( ntk, driver_type::none );

    /* opposites are filled for nodes with mixed driver types, since they have
       two nodes in the network. */
    std::unordered_map<node, signal<klut_network>> opposites;

    /* initial driver types */
    ntk.foreach_po( [&]( auto const& f ) {
      switch ( node_driver_type[f] )
      {
      case driver_type::none:
        node_driver_type[f] = ntk.is_complemented( f ) ? driver_type::neg : driver_type::pos;
        break;
      case driver_type::pos:
        node_driver_type[f] = ntk.is_complemented( f ) ? driver_type::mixed : driver_type::pos;
        break;
      case driver_type::neg:
        node_driver_type[f] = ntk.is_complemented( f ) ? driver_type::neg : driver_type::mixed;
        break;
      case driver_type::mixed:
      default:
        break;
      }
    } );

    /* TODO: change -- it could be that internal nodes also point to an output driver node */
    ntk.foreach_node( [&]( auto const n ) {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) || node_match[ntk.node_to_index( n )].map_refs == 0 )
        return;

      auto const& best_cut = cuts[ntk.node_to_index( n )][0];
      for ( auto const& l : best_cut )
      {
        if ( node_driver_type[ntk.index_to_node( l )] == driver_type::neg )
        {
          node_driver_type[ntk.index_to_node( l )] = driver_type::mixed;
        }
      }
    } );

    /* constants */
    auto add_constant_to_map = [&]( bool value ) {
      const auto n = ntk.get_node( ntk.get_constant( value ) );
      switch ( node_driver_type[n] )
      {
      default:
      case driver_type::none:
      case driver_type::pos:
        node_to_signal[n] = res.get_constant( value );
        break;

      case driver_type::neg:
        node_to_signal[n] = res.get_constant( !value );
        break;

      case driver_type::mixed:
        node_to_signal[n] = res.get_constant( value );
        opposites[n] = res.get_constant( !value );
        break;
      }
    };

    add_constant_to_map( false );
    if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
    {
      add_constant_to_map( true );
    }

    /* primary inputs */
    ntk.foreach_pi( [&]( auto n ) {
      signal<klut_network> res_signal;
      switch ( node_driver_type[n] )
      {
      default:
      case driver_type::none:
      case driver_type::pos:
        res_signal = res.create_pi();
        node_to_signal[n] = res_signal;
        break;

      case driver_type::neg:
        res_signal = res.create_pi();
        node_to_signal[n] = res.create_not( res_signal );
        break;

      case driver_type::mixed:
        res_signal = res.create_pi();
        node_to_signal[n] = res_signal;
        opposites[n] = res.create_not( node_to_signal[n] );
        break;
      }
    } );

    /* TODO: add sequential compatibility */
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_ci( n ) || ntk.is_constant( n ) )
        continue;

      const auto index = ntk.node_to_index( n );
      if ( node_match[index].map_refs == 0 )
        continue;

      auto const& best_cut = cuts[index][0];

      /* if wide cut, perform decomposition */
      kitty::dynamic_truth_table tt;
      std::vector<signal<klut_network>> children;
      // if ( best_cut.size() > ps.cut_enumeration_ps.cut_size )
      // {
      //   std::tie( tt, children ) = create_lut_decompose( res, color_ntk, n, node_to_signal );
      // }
      // else
      // {
      if constexpr ( has_foreach_choice_v<Ntk> )
        std::tie( tt, children ) = create_lut_choice( n, node_to_signal );
      else
        std::tie( tt, children ) = create_lut( n, node_to_signal );
      // }

      switch ( node_driver_type[n] )
      {
      default:
      case driver_type::none:
      case driver_type::pos:
        node_to_signal[n] = res.create_node( children, tt );
        break;

      case driver_type::neg:
        node_to_signal[n] = res.create_node( children, ~tt );
        break;

      case driver_type::mixed:
        node_to_signal[n] = res.create_node( children, tt );
        opposites[n] = res.create_node( children, ~tt );
        break;
      }
    }

    /* outputs */
    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) && node_driver_type[f] == driver_type::mixed )
        res.create_po( opposites[ntk.get_node( f )] );
      else
        res.create_po( node_to_signal[f] );
    } );

    st.area = area;
    st.delay = delay;
    st.edges = edges;

    return res;
  }

  inline lut_info create_lut( node const& n, node_map<signal<klut_network>, Ntk>& node_to_signal )
  {
    auto const& best_cut = cuts[ntk.node_to_index( n )][0];

    std::vector<signal<klut_network>> children;
    for ( auto const& l : best_cut )
    {
      children.push_back( node_to_signal[ntk.index_to_node( l )] );
    }

    kitty::dynamic_truth_table tt;

    if constexpr ( StoreFunction )
    {
      tt = truth_tables[best_cut->func_id];
    }
    else
    {
      /* recursively compute the function for each choice until success */
      ntk.incr_trav_id();
      unordered_node_map<kitty::dynamic_truth_table, Ntk> node_to_value( ntk );

      /* add constants */
      node_to_value[ntk.get_node( ntk.get_constant( false ) )] = kitty::dynamic_truth_table( best_cut.size() );
      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
      if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
      {
        node_to_value[ntk.get_node( ntk.get_constant( true ) )] = ~kitty::dynamic_truth_table( best_cut.size() );
        ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
      }

      /* add leaves */
      uint32_t ctr = 0;
      for ( uint32_t leaf : best_cut )
      {
        kitty::dynamic_truth_table tt_leaf( best_cut.size() );
        kitty::create_nth_var( tt_leaf, ctr++ );
        node_to_value[ntk.index_to_node( leaf )] = tt_leaf;
        ntk.set_visited( ntk.index_to_node( leaf ), ntk.trav_id() );
      }

      /* recursively compute the function */     
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        compute_function_rec( ntk.get_node( f ), node_to_value );
      } );

      std::vector<kitty::dynamic_truth_table> tts;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        tts.push_back( node_to_value[ntk.get_node( f )] );
      } );
      tt = ntk.compute( n, tts.begin(), tts.end() );
    }

    return { tt, children };
  }

  inline lut_info create_lut_choice( node const& n, node_map<signal<klut_network>, Ntk>& node_to_signal )
  {
    auto const& best_cut = cuts[ntk.node_to_index( n )][0];

    std::vector<signal<klut_network>> children;
    for ( auto const& l : best_cut )
    {
      children.push_back( node_to_signal[ntk.index_to_node( l )] );
    }

    kitty::dynamic_truth_table tt;

    if constexpr ( StoreFunction )
    {
      tt = truth_tables[best_cut->func_id];
    }
    else
    {
      /* recursively compute the function for each choice until success */
      ntk.incr_trav_id();
      unordered_node_map<kitty::dynamic_truth_table, Ntk> node_to_value( ntk );

      /* add constants */
      node_to_value[ntk.get_node( ntk.get_constant( false ) )] = kitty::dynamic_truth_table( best_cut.size() );
      ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
      if ( ntk.get_node( ntk.get_constant( false ) ) != ntk.get_node( ntk.get_constant( true ) ) )
      {
        node_to_value[ntk.get_node( ntk.get_constant( true ) )] = ~kitty::dynamic_truth_table( best_cut.size() );
        ntk.set_visited( ntk.get_node( ntk.get_constant( true ) ), ntk.trav_id() );
      }

      /* add leaves */
      uint32_t ctr = 0;
      for ( uint32_t leaf : best_cut )
      {
        kitty::dynamic_truth_table tt_leaf( best_cut.size() );
        kitty::create_nth_var( tt_leaf, ctr++ );
        node_to_value[ntk.index_to_node( leaf )] = tt_leaf;
        ntk.set_visited( ntk.index_to_node( leaf ), ntk.trav_id() );
      }

      /* recursively try to compute the function for each choice */
      bool global_success = false;
      ntk.foreach_choice( n, [&]( auto const& g ) {
        bool success = true;
        ntk.foreach_fanin( g, [&]( auto const& f ) {
          success &= compute_function_choice_rec( ntk.get_node( f ), node_to_value );
          return success;
        } );

        /* move to the next choice */
        if ( !success )
          return true;

        std::vector<kitty::dynamic_truth_table> tts;
        ntk.foreach_fanin( g, [&]( auto const& f ) {
          tts.push_back( node_to_value[ntk.get_node( f )] );
        } );
        tt = ntk.compute( g, tts.begin(), tts.end() );

        /* get the correct phase */
        signal<Ntk> repr = ntk.get_choice_representative_signal( g );
        if ( ntk.is_complemented( repr ) )
          tt = ~tt;

        global_success = true;
        return false;
      } );
      assert( global_success );
    }

    return { tt, children };
  }

  void compute_function_rec( node const& n, unordered_node_map<kitty::dynamic_truth_table, Ntk>& node_to_value )
  {
    if ( ntk.visited( n ) == ntk.trav_id() )
    {
      assert( node_to_value.has( n ) );
      return;
    }

    assert( !ntk.is_pi( n ) );
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      compute_function_rec( ntk.get_node( f ), node_to_value );
    } );
    
    /* compute the function */
    std::vector<kitty::dynamic_truth_table> tts;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      tts.push_back( node_to_value[ntk.get_node( f )] );
    } );

    node_to_value[n] = ntk.compute( n, tts.begin(), tts.end() );
  }

  bool compute_function_choice_rec( node const& n, unordered_node_map<kitty::dynamic_truth_table, Ntk>& node_to_value )
  {
    bool choice_success = false;
    ntk.foreach_choice( n, [&]( auto const& g ) {
      if ( ntk.visited( g ) == ntk.trav_id() )
      {
        choice_success = node_to_value.has( g );
        return false;
      }

      if ( ntk.is_pi( g ) )
        return false;

      ntk.set_visited( g, ntk.trav_id() );

      bool success = true;
      ntk.foreach_fanin( g, [&]( auto const& f ) {
        success &= compute_function_choice_rec( ntk.get_node( f ), node_to_value );
        return success;
      } );

      if ( !success )
        return true;
      
      /* compute the function */
      choice_success = true;
      std::vector<kitty::dynamic_truth_table> tts;
      ntk.foreach_fanin( g, [&]( auto const& f ) {
        tts.push_back( node_to_value[ntk.get_node( f )] );
      } );
      auto tt = ntk.compute( g, tts.begin(), tts.end() );

      /* get the correct phase */
      signal<Ntk> repr = ntk.get_choice_representative_signal( g );
      if ( ntk.is_complemented( repr ) )
        node_to_value[n] = ~tt;
      else
        node_to_value[n] = tt;
      return false;
    } );

    return choice_success;
  }

  void derive_mapping()
  {
    ntk.clear_mapping();

    for ( auto const& n : topo_order )
    {
      if ( ntk.is_ci( n ) || ntk.is_constant( n ) )
        continue;

      const auto index = ntk.node_to_index( n );
      if ( node_match[index].map_refs == 0 )
        continue;

      std::vector<node> nodes;
      auto const& best_cut = cuts[index][0];

      for ( auto const& l : best_cut )
      {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      if constexpr ( StoreFunction )
      {
        ntk.set_cell_function( n, truth_tables[best_cut->func_id] );
      }
    }

    st.area = area;
    st.delay = delay;
    st.edges = edges;
  }
#pragma endregion

private:
  Ntk& ntk;
  lut_map_params const& ps;
  lut_map_stats& st;

  uint32_t iteration{ 0 };        /* current mapping iteration */
  uint32_t area_iteration{ 0 };   /* current area iteration */
  uint32_t delay{ 0 };            /* current delay of the mapping */
  uint32_t area{ 0 };             /* current area of the mapping */
  uint32_t edges{ 0 };            /* current edges of the mapping */
  uint32_t cuts_total{ 0 };       /* current computed cuts */
  uint32_t multi_gates{ 0 };       /* number of detected multi-input ANDs */
  const float epsilon{ 0.005f };  /* epsilon */
  LUTCostFn lut_cost{};

  std::vector<node> topo_order;
  std::vector<node> topo_order2;
  std::vector<uint32_t> tmp_visited;
  std::vector<node_lut> node_match;

  std::deque<cut_set_t> cuts;   /* compressed representation of cuts */
  wide_cut_set_t wide_set;      /* compressed representation of wide cuts */
  cut_merge_t lcuts;            /* cut merger container */
  wide_cut_merge_t lwidecuts;   /* wide cut merger container */
  tt_cache truth_tables;        /* cut truth tables */
  cost_cache truth_tables_cost; /* truth tables cost */

  adaptive_buckets_t dcuts1;        /* cut set for multi-input modes decomposition */
  adaptive_buckets_t dcuts2;        /* cut merger for multi-input modes decomposition */
  std::vector<uint32_t> bins;       /* bin containers for bin-packing */
  std::vector<bin_data> bins_ext;   /* delay for bin containers */
  std::vector<uint32_t> packs;      /* packs for bin-packing */
  pack_index_list pack_indexes;     /* indexes used to sort packs */
  bin_merger_t lbins;               /* merger containers of bins */
};
#pragma endregion

} /* namespace detail */

/*! \brief LUT mapper.
 *
 * This function implements a LUT mapping algorithm.  It is controlled by one
 * template argument `ComputeTruth` (defaulted to `false`) which controls
 * whether the LUT function is computed or the mapping is structural. In the
 * former case, truth tables are computed during cut enumeration,
 * which requires more runtime.
 * 
 * This function returns a k-LUT network.
 *
 * The template `LUTCostFn` sets the cost function to evaluate depth and
 * size of a truth table given its support size if `ComputeTruth` is set
 * to false, or its function if `ComputeTruth` is set to true.
 * 
 * This implementation offers more options such as delay oriented mapping
 * and edges minimization compared to the command `lut_mapping`.
 *
 * **Required network functions:**
 * - `size`
 * - `is_ci`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_co`
 * - `foreach_node`
 * - `fanout_size`
 */
template<class Ntk, bool ComputeTruth = false, class LUTCostFn = lut_unitary_cost>
klut_network lut_map( Ntk& ntk, lut_map_params ps = {}, lut_map_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_ci_v<Ntk>, "Ntk does not implement the foreach_ci method" );
  static_assert( has_foreach_co_v<Ntk>, "Ntk does not implement the foreach_co method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );

  lut_map_stats st;
  klut_network klut;

  if ( ps.multi_decomposition )
  {
    /* copy network */
    Ntk ntk_copy = cleanup_dangling( ntk );
    choice_view<Ntk> choice_ntk{ ntk_copy };

    /* deactivate cut expansion */
    ps.cut_expansion = false;

    detail::lut_map_impl<choice_view<Ntk>, false, LUTCostFn> p( choice_ntk, ps, st );
    klut = p.run();
  }
  else
  {
    detail::lut_map_impl<Ntk, ComputeTruth, LUTCostFn> p( ntk, ps, st );
    klut = p.run();
  }

  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }

  return klut;
}

/*! \brief LUT mapper inplace.
 *
 * This function implements a LUT mapping algorithm.  It is controlled by one
 * template argument `StoreFunction` (defaulted to `false`) which controls
 * whether the LUT function is stored in the mapping. In that case
 * truth tables are computed during cut enumeration, which requires more
 * runtime.
 * 
 * The input network must be wrapped in a `mapping_view`. The computed mapping
 * is stored in the view. In this version, some features of the mapper are
 * disabled, such as on-the-fly decompositions, due to incompatibility.
 *
 * The template `LUTCostFn` sets the cost function to evaluate depth and
 * size of a truth table given its support size, if `StoreFunction` is set
 * to false, or its function, if `StoreFunction` is set to true.
 * 
 * This implementation offers more options such as delay oriented mapping
 * and edges minimization compared to the command `lut_mapping`.
 *
 * **Required network functions:**
 * - `size`
 * - `is_ci`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_co`
 * - `foreach_node`
 * - `fanout_size`
 * - `clear_mapping`
 * - `add_to_mapping`
 * - `set_lut_function` (if `StoreFunction` is true)
 *
 *
   \verbatim embed:rst

   .. note::

      The implementation of this algorithm was inspired by the LUT
      mapping command ``&if`` in ABC.
   \endverbatim
 */
template<class Ntk, bool StoreFunction = false, class LUTCostFn = lut_unitary_cost>
void lut_map_inplace( Ntk& ntk, lut_map_params const& ps = {}, lut_map_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_ci_v<Ntk>, "Ntk does not implement the foreach_ci method" );
  static_assert( has_foreach_co_v<Ntk>, "Ntk does not implement the foreach_co method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );

  lut_map_stats st;
  detail::lut_map_impl<Ntk, StoreFunction, LUTCostFn> p( ntk, ps, st );
  p.run_inplace();

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
