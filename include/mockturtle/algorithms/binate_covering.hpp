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
  \file binate_covering.hpp
  \brief Binate covering

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>

#include <fmt/format.h>
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>

#include "../utils/stopwatch.hpp"
#include "cut_enumeration.hpp"

namespace mockturtle
{

/*! \brief Parameters for binate_covering.
 *
 * The data structure `binate_covering_params` holds configurable parameters
 * with default arguments for `binate_covering`.
 */
struct binate_covering_params
{
  binate_covering_params()
  {
    cut_enumeration_ps.cut_size = 4;
    cut_enumeration_ps.cut_limit = 8;
    cut_enumeration_ps.minimize_truth_table = false;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut size is 4, the default cut limit is 8.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Debug mode. */
  bool debug{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for binate_covering.
 *
 * The data structure `binate_covering_stats` provides data collected by running
 * `binate_covering`.
 */
struct binate_covering_stats
{
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

template<class Ntk, bool StoreFunction>
class binate_covering_impl
{
public:
  using network_cuts_t = network_cuts<Ntk, StoreFunction, empty_cut_data>;
  using cut_t = typename network_cuts_t::cut_t;
  using covering_matrix_t = std::vector<kitty::partial_truth_table>;

public:
  binate_covering_impl( Ntk& ntk, binate_covering_params const& ps, binate_covering_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        cov(),
        cov_t(),
        map_refs( ntk.size(), 0 ),
        cuts( cut_enumeration<Ntk, StoreFunction>( ntk, ps.cut_enumeration_ps ) )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    /* initialize data structure */
    init();

    /* add cuts to the covering matrix */
    populate();

    /* solve the binate covering problem */
    solve();

    /* write solution */
    derive_mapping();
  }

private:
  void init()
  {
    /* get sizes */
    num_rows = ntk.num_gates();
    num_columns = cuts.total_cuts();
    offset = ntk.size() - ntk.num_gates();

    /* allocate matrix */
    cov = covering_matrix_t( num_rows );
    cov_t = covering_matrix_t( num_columns );

    /* allocate columns */
    for ( uint32_t i = 0; i < num_rows; ++i )
    {
      cov[i] = kitty::partial_truth_table( num_columns );
    }

    /* allocate columns */
    for ( uint32_t i = 0; i < num_columns; ++i )
    {
      cov_t[i] = kitty::partial_truth_table( num_rows );
    }

    /* allocate solution */
    best_solution = kitty::partial_truth_table( num_columns );
  }

  void populate()
  {
    uint32_t pcol = 0;
    ntk.foreach_gate( [&]( auto const& n ) {
      uint32_t index = ntk.node_to_index( n );
      for ( cut_t const* cut : cuts.cuts( index ) )
      {
        /* exclude trivial cut */
        if ( cut->size() == 1 && *( cut->begin() ) == index )
          continue;

        add_cut( n, *cut, pcol++ );
      }
    } );

    if ( ps.debug )
    {
      print_cov();
    }
  }

  void solve()
  {
    /* TO implement */
    return;
  }

#pragma region Initialization routines
  void add_cut( node<Ntk> const& n, cut_t const& cut, uint32_t pcol )
  {
    /* add the volume of the cut */

    ntk.incr_trav_id();
    for ( auto l : cut )
    {
      ntk.set_visited( ntk.index_to_node( l ), ntk.trav_id() );
    }

    add_cut_volume_rec( n, pcol );
  }

  void add_cut_volume_rec( node<Ntk> const& n, uint32_t pcol )
  {
    if ( ntk.is_pi( n ) )
      return;

    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    /* mark as visited */
    ntk.set_visited( n, ntk.trav_id() );

    /* add to volume */
    cov_table_add_bit( ntk.node_to_index( n ), pcol );
    cov_t_table_add_bit( pcol, ntk.node_to_index( n ) );

    /* recur */
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      add_cut_volume_rec( ntk.get_node( f ), pcol );
    } );
  }
#pragma endregion

  void derive_mapping()
  {
    ntk.clear_mapping();

    ntk.foreach_gate( [&]( auto const& n ) {
      const auto index = ntk.node_to_index( n );
      if ( map_refs[index] == 0 )
        return true;

      std::vector<node<Ntk>> nodes;
      for ( auto const& l : cuts.cuts( index ).best() )
      {
        nodes.push_back( ntk.index_to_node( l ) );
      }
      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );

      if constexpr ( StoreFunction )
      {
        ntk.set_cell_function( n, cuts.truth_table( cuts.cuts( index ).best() ) );
      }
      return true;
    } );
  }

  void print_cov()
  {
    ntk.foreach_gate( [&]( auto const& n, uint32_t i ) {
      std::cout << fmt::format( "n{}\t : ", ntk.node_to_index( n ) );
      kitty::print_binary( cov[i] );
      std::cout << "\n";
    } );
  }
private:
  inline void cov_table_add_bit( uint32_t r, uint32_t c ) { cov[r - offset]._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void cov_t_table_add_bit( uint32_t r, uint32_t c ) { cov_t[r]._bits[( c - offset ) >> 6] |= UINT64_C( 1 ) << ( ( c - offset ) & 0x3f ); }

private:
  Ntk& ntk;
  binate_covering_params const& ps;
  binate_covering_stats& st;

  uint32_t num_rows{ 0 };     /* number of rows */
  uint32_t num_columns{ 0 };  /* number of columns */
  uint32_t offset{ 0 };       /* access offset */
  uint32_t area{ 0 };         /* current area of the mapping */

  uint32_t best_cost{ UINT32_MAX };          /* tracks the cost of the best solution */
  kitty::partial_truth_table best_solution;  /* tracks the best solution: set of columns */

  covering_matrix_t cov;      /* covering matrix */
  covering_matrix_t cov_t;    /* covering matrix transposed */
  std::vector<uint32_t> map_refs;
  
  network_cuts_t cuts;
};

}; /* namespace detail */


template<class Ntk, bool StoreFunction = false>
void binate_covering( Ntk& ntk, binate_covering_params const& ps = {}, binate_covering_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_ci_v<Ntk>, "Ntk does not implement the is_ci method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_co_v<Ntk>, "Ntk does not implement the foreach_co method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_clear_mapping_v<Ntk>, "Ntk does not implement the clear_mapping method" );
  static_assert( has_add_to_mapping_v<Ntk>, "Ntk does not implement the add_to_mapping method" );
  static_assert( !StoreFunction || has_set_cell_function_v<Ntk>, "Ntk does not implement the set_cell_function method" );

  binate_covering_stats st;
  detail::binate_covering_impl<Ntk, StoreFunction> p( ntk, ps, st );
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
