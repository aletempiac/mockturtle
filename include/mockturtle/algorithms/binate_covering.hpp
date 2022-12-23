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

#include <chrono>
#include <cstdint>

#include <fmt/format.h>
#include <kitty/partial_truth_table.hpp>
#include <kitty/print.hpp>

#include "../utils/stopwatch.hpp"
#include "cut_enumeration.hpp"
#include "cut_enumeration/mf_cut.hpp"

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

  /*! \brief Initial bound for the search. */
  uint32_t bound{ UINT32_MAX };

  /*! \brief Timeout limit (seconds) */
  float timeout{ 10.0 };

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

template<class Ntk, bool StoreFunction, typename CutData>
class binate_covering_impl
{
public:
  using clock = typename std::chrono::steady_clock;
  using time_point = typename clock::time_point;
  using network_cuts_t = network_cuts<Ntk, StoreFunction, CutData>;
  using cut_t = typename network_cuts_t::cut_t;
  using covering_matrix_t = std::vector<kitty::partial_truth_table>;

public:
  binate_covering_impl( Ntk& ntk, binate_covering_params const& ps, binate_covering_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st ),
        cov(),
        cov_t(),
        cuts( cut_enumeration<Ntk, StoreFunction, CutData>( ntk, ps.cut_enumeration_ps ) ),
        time_begin( clock::now() )
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
    constraints = covering_matrix_t( num_columns );

    /* allocate columns */
    for ( uint32_t i = 0; i < num_rows; ++i )
    {
      cov[i] = kitty::partial_truth_table( num_columns );
    }

    /* allocate columns */
    for ( uint32_t i = 0; i < num_columns; ++i )
    {
      cov_t[i] = kitty::partial_truth_table( num_rows );
      constraints[i] = kitty::partial_truth_table( num_rows );
    }

    /* allocate solution */
    current_solution = kitty::partial_truth_table( num_rows );
    current_coverage = kitty::partial_truth_table( num_rows );
    current_constraints = kitty::partial_truth_table( num_rows );
    best_solution = kitty::partial_truth_table( num_rows );

    /* get initial bounds */
    best_cost = ps.bound;
    if ( best_cost >= UINT32_MAX - num_rows )
      best_cost = UINT32_MAX - num_rows;
    upper_bound = num_rows;
  }

  void populate()
  {
    uint32_t pcol = 0;

    /* add cuts */
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

    /* add constraints */
    ntk.foreach_po( [&]( auto const& f ) {
      uint32_t index = ntk.node_to_index( ntk.get_node( f ) );
      if ( index >= offset )
        constraints_add_bit( ntk.node_to_index( ntk.get_node( f ) ) - offset );
    } );

    // if ( ps.debug )
    // {
    //   print_cov();
    // }
  }

  void solve()
  {
    bool res = solve_rec( num_rows - 1, num_columns - 1 );

    if ( ps.verbose )
    {
      if ( res )
        std::cout << fmt::format( "Best optimal solution with cost {}\n", best_cost );
      else
        std::cout << fmt::format( "Best sub-optimal solution with cost {}\n", best_cost );
    }
    if ( ps.debug )
    {
      kitty::print_binary( best_solution );
      std::cout << "\n";
    }
  }

#pragma region Solver
  bool solve_rec( uint32_t row_index, int col_index )
  {
    /* solution found */
    uint32_t rows_covered = kitty::count_ones( current_coverage );
    if ( rows_covered == num_rows )
    {
      evaluate_solution();
      return is_not_timeout();
    }

    /* bound branching */
    if ( current_cost >= best_cost - 1 )
      return is_not_timeout();
    // if ( current_cost >= best_cost + rows_covered - num_rows )
    //   return is_not_timeout();

    /* get current step choices */
    int cut_size = cuts.cuts( row_index + offset ).size() - 1;

    /* find row to cover */
    while ( !constraints_has_bit( row_index ) )
    {
      --row_index;
      col_index -= cut_size;
      cut_size = cuts.cuts( row_index + offset ).size() - 1;
    }
  
    /* save current solution; TODO: find alternative for efficiency */
    auto cov_step = current_coverage;
    auto constraints_step = current_constraints;

    /* select current bit in current selection */
    solution_add_bit( row_index );

    /* select a column and recur */
    for ( int i = col_index; i > col_index - cut_size; --i )
    {
      /* add selection and column to the solution */
      coverage_add_column( i );
      constraints_add_column( i );

      /* step evaluation */
      evaluate_step();

      /* branch */
      if ( !solve_rec( row_index - 1, col_index - cut_size ) )
        return false;

      /* remove column from solution */
      undo_step();
      current_coverage = cov_step;
      current_constraints = constraints_step;
    }

    /* remove bit from solution */
    solution_remove_bit( row_index );

    return true;
  }

  inline void evaluate_step()
  {
    current_cost += 1;
  }

  inline void undo_step()
  {
    current_cost -= 1;
  }

  inline void evaluate_solution()
  {
    if ( current_cost >= best_cost )
      return;

    best_cost = current_cost;
    best_solution = current_solution;

    if ( ps.debug )
    {
      std::cout << fmt::format( "New solution with cost {}\n", best_cost );
    }
  }
#pragma endregion;

#pragma region Initialization routines
  void add_cut( node<Ntk> const& n, cut_t const& cut, uint32_t pcol )
  {
    /* add the volume of the cut */
    ntk.incr_trav_id();
    for ( auto l : cut )
    {
      ntk.set_visited( ntk.index_to_node( l ), ntk.trav_id() );
      if ( l >= offset )
        constraints_table_add_bit( pcol, l - offset );
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
    cov_table_add_bit( ntk.node_to_index( n ) - offset, pcol );
    cov_t_table_add_bit( pcol, ntk.node_to_index( n ) - offset );

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
      if ( !best_solution_has_bit( index - offset ) )
        return true;

      ntk.incr_trav_id();
      std::vector<node<Ntk>> nodes;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        collect_leaves_rec( ntk.get_node( f ), nodes );
      } );

      ntk.add_to_mapping( n, nodes.begin(), nodes.end() );
      return true;
    } );
  }

  void collect_leaves_rec( node<Ntk> const& n, std::vector<node<Ntk>>& leaves )
  {
    const auto index = ntk.node_to_index( n );

    if ( ntk.visited( n ) == ntk.trav_id() )
      return;
    
    ntk.set_visited( n, ntk.trav_id() );

    if ( ntk.is_pi( n ) || best_solution_has_bit( index - offset ) )
    {
      leaves.push_back( n );
      return;
    }

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      collect_leaves_rec( ntk.get_node( f ), leaves );
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
  inline void cov_table_add_bit( uint32_t r, uint32_t c ) { cov[r]._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void cov_t_table_add_bit( uint32_t r, uint32_t c ) { cov_t[r]._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void constraints_table_add_bit( uint32_t r, uint32_t c ) { constraints[r]._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void solution_add_bit( uint32_t c ) { current_solution._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void constraints_add_bit( uint32_t c ) { current_constraints._bits[c >> 6] |= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline void solution_remove_bit( uint32_t c ) { current_solution._bits[c >> 6] ^= UINT64_C( 1 ) << ( c & 0x3f ); }
  inline bool solution_has_bit( uint32_t c ) { return ( current_solution._bits[c >> 6] & ( UINT64_C( 1 ) << ( c & 0x3f ) ) ) > 0; }
  inline bool coverage_has_bit( uint32_t c ) {  return ( current_coverage._bits[c >> 6] & ( UINT64_C( 1 ) << ( c & 0x3f ) ) ) > 0; }
  inline bool constraints_has_bit( uint32_t c ) {  return ( current_constraints._bits[c >> 6] & ( UINT64_C( 1 ) << ( c & 0x3f ) ) ) > 0; }
  inline bool best_solution_has_bit( uint32_t c ) { return ( best_solution._bits[c >> 6] & ( UINT64_C( 1 ) << ( c & 0x3f ) ) ) > 0; }
  inline void coverage_add_column( uint32_t r ) { current_coverage |= cov_t[r]; }
  inline void constraints_add_column( uint32_t r ) { current_constraints |= constraints[r]; }
  inline bool is_not_timeout() { return std::chrono::duration_cast<std::chrono::duration<float>>( clock::now() - time_begin ).count() < ps.timeout; }

private:
  Ntk& ntk;
  binate_covering_params const& ps;
  binate_covering_stats& st;

  uint32_t num_rows{ 0 };     /* number of rows */
  uint32_t num_columns{ 0 };  /* number of columns */
  uint32_t offset{ 0 };       /* access offset */
  uint32_t area{ 0 };         /* current area of the mapping */

  time_point time_begin;      /* tracks time */

  uint32_t upper_bound{ UINT32_MAX };             /* upper bound */
  uint32_t best_cost{ UINT32_MAX };               /* tracks the cost of the best solution */
  kitty::partial_truth_table best_solution;       /* tracks the best solution: set of rows */
  uint32_t current_cost{ 0 };                     /* tracks the current cost */
  kitty::partial_truth_table current_solution;    /* tracks the current solution */
  kitty::partial_truth_table current_coverage;    /* tracks the current solution node coverage */
  kitty::partial_truth_table current_constraints; /* tracks the current solution node coverage */

  covering_matrix_t cov;            /* covering matrix */
  covering_matrix_t cov_t;          /* covering matrix transposed */
  covering_matrix_t constraints;    /* constrain matrix */

  network_cuts_t cuts;
};

}; /* namespace detail */

template<class Ntk, bool StoreFunction = false, typename CutData = cut_enumeration_mf_cut>
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
  detail::binate_covering_impl<Ntk, StoreFunction, CutData> p( ntk, ps, st );
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
