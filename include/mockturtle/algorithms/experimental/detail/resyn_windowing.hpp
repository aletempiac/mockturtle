/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file resyn_windowing.hpp
  \brief A windowing engine for rewriting

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <iostream>
#include <optional>
#include <vector>

#include "../../../utils/algorithm.hpp"

namespace mockturtle
{

struct resyn_windowing_params
{
  /*! \brief Maximum number of gates to include in a window. */
  uint32_t max_gates{ 10 };

  /*! \brief Maximum fanout of a node to expand. */
  uint32_t skip_fanout_limit{ 5 };
};

namespace detail
{

template<class Ntk>
class resyn_windowing
{
private:
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit resyn_windowing( Ntk const& ntk, resyn_windowing_params const& ps )
      : ntk( ntk ), ps( ps ), leaves(), roots(), gates(), candidates()
  {
    leaves.reserve( ps.max_gates );
    roots.reserve( ps.max_gates );
    gates.reserve( ps.max_gates );
    candidates.reserve( ps.max_gates );
  }

  void compute_window( node const& pivot )
  {
    leaves.clear();
    roots.clear();
    gates.clear();

    /* add pivot to gates */
    ntk.incr_trav_id();
    ntk.set_visited( pivot, ntk.trav_id() );
    gates.push_back( pivot );

    if ( ps.max_gates < 2 )
      return;

    /* decrement fanout size of leaves */
    ntk.foreach_fanin( pivot, [&]( auto const& f ) {
      ntk.decr_fanout_size( ntk.get_node( f ) );
    } );

    /* increment traverse ID */
    ntk.incr_trav_id();

    /* TODO: modify the code to include boundary buffers and inverters in the window */
    /* add iteratively nodes to the window */
    std::optional<node> next;
    while ( gates.size() < ps.max_gates && ( next = find_next_pivot() ) )
    {
      assert( ntk.visited( *next ) < ntk.trav_id() - 1 );
      gates.push_back( *next );
      ntk.set_visited( *next, ntk.trav_id() - 1 );

      /* decrement fanout size of leaves */
      ntk.foreach_fanin( *next, [&]( auto const& f ) {
        ntk.decr_fanout_size( ntk.get_node( f ) );
      } );
    }

    /* restore fanout counts */
    for ( node const& n : gates )
    {
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        ntk.incr_fanout_size( ntk.get_node( f ) );
      } );
    }

    /* collect roots, leaves, and gates in topo order */
    collect_roots();
    collect_nodes();

    assert( gates.size() <= ps.max_gates );
    return;
  }

  std::vector<node> const& get_gates() const
  {
    return gates;
  }

  std::vector<node> const& get_leaves() const
  {
    return leaves;
  }

  std::vector<signal> const& get_roots() const
  {
    return roots;
  }

  uint32_t const num_gates() const
  {
    return gates.size();
  }

  uint32_t const num_leaves() const
  {
    return leaves.size();
  }

  uint32_t const num_roots() const
  {
    return leaves.size();
  }

  std::array<uint64_t, 2> const& get_hash() const
  {
    return hash;
  }

  void report_info( std::ostream& out = std::cout )
  {
    out << fmt::format( "[i] W: I = {};\t G = {};\t O = {}\n", leaves.size(), gates.size(), roots.size() );
  }

private:
  std::optional<node> find_next_pivot()
  {
    candidates.clear();
    auto best = candidates.cend();

    /* marks meaning */
    /* visited == trav_id - 1 --> in the window */
    /* visited == trav_id --> candidate in the TFI of the window */

    /* look for MFFC nodes */
    do
    {
      for ( node const& n : gates )
      {
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          node const& g = ntk.get_node( f );
          if ( ntk.visited( g ) < ntk.trav_id() - 1 && ntk.fanout_size( g ) == 0 && !ntk.is_ci( g ) )
          {
            candidates.push_back( g );
            ntk.set_visited( g, ntk.trav_id() );
          }
        } );
      }

      if ( !candidates.empty() )
      {
        /* select the best candidate */
        best = max_element_unary(
            candidates.begin(), candidates.end(),
            [&]( auto const& cand ) {
              auto cnt{ 0 };
              ntk.foreach_fanin( cand, [&]( auto const& f ) {
                cnt += ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() ? 1 : 0;
              } );
              return cnt;
            },
            -1 );
        break;
      }

      /* add all the input candidates */
      for ( node const& n : gates )
      {
        ntk.foreach_fanin( n, [&]( auto const& f ) {
          node const& g = ntk.get_node( f );
          if ( ntk.visited( g ) < ntk.trav_id() - 1 && !ntk.is_ci( g ) )
          {
            candidates.push_back( g );
            ntk.set_visited( g, ntk.trav_id() );
          }
        } );
      }

      /* add all the output candidates */
      for ( node const& n : gates )
      {
        if ( ntk.fanout_size( n ) == 0 || ntk.fanout_size( n ) > ps.skip_fanout_limit )
          continue;

        /* Output candidate member of the MFFC */
        auto const& fanout_v = ntk.fanout( n );
        if ( ntk.fanout_size( n ) == 1 && fanout_v.size() == 1 && ntk.visited( fanout_v.front() ) < ntk.trav_id() - 1 )
        {
          candidates.push_back( fanout_v.front() );
          best = candidates.cend() - 1;
          break;
        }

        /* does not mark the nodes (TFI has the priority to avoid reconvergences passing from outside the window) */
        std::copy_if( fanout_v.begin(), fanout_v.end(),
                      std::back_inserter( candidates ),
                      [&]( auto const& g ) {
                        return ntk.visited( g ) < ntk.trav_id() - 1;
                      } );
      }

      if ( !candidates.empty() )
      {
        /* select the best candidate */
        best = max_element_unary(
            candidates.begin(), candidates.end(),
            [&]( auto const& cand ) {
              auto cnt{ 0 };
              ntk.foreach_fanin( cand, [&]( auto const& f ) {
                cnt += ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() ? 1 : 0;
              } );
              return cnt;
            },
            -1 );
      }
    } while ( false );

    if ( candidates.empty() )
    {
      return std::nullopt;
    }

    assert( best != candidates.end() );

    /* reset the marks */
    for ( auto const& n : candidates )
    {
      ntk.set_visited( n, ntk.trav_id() - 2 );
    }

    return *best;
  }

  void collect_roots()
  {
    hash[0] = 0;

    for ( node const& n : gates )
    {
      /* equivalent to is_po() */
      if ( ntk.fanout_size( n ) != ntk.fanout( n ).size() )
      {
        roots.push_back( ntk.make_signal( n ) );
        hash[0] |= UINT64_C( 1 ) << ( ntk.node_to_index( n ) % 64 );
        continue;
      }

      ntk.foreach_fanout( n, [&]( auto const& g ) {
        if ( ntk.visited( g ) < ntk.trav_id() - 1 )
        {
          roots.push_back( ntk.make_signal( n ) );
          return false;
        }
        return true;
      } );
    }
  }

  void collect_nodes()
  {
    auto prev_size = gates.size();
    gates.clear();
    hash[1] = 0;

    /* collects gates and leaves */
    for ( signal const& s : roots )
    {
      node const& n = ntk.get_node( s );
      if ( ntk.visited( n ) == ntk.trav_id() )
        continue;

      collect_nodes_rec( n );
    }

    assert( gates.size() == prev_size );
  }

  void collect_nodes_rec( node const& n )
  {
    assert( ntk.visited( n ) != ntk.trav_id() );
    ntk.set_visited( n, ntk.trav_id() );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      node const& g = ntk.get_node( f );

      if ( ntk.visited( g ) < ntk.trav_id() - 1 )
      {
        /* leaf */
        leaves.push_back( g );
        ntk.set_visited( g, ntk.trav_id() );
        hash[1] |= UINT64_C( 1 ) << ( ntk.node_to_index( g ) % 64 );
      }
      else if ( ntk.visited( g ) == ntk.trav_id() - 1 )
      {
        /* gate */
        collect_nodes_rec( g );
      }
    } );

    gates.push_back( n );
  }

private:
  Ntk const& ntk;
  resyn_windowing_params const& ps;

  std::vector<node> leaves;
  std::vector<signal> roots;
  std::vector<node> gates;
  std::vector<node> candidates;

  std::array<uint64_t, 2> hash;
};

} /* namespace detail */

} /* namespace mockturtle */