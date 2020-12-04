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
  \file topo_view.hpp
  \brief Reimplements foreach_node to guarantee topological order

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cassert>
#include <vector>
#include <iostream>

#include "../networks/detail/foreach.hpp"
#include "../traits.hpp"

namespace mockturtle
{

/*! \brief Ensures topological order for of all nodes reachable from the outputs.
 *
 * Overrides the interface methods `foreach_node`, `foreach_gate`,
 * `size`, `num_gates`.
 *
 * This class computes *on construction* a topological order of the nodes which
 * are reachable from the outputs.  Constant nodes and primary inputs will also
 * be considered even if they are not reachable from the outputs.  Further,
 * constant nodes and primary inputs will be visited first before any gate node
 * is visited.  Constant nodes precede primary inputs, and primary inputs are
 * visited in the same order in which they were created.
 *
 * Since the topological order is computed only once when creating an instance,
 * this view disables changes to the network interface.  Also, since only
 * reachable nodes are traversed, not all network nodes may be called in
 * `foreach_node` and `foreach_gate`.
 *
 * **Required network functions:**
 * - `get_constant`
 * - `foreach_pi`
 * - `foreach_po`
 * - `foreach_fanin`
 * - `incr_trav_id`
 * - `set_visited`
 * - `trav_id`
 * - `visited`
 *
 * Example
 *
   \verbatim embed:rst

   .. code-block:: c++

      // create network somehow; aig may not be in topological order
      aig_network aig = ...;

      // create a topological view on the network
      topo_view aig_topo{aig};

      // call algorithm that requires topological order
      cut_enumeration( aig_topo );
   \endverbatim
 */
template<class Ntk>
class eq_classes
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  /*! \brief Default constructor.
   *
   * Constructs topological view on another network.
   */
  eq_classes( Ntk const& ntk ) : ntk( ntk )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
    static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_pi method" );

    init_eqclass();
  }


  void create_repr( node const& n1, node const& n2 )
  {
    if ( ntk.node_to_index( n1 ) == ntk.node_to_index( n2 ) ) {
      return;
    }

    auto rep1 = get_eqrepr( n1 );
    auto rep2 = get_eqrepr( n2 );

    if ( ntk.node_to_index( rep1 ) == ntk.node_to_index( rep2 ) ) {
      /* already in the same class */
      return;
    }

    /* merge the eq lists, set the representative as the node with lowest ID */
    if ( ntk.node_to_index( rep1 ) < ntk.node_to_index( rep2 ) )
    {
      eqrep[ntk.node_to_index( rep2 )] = ntk.get_node( eqnodes[ntk.node_to_index( rep1 )] );
      eqnodes[ntk.node_to_index( rep1 )] = eqnodes[ntk.node_to_index( rep2 )];
      eqnodes[ntk.node_to_index( rep2 )] = ntk.make_signal( rep2 );
    }
    else
    {
      eqrep[ntk.node_to_index( rep1 )] = ntk.get_node( eqnodes[ntk.node_to_index( rep2 )] );
      eqnodes[ntk.node_to_index( rep2 )] = eqnodes[ntk.node_to_index( rep1 )];
      eqnodes[ntk.node_to_index( rep1 )] = ntk.make_signal( rep1 );
    }
  }


  void create_repr( node const& n1, signal const& s2 )
  {
    auto const n2 = ntk.get_node( s2 );
    auto const id1 = ntk.node_to_index( n1 );
    auto const id2 = ntk.node_to_index( n2 );

    //std::cout << ntk.is_complemented( s2 ) << std::endl;

    if ( id1 == id2 ) {
      return;
    }

    auto rep1 = get_eqrepr( n1 );
    auto rep2 = get_eqrepr( n2 );

    if ( ntk.node_to_index( rep1 ) == ntk.node_to_index( rep2 ) ) {
      /* already in the same class */
      return;
    }

    /* merge the eq lists, set the representative as the node with lowest ID */
    if ( ntk.node_to_index( rep1 ) < ntk.node_to_index( rep2 ) )
    {
      /* before merging, complement nodes accordingly to the new representative phase, if needed */
      bool inv = false;
      if ( ( ( eqrep[id1] != n1 && ntk.is_complemented( eqnodes[id1] ) ) != ntk.is_complemented( s2 ) ) !=
             ( eqrep[id2] != n2 && ntk.is_complemented( eqnodes[id2] ) ) )
      {
        inv_eqnodes( rep2 );
        inv = true;
      }
      eqrep[ntk.node_to_index( rep2 )] = ntk.get_node( eqnodes[ntk.node_to_index( rep1 )] );
      eqnodes[ntk.node_to_index( rep1 )] = eqnodes[ntk.node_to_index( rep2 )];
      /* store the right phase */
      eqnodes[ntk.node_to_index( rep2 )] = ntk.make_signal( rep2 ) ^ inv;
    }
    else
    {
      bool inv = false;
      if ( ( ( eqrep[id1] != n1 && ntk.is_complemented( eqnodes[id1] ) ) != ntk.is_complemented( s2 ) ) !=
             ( eqrep[id2] != n2 && ntk.is_complemented( eqnodes[id2] ) ) )
      {
        inv_eqnodes( rep1 );
        inv = true;
      }
      eqrep[ntk.node_to_index( rep1 )] = ntk.get_node( eqnodes[ntk.node_to_index( rep2 )] );
      eqnodes[ntk.node_to_index( rep2 )] = eqnodes[ntk.node_to_index( rep1 )];
      eqnodes[ntk.node_to_index( rep1 )] = ntk.make_signal( rep1 ) ^ inv;
    }
  }


  node get_eqrepr( node const& n ) const
  {
    assert( ntk.node_to_index( n ) < ntk.size() );

    auto rep = eqrep[ntk.node_to_index( n )];
    while ( rep != eqrep[ntk.node_to_index( rep )] ) {
      rep = eqrep[ntk.node_to_index( rep )];
    }
    return rep;
  }


  std::vector<node> get_eqnodes( node const& n ) const
  {
    assert( ntk.node_to_index( n ) < ntk.size() );

    std::vector<node> eqnd;
    auto p = n;

    while ( ntk.node_to_index( p ) != eqrep[ntk.node_to_index( p )] ) {
      p = eqrep[ntk.node_to_index( p )];
      eqnd.push_back( p );
    }
    p = ntk.get_node( eqnodes[ntk.node_to_index( p )] );
    while ( ntk.node_to_index( p ) != ntk.node_to_index( n ) ) {
      eqnd.push_back( p );
      p = eqrep[ntk.node_to_index( p )];
    }
    return eqnd;
  }


  bool is_eqrepr( node const& n ) const
  {
    return eqrep[ntk.node_to_index( n )] == n;
  }


  signal get_eqrepr_signal( node const& n ) const
  {
    auto repr = ntk.make_signal( get_eqrepr( n ) );

    if ( ntk.get_node( repr ) == n ) {
      return repr;
    }

    return repr ^ ntk.is_complemented( eqnodes[ntk.node_to_index( n )] );
  }


  signal get_eqrepr_signal( signal const& sig ) const
  {
    auto n = ntk.get_node( sig );
    auto repr = get_eqrepr( n );

    if ( repr == n ) {
      return sig;
    }

    bool c = ntk.is_complemented( eqnodes[ntk.node_to_index( n )] ) != ntk.is_complemented( sig );
    return ntk.make_signal( repr ) ^ c;
  }


  void print_eqclasses( std::ostream& os = std::cout ) const
  {
    ntk.foreach_gate( [&]( auto const& n ) {
      if ( eqrep[ntk.node_to_index( n )] == n ) {
        auto p = ntk.get_node( eqnodes[ntk.node_to_index( n )] );
        if ( p == n ) {
          return true;
        }
        while ( ntk.node_to_index( p ) != ntk.node_to_index( n ) ) {
          if ( !ntk.is_dead( p ) ) {
            std::cout << "fail: ";
          }
          os << p << "(" << ntk.is_complemented( eqnodes[p] ) << ") ";
          p = eqrep[ntk.node_to_index( p )];
        }
        assert( !ntk.is_dead( p ) );
        os << p << std::endl;
      }
      return true;
    } );
  }


  template<typename Fn>
  void foreach_node_in_eqclass( node const& n, Fn&& fn ) const
  {
    auto p = n;
    if ( !fn( p ) ) {
      return;
    }
    while ( ntk.node_to_index( p ) != eqrep[ntk.node_to_index( p )] ) {
      p = eqrep[ntk.node_to_index( p )];
      if ( !fn( p ) ) {
        return;
      }
    }
    p = ntk.get_node( eqnodes[ntk.node_to_index( p )] );
    while ( ntk.node_to_index( p ) != ntk.node_to_index( n ) ) {
      if ( !fn( p ) ) {
        return;
      }
      p = eqrep[ntk.node_to_index( p )];
    }
  }


private:
  void init_eqclass()
  {
    eqrep.reserve( ntk.size() );
    eqnodes.reserve( ntk.size() );
    ntk.foreach_node( [&]( auto n ) {
      eqrep[ntk.node_to_index( n )] = n;
      eqnodes[ntk.node_to_index( n )] = ntk.make_signal( n );
    } );
  }


  void inv_eqnodes( node const& rep )
  {
    assert( ntk.node_to_index( rep ) < ntk.size() );

    auto p = ntk.get_node( eqnodes[rep] );

    while ( ntk.node_to_index( p ) != eqrep[ntk.node_to_index( p )] ) {
      eqnodes[p] = !eqnodes[p];
      p = eqrep[ntk.node_to_index( p )];
    }
  }

private:
  Ntk const& ntk;
  std::vector<node> eqrep;
  std::vector<signal> eqnodes;
};


template<class T>
eq_classes(T const&) -> eq_classes<T>;

} // namespace mockturtle
