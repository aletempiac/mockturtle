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
  \file gates.hpp
  \brief Implements utilities to manage and create gates

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cassert>
#include <vector>
#include <unordered_map>
#include <optional>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>

#include "../traits.hpp"

namespace mockturtle::detail
{

struct gate
{
  std::string name;
  double area;
  double delay;
  kitty::dynamic_truth_table tt;

  uint8_t n_inputs;

  gate() : area( 0.0f ), delay( 0.0f ) {}
};


struct supergate
{
  struct gate *root;
  kitty::dynamic_truth_table tt;

  unsigned n_inputs : 3;
  
  float area;
  float delay;

  supergate() : area( 0.0f ), delay( 0.0f ) {}
};


template<typename Ntk, unsigned NInputs>
struct exact_supergate
{
  signal<Ntk> const root;
  kitty::static_truth_table<NInputs> const& tt;

  unsigned n_inputs : 3;
  
  float area;
  float worstDelay;
  float tdelay[NInputs];

  exact_supergate( signal<Ntk> const root, kitty::static_truth_table<NInputs> const& tt )
  : root( root ),
    tt( tt ),
    n_inputs( 0u ),
    area( 0.0f ),
    worstDelay( 0.0f )
    {
      for ( auto i = 0u; i < NInputs; i++ )
      {
        tdelay[i] = 0.0f;
      }
    }
};


class library
{
  library()
  {
  }

};


struct exact_library_params
{
  float area_gate{1.0f};
  float area_inverter{0.0f};
  float delay_gate{1.0f};
  float delay_inverter{0.0f};

  bool verbose{false};
};


template<typename Ntk, class RewritingFn, unsigned NInputs>
class exact_library
{

public:
  exact_library( RewritingFn const& rewriting_fn, exact_library_params const& ps = {} )
  : database(),
    rewriting_fn( rewriting_fn ),
    ps( ps ),
    super_lib()

  {
    generate_library();
  }

  const std::vector<exact_supergate<Ntk, NInputs>>* get_supergates( kitty::static_truth_table<NInputs> const& tt ) const
  {
    auto match = super_lib.find( tt );
    if ( match != super_lib.end() )
      return &match->second;
    return NULL;
  }

  const Ntk &get_database() const
  {
    return database;
  }


private:
  void generate_library()
  {
    std::vector<signal<Ntk>> pis;
    for ( auto i = 0u; i < NInputs; ++i )
    {
      pis.push_back( database.create_pi() );
    }

    /* Compute NPN classes */
    std::unordered_set<kitty::static_truth_table<NInputs>, kitty::hash<kitty::static_truth_table<NInputs>>> classes;
    kitty::static_truth_table<NInputs> tt;
    do
    {
      const auto res = kitty::exact_npn_canonization( tt );
      classes.insert( std::get<0>( res ) );
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );

    /* Constuct supergates */
    for ( auto const &entry : classes )
    {
      std::vector<exact_supergate<Ntk, NInputs>> supergates;

      const auto add_supergate = [&]( auto const& f_new ) {
        exact_supergate<Ntk, NInputs> sg( f_new, entry );
        compute_info( sg );
        supergates.push_back( sg );
        database.create_po( f_new );
        return true;
      };

      kitty::dynamic_truth_table function = kitty::extend_to( entry, NInputs );
      rewriting_fn( database, function, pis.begin(), pis.end(), add_supergate );
      super_lib.insert( {entry, supergates} );
    }

    if ( ps.verbose )
    {
      for ( auto const &entry : classes )
      {
        kitty::print_hex( entry );
        std::cout << ": ";

        for ( auto const&  gate : super_lib[entry] )
        {
          printf( "%.0f,%.0f,%d ", gate.worstDelay, gate.area, gate.n_inputs );
        }
        std::cout << std::endl;
      }
    }
  }


  void compute_info( exact_supergate<Ntk, NInputs> &sg )
  {
    sg.area = compute_info_rec( sg, sg.root, 0.0f );
    for( auto i = 0u; i < NInputs; ++i )
    {
      if ( sg.tdelay[i] != 0.0f )
        sg.n_inputs++;
    }
  }


  float compute_info_rec( exact_supergate<Ntk, NInputs> &sg, signal<Ntk> const& root, float delay )
  {
    auto n = database.get_node( root );

    if ( database.is_constant( n ) )
      return 0.0f;

    float area = 0.0f;
    float tdelay = delay;

    if ( database.is_complemented( root ) )
    {
      area += ps.area_inverter;
      tdelay -= ps.delay_inverter;
    }

    if ( database.is_pi( n ) )
    {
      sg.tdelay[database.index_to_node( n ) - 1u] = std::min(sg.tdelay[database.index_to_node( n ) - 1u], tdelay);
      sg.worstDelay = std::min(sg.worstDelay, tdelay);
      return area;
    }

    area += ps.area_gate;
    tdelay -= ps.delay_gate;

    database.foreach_fanin( n, [&]( auto const& child ) {
      area += compute_info_rec( sg, child, tdelay );
    } );

    return area;
  }

private:
  Ntk database;
  RewritingFn const& rewriting_fn;
  exact_library_params const& ps;
  std::unordered_map<kitty::static_truth_table<NInputs>, std::vector<exact_supergate<Ntk, NInputs>>, kitty::hash<kitty::static_truth_table<NInputs>>> super_lib;
};

}
