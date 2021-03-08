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

#include "../io/genlib_reader.hpp"
#include "../traits.hpp"

namespace mockturtle
{

// struct gate
// {
//   std::string name;
//   std::string expression;
//   double area;
//   double delay;
//   kitty::dynamic_truth_table function;

//   uint8_t n_inputs;

//   gate() : area( 0.0f ), delay( 0.0f ) {}
// };

template<unsigned NInputs>
struct supergate
{
  struct gate *root;

  /* area */
  float area;
  /* worst delay */
  float worstDelay;
  /* permuted pin-to-pin delay w.r.t NP-repr. */
  float tdelay[NInputs];

  /* permutation vector to np representative*/
  std::array<uint8_t, NInputs> permutation;
  /* permutated negation */
  uint8_t polarity{0};
};


template<unsigned NInputs = 4u>
class tech_library
{

public:
  tech_library( std::vector<gate> gates )
    : inv_area( 0 ),
      inv_delay( 0 ),
      _gates( gates ),
      _super_lib()
  {
    generate_library();
  }

  const std::vector<supergate<NInputs>>* get_supergates( kitty::static_truth_table<NInputs> const& tt ) const
  {
    auto match = _super_lib.find( tt );
    if ( match != _super_lib.end() )
      return &match->second;
    return NULL;
  }

  const std::tuple<float, float> get_inverter_info() const
  {
    return std::make_pair( inv_area, inv_delay );
  }

private:
  void generate_library()
  {
    for ( auto& gate : _gates )
    {
      if ( gate.function.num_vars() == 1 )
      {
        if ( kitty::is_const0( kitty::cofactor1( gate.function, 0 ) ) )
        {
          inv_area = gate.area;
          inv_delay = gate.delay;
        }
      }

      /* NPN canonization of the function */
      const auto tt = kitty::extend_to<NInputs>( gate.function );
      std::cout << gate.name << ": " << std::endl;
      kitty::print_hex( tt );
      std::cout << std::endl;
      auto [tt_np, neg, perm] = kitty::exact_npn_canonization( tt );
      /* from NPN class to NP */
      if ( ( ( neg >> NInputs ) & 1 ) == 1 )
      {
        tt_np = ~tt_np;
        neg ^= 1 << NInputs;
      }
      kitty::print_hex( tt_np );
      std::cout << std::endl;
      std::cout << std::endl;
      /* initialize gate info */
      supergate<NInputs> sg;
      sg.root = &gate;
      sg.area = gate.area;
      sg.worstDelay = gate.delay;
      sg.polarity = 0;
      /* a permutation index points to the new input at that position */
      for ( auto i = 0u; i < perm.size() && i < NInputs; ++i )
      {
        sg.permutation[perm[i]] = i;
        sg.tdelay[perm[i]] = gate.delay;
        sg.polarity |= ( ( neg >> perm[i] ) & 1 ) << i;
      }
      _super_lib[tt_np].push_back( sg );
    }

    printf( "inv(%.2f, %.2f)\n", inv_delay, inv_area );
    for ( auto const& entry : _super_lib )
    {
      kitty::print_hex( entry.first );
      std::cout << ": ";
      for ( auto const& gate : entry.second )
      {
        printf( "%s(%.2f, %.2f) ", gate.root->name.c_str(), gate.worstDelay, gate.area );
      }
      std::cout << std::endl;
    }
  }

private:
  float inv_area;
  float inv_delay;
  std::vector<gate> _gates;
  std::unordered_map<kitty::static_truth_table<NInputs>, std::vector<supergate<NInputs>>, kitty::hash<kitty::static_truth_table<NInputs>>> _super_lib;
};


template<typename Ntk, unsigned NInputs>
struct exact_supergate
{
  signal<Ntk> const root;

  /* number of inputs of the supergate */
  unsigned n_inputs : 3;
  /* saved polarities for inputs and/or outputs */
  unsigned polarity : 8;
  
  /* area */
  float area;
  /* worst delay */
  float worstDelay;
  /* pin-to-pin delay */
  float tdelay[NInputs];

  exact_supergate( signal<Ntk> const root )
  : root( root ),
    n_inputs( 0u ),
    polarity( 0u ),
    area( 0.0f ),
    worstDelay( 0.0f )
    {
      for ( auto i = 0u; i < NInputs; i++ )
      {
        tdelay[i] = 0.0f;
      }
    }
};

struct exact_library_params
{
  float area_gate{1.0f};
  float area_inverter{0.2f};
  float delay_gate{1.0f};
  float delay_inverter{1.0f};

  bool np_classification{true};
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

  const std::tuple<float, float> get_inverter_info() const
  {
    return std::make_pair( ps.area_inverter, ps.delay_inverter );
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
      std::vector<exact_supergate<Ntk, NInputs>> supergates_pos;
      std::vector<exact_supergate<Ntk, NInputs>> supergates_neg;
      auto const not_entry = ~entry;

      const auto add_supergate = [&]( auto const& f_new ) {
        bool complemented = database.is_complemented( f_new );
        auto f = f_new;
        if ( ps.np_classification && complemented ) {
          f = !f;
        }
        exact_supergate<Ntk, NInputs> sg( f );
        compute_info( sg );
        if ( ps.np_classification && complemented )
        {
          supergates_neg.push_back( sg );
        }
        else
        {
          supergates_pos.push_back( sg );
        }
        database.create_po( f );
        return true;
      };

      kitty::dynamic_truth_table function = kitty::extend_to( entry, NInputs );
      rewriting_fn( database, function, pis.begin(), pis.end(), add_supergate );
      if ( supergates_pos.size() > 0 )
        super_lib.insert( {entry, supergates_pos} );
      if ( ps.np_classification && supergates_neg.size() > 0 )
        super_lib.insert( {not_entry, supergates_neg} );
    }

    if ( ps.verbose )
    {
      std::cout << "Classified in " << super_lib.size() << " entries" << std::endl;
      for ( auto const &pair : super_lib )
      {
        kitty::print_hex( pair.first );
        std::cout << ": ";

        for ( auto const&  gate : pair.second )
        {
          printf( "%.2f,%.2f,%d ", gate.worstDelay, gate.area, gate.n_inputs );
        }
        std::cout << std::endl;
      }
    }
  }

  /* Computes delay and area info */
  void compute_info( exact_supergate<Ntk, NInputs> &sg )
  {
    database.incr_trav_id();
    /* info does not consider input and output inverters */
    bool compl_root = database.is_complemented( sg.root );
    auto const root = compl_root ? !sg.root : sg.root;
    sg.area = compute_info_rec( sg, root, 0.0f );

    /* output polarity */
    sg.polarity |= ( unsigned( compl_root ) ) << NInputs;
    /* number of inputs */
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

    if ( database.is_pi( n ) )
    {
      sg.tdelay[database.index_to_node( n ) - 1u] = std::min(sg.tdelay[database.index_to_node( n ) - 1u], tdelay);
      sg.worstDelay = std::min(sg.worstDelay, tdelay);
      sg.polarity |= ( unsigned( database.is_complemented( root ) ) ) << ( database.index_to_node( n ) - 1u );
      return area;
    }

    tdelay -= ps.delay_gate;

    /* add gate area once */
    if ( database.visited( n ) != database.trav_id() )
    {
      area += ps.area_gate;
      database.set_value( n, 0u );
      database.set_visited( n, database.trav_id() );
    }

    if ( database.is_complemented( root ) )
    {
      tdelay -= ps.delay_inverter;
      /* add inverter area only once (shared by fanout) */
      if ( database.value( n ) == 0u )
      {
        area += ps.area_inverter;
        database.set_value( n, 1u );
      }
    }

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
