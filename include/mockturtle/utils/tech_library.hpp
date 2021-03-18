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

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>

#include "../io/genlib_reader.hpp"
#include "../traits.hpp"

namespace mockturtle
{

struct tech_library_params
{
  /*! \brief reports np enumerations */
  bool verbose{false};

  /*! \brief reports all the entries in the library */
  bool very_verbose{false};
};


template<unsigned NInputs>
struct supergate
{
  struct gate const* root;

  /* area */
  float area;
  /* worst delay */
  float worstDelay;
  /* pin-to-pin delay */
  std::array<float, NInputs> tdelay;

  /* np permutation vector */
  std::vector<uint8_t> permutation;

  /* pin negations */
  uint8_t polarity{0};
};


template<unsigned NInputs = 5u>
class tech_library
{
  using lib_t = std::unordered_map<kitty::static_truth_table<NInputs>, std::vector<supergate<NInputs>>, kitty::hash<kitty::static_truth_table<NInputs>>>;

public:
  tech_library( std::vector<gate> const gates, tech_library_params const ps = {} )
    : _gates( gates ),
      _ps ( ps ),
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

  const std::tuple<float, float, uint32_t> get_inverter_info() const
  {
    return std::make_tuple( inv_area, inv_delay, inv_id );
  }

  unsigned max_gate_size()
  {
    return max_size;
  }

  const std::vector<gate> get_gates() const
  {
      return _gates;
  }

private:
  void generate_library()
  {
    for ( auto& gate : _gates )
    {
      if ( gate.function.num_vars() == 1 )
      {
        /* extract inverter delay and area */
        if ( kitty::is_const0( kitty::cofactor1( gate.function, 0 ) ) )
        {
          inv_area = gate.area;
          inv_delay = gate.delay;
          inv_id = gate.id;
        }
      }
      if ( gate.function.num_vars() > NInputs )
      {
        std::cerr << "WARNING: gate " << gate.name << " IGNORED, too many variables for the library settings" << std::endl;
        continue;
      }

      max_size = std::max( max_size, gate.num_vars );

      uint32_t np_count = 0;

      const auto on_np = [&]( auto const& tt, auto neg, auto const& perm ) {
        supergate<NInputs> sg;
        sg.root = &gate;
        sg.area = gate.area;
        sg.worstDelay = gate.delay;
        sg.polarity = 0;
        sg.permutation = perm;

        for ( auto i = 0u; i < perm.size() && i < NInputs; ++i )
        {
          sg.tdelay[i] = gate.delay;  /* if pin-to-pin delay change to: gate.delay[perm[i]] */
          sg.polarity |= ( ( neg >> perm[i] ) & 1 ) << i; /* permutate input negation to match the right pin */
        }
        for ( auto i = perm.size(); i < NInputs; ++i )
        {
          sg.tdelay[i] = 0; /* added for completeness but not necessary */
        }

        const auto static_tt = kitty::extend_to<NInputs>( tt );

        auto& v = _super_lib[static_tt];

        /* ordered insert by ascending area and number of input pins */
        auto it = std::lower_bound( v.begin(), v.end(), sg, [&]( auto const& s1, auto const& s2 ) {
          if ( s1.area < s2.area )
            return true;
          if ( s1.area > s2.area )
            return false;
          if ( s1.root->num_vars < s2.root->num_vars )
            return true;
          if ( s1.root->num_vars > s2.root->num_vars )
            return true;
          return s1.root->id < s2.root->id;
        } );

        bool to_add = true;
        /* search for duplicated element due to symmetries */
        while ( it != v.end() )
        {
          if ( sg.root->id == it->root->id )
          {
            /* if already in the library exit, else ignore permutations if with equal delay cost */
            if ( sg.polarity == it->polarity && sg.tdelay == it->tdelay )
            {
              to_add = false;
              break;
            }
          }
          else
          {
            break;
          }
          ++it;
        }

        if ( to_add )
        {
          v.insert( it, sg );
          ++np_count;
        }

        /* check correct results */
        // assert( gate.function == create_from_npn_config( std::make_tuple( tt, neg, sg.permutation ) ) );
      };

      /* NP enumeration of the function */
      const auto tt = gate.function;
      kitty::exact_np_enumeration( tt, on_np );

      if ( _ps.verbose )
      {
        std::cout << "Gate " << gate.name << ", num_vars = " << gate.num_vars << ", np entries = " << np_count << std::endl;
      }
    }

    if ( _ps.very_verbose )
    {
      for ( auto const& entry : _super_lib )
      {
        kitty::print_hex( entry.first );
        std::cout << ": ";
        for ( auto const& gate : entry.second )
        {
          printf( "%s(d:%.2f, a:%.2f, p:%d) ", gate.root->name.c_str(), gate.worstDelay, gate.area, gate.polarity );
        }
        std::cout << std::endl;
      }
    }
  }

private:
  /* inverter info */
  float inv_area{0.0};
  float inv_delay{0.0};
  uint32_t inv_id{0};

  unsigned max_size{0}; /* max #fanins of the gates in the library */

  std::vector<gate> const _gates; /* collection of gates */
  tech_library_params const _ps;
  lib_t _super_lib; /* library of enumerated gates */
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
  float area_gate{2.5f};
  float area_inverter{1.0f};
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
          printf( "%.2f,%.2f,%d,%d,:", gate.worstDelay, gate.area, gate.polarity, gate.n_inputs );
          for ( auto j = 0u; j < NInputs; ++j )
            printf( "%.2f/", gate.tdelay[j] );
          std::cout << " ";
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
