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
  \file tech_library.hpp
  \brief Implements utilities to enumerates gates for technology mapping

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
#include <lorina/lorina.hpp>

#include "../io/genlib_reader.hpp"
#include "../utils/super.hpp"
#include "../traits.hpp"
#include "../utils/truth_table_cache.hpp"

namespace mockturtle
{

/*
std::string mcnc_library =  "GATE   inv1    1	O=!a;		        PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                            "#GATE  inv2	  2	O=!a;		        PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                            "#GATE  inv3	  3	O=!a;		        PIN * INV 3 999 1.1 0.09 1.1 0.09\n"
                            "#GATE  inv4	  4	O=!a;		        PIN * INV 4 999 1.2 0.07 1.2 0.07\n"
                            "GATE   nand2	  2	O=!(ab);		    PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                            "GATE   nand3	  3	O=!(abc);	      PIN * INV 1 999 1.1 0.3 1.1 0.3\n"
                            "GATE   nand4   4	O=!(abcd);	    PIN * INV 1 999 1.4 0.4 1.4 0.4\n"
                            "GATE   nor2	  2	O=!{ab};		    PIN * INV 1 999 1.4 0.5 1.4 0.5\n"
                            "GATE   nor3	  3	O=!{abc};	      PIN * INV 1 999 2.4 0.7 2.4 0.7\n"
                            "GATE   nor4	  4	O=!{abcd};	    PIN * INV 1 999 3.8 1.0 3.8 1.0\n"
                            "GATE   and2	  3	O=(ab);		      PIN * NONINV 1 999 1.9 0.3 1.9 0.3\n"
                            "GATE   or2		  3	O={ab};		      PIN * NONINV 1 999 2.4 0.3 2.4 0.3\n"
                            "GATE   xor2a	  5	O=[ab];     	  PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                            "#GATE  xor2b	  5	O=[ab];     	  PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                            "GATE   xnor2a	5	O=![ab];		    PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                            "#GATE  xnor2b	5	O=![ab];		    PIN * UNKNOWN 2 999 2.1 0.5 2.1 0.5\n"
                            "GATE   aoi21	  3	O=!{(ab)c};	    PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                            "GATE   aoi22	  4	O=!{(ab)(cd)};	PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                            "GATE   oai21	  3	O=!({ab}c);	    PIN * INV 1 999 1.6 0.4 1.6 0.4\n"
                            "GATE   oai22	  4	O=!({ab}{cd});	PIN * INV 1 999 2.0 0.4 2.0 0.4\n"
                            "GATE   buf    	2	O=a;        	  PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                            "GATE   zero	  0	O=0;\n"
                            "GATE   one		  0	O=1;";
*/

struct tech_library_params
{
  /*! \brief reports np enumerations */
  bool verbose{false};

  /*! \brief allows computing supergates */
  bool compute_supergates{false};

  /*! \brief Pruning with respect to the fanin of the root_gate while computing supergates */
  uint32_t prune_limit{3};

  /*! \brief Max delay which the combinational supergates can have */
  float max_delay{0.0f};

  /*! brief no of levels for computing supergates */
  bool levels{1};

  /* brief no of fanins for the supergates */
  uint8_t supergates_fanins{5};

  /* brief count of total supergates */
  uint32_t total_supergates{10000};

  /*! \brief reports all the entries in the library */
  bool very_verbose{false};
};

template<unsigned NInputs>
struct comb_supergate
{
    uint32_t id;
    uint32_t root_id;
    uint32_t num_vars{0};
    bool is_comb_supergate{false};
    double area{0.0f};
    /* worst delay */
    float worst_delay{0.0f};
    ///* pin-to-pin delay */
    std::array<float, NInputs> tdelay;
    std::vector<int32_t> fanin_list{};
    kitty::dynamic_truth_table function;
};



template<unsigned NInputs>
struct supergate
{
  struct gate const* root;

  /* pointer to the comb_supergate */
  struct comb_supergate<NInputs> const* root_cg;

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
  using supergates_list_t = std::vector<supergate<NInputs>>;
  using tt_hash = kitty::hash<kitty::static_truth_table<NInputs>>;
  using lib_t = std::unordered_map<kitty::static_truth_table<NInputs>, supergates_list_t, tt_hash>;
  using list_supergate_t = std::vector<comb_supergate<NInputs>>;

public:
  explicit tech_library( std::vector<gate> const& gates, std::vector<mockturtle::map_superGate>vec_superGates, mockturtle::super_info const& vals, tech_library_params const ps = {} )
    : _gates( gates ),
      _superGates( vec_superGates ),
      _ps ( ps ),
      _lsg(),
      _super_lib()
  {
    if( vec_superGates.size() != 0 )
    {
        superGate_library<NInputs> sg_lib( gates, vals, vec_superGates );
    }
    else
        generate_library();
  }

  const supergates_list_t* get_supergates( kitty::static_truth_table<NInputs> const& tt ) const
  {
    auto match = _super_lib.find( tt );
    if ( match != _super_lib.end() )
      return &match->second;
    return nullptr;
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

  void compute_supergates( uint32_t levels = 1 )
  {
      assert( levels != 0 );

      
      /* TODO Sort the gates in order of increasing delay */
      

       
      /* computing level0 combination supergates */
			for (auto const& g: _gates)
      {
          auto nfanins = g.num_vars; 
          auto index = g.id; 
          auto root_area = g.area;
          auto root_delay = g.compute_worst_delay();

          if ( nfanins < 1 )
              continue;
          /* Root gates with more than prune_limit fanins would not be considered */
          if ( nfanins >= _ps.prune_limit )
              continue;

          /* have to update the csg_limit with filtering methods 
           * Only those gates which are possible with the current
           * root gates shall be considered 
           */
          list_supergate_t csg_list;
          for(auto l: _lsg)
          {
              if ( l.num_vars >= _ps.prune_limit )
              {
                  continue;
              }
              if( l.worst_delay + root_delay > _ps.max_delay )
              {
                  csg_list.emplace_back( l );
              }

          }
          std::cout << "CSG size " << csg_list.size() << "\t lsg size " << _lsg.size() << std::endl;
          switch( nfanins )
          {
              case 0: // Constants 
                  assert( 0 );
                  break;
              case 1: // Inverters 
                  for ( auto const& g0: csg_list )
                  {
                      /* skipping constants  and inverters */
                      if ( g0.num_vars <= 1)
                          continue;
                      if( compute_worst_delay( g0 ) == 0)
                          continue;

                      comb_supergate<NInputs> c0;
                      c0.num_vars = nfanins + g0.num_vars - 1;
                      if ( c0.num_vars > NInputs )
                          continue;

                      c0.root_id            = index;
                      c0.id                 = _lsg.size();
                      c0.is_comb_supergate  = true;
                      c0.fanin_list.emplace_back( g0.id );
                      c0.worst_delay        = compute_worst_delay( c0 ); 
                      c0.area               = compute_area ( c0 );
                      c0.function           = compute_tt(c0);

                      if ( !is_dominated( c0 ) )
                          continue;
                      _lsg.emplace_back( g0 );
                  }
                  break;
              case 2:  /* Two input gates */
                  /* when constants or the input variable are considered as one of the inputs, 
                   * then following considerations are made
                   * 0th variable will be given -1, 1th variable will be given -2 and so on 
                   * and so forth*/
                  for (auto const& g0: csg_list)
                  {
                      if ( g0.num_vars < 1)
                          continue;
                      comb_supergate<NInputs> c0;
                      c0.num_vars = nfanins + g0.num_vars - 1;
                      if ( c0.num_vars > NInputs )
                          continue;
                      c0.fanin_list.emplace_back( g0.id );
                      c0.root_id = index;
                      c0.id = _lsg.size();
                      c0.is_comb_supergate = true;
                      c0.fanin_list.emplace_back( -1 );
                      c0.worst_delay        = compute_worst_delay( c0 ); 
                      c0.area               = compute_area( c0 );
                      c0.function           = compute_tt(c0);

                      if ( !is_dominated( c0 ) )
                          continue;
                      _lsg.emplace_back( c0 );

                      for (auto const& g1: csg_list)
                      {
                          if ( g1.num_vars < 1)
                              continue;
                          comb_supergate<NInputs> c1;
                          c1.num_vars = nfanins + g0.num_vars + g1.num_vars - 2;
                          if ( c1.num_vars > NInputs )
                              break;
                          c1.fanin_list.emplace_back( g0.id );
                          c1.fanin_list.emplace_back( g1.id );
                          c1.root_id = index;
                          c1.id = _lsg.size();
                          c1.is_comb_supergate = true;
                          c1.worst_delay        = compute_worst_delay( c1 );
                          c1.area               = compute_area( c1 );
                          c1.function           = compute_tt( c1 );
                          if ( !is_dominated( c1 ) )
                              continue;
                          _lsg.emplace_back( c1 );
                      }

                  }
                  break;
              case 3: /* Three input gates */  
                  std::cout << "three input gates " << std::endl;
                  for (auto const& g0: _lsg)
                  {
                      if ( g0.num_vars < 1)
                          continue;
                      comb_supergate<NInputs> c0;
                      c0.num_vars = nfanins + g0.num_vars - 1;
                      if ( c0.num_vars > NInputs )
                          break;
                      c0.fanin_list.emplace_back( g0.id );
                      c0.root_id = index;
                      c0.id = _lsg.size();
                      c0.is_comb_supergate = true;

                      c0.fanin_list.emplace_back( -1 );
                      c0.fanin_list.emplace_back( -2 );
                      _lsg.emplace_back( c0 );

                      for (auto const& g1: _lsg)
                      {
                          if ( g1.num_vars < 1)
                              continue;
                          comb_supergate<NInputs> c1;
                          c1.num_vars = nfanins + g0.num_vars + g1.num_vars - 2;
                          if( c1.num_vars >  NInputs ) 
                              break;
                          c1.fanin_list.emplace_back( g0.id );
                          c1.fanin_list.emplace_back( g1.id );
                          c1.fanin_list.emplace_back( -1 );
                          c1.root_id = index; 
                          c1.id = _lsg.size();
                          c1.is_comb_supergate = true;
                          _lsg.emplace_back( c1 );

                          for (auto const& g2: _lsg)
                          {
                              if ( g2.num_vars < 1)
                                  continue;
                              comb_supergate<NInputs> c2;
                              c2.num_vars = nfanins + g0.num_vars + g1.num_vars + g2.num_vars - 3;
                              if ( c2.num_vars > NInputs )
                                  break;
                              c2.fanin_list.emplace_back( g0.id );
                              c2.fanin_list.emplace_back( g1.id );
                              c2.fanin_list.emplace_back( g2.id );
                              c2.root_id = index; 
                              c2.id = _lsg.size();
                              c2.is_comb_supergate = true;
                              _lsg.emplace_back( c2 );
                          }
                      }

                  }
                  break;
              case 4: /* Four input gates */
                  std::cout << "four input gates " << std::endl;
                  for (auto const& g0: _lsg )
                  {
                      if ( g0.num_vars < 1)
                          continue;
                      comb_supergate<NInputs> c0;
                      c0.num_vars = nfanins + g0.num_vars - 1;
                      if ( c0.num_vars > NInputs )
                          break;
                      c0.fanin_list.emplace_back( g0.id );
                      c0.root_id = index;
                      c0.id = _lsg.size();
                      c0.is_comb_supergate = true;

                      c0.fanin_list.emplace_back( -1 );
                      c0.fanin_list.emplace_back( -2 );
                      c0.fanin_list.emplace_back( -3 );
                      _lsg.emplace_back( c0 );

                      for (auto const& g1: _lsg)
                      {
                          if ( g1.num_vars < 1)
                              continue;
                          comb_supergate<NInputs> c1;
                          c1.num_vars = nfanins + g0.num_vars + g1.num_vars - 2;
                          if ( c1.num_vars > NInputs )
                              break;
                          c1.fanin_list.emplace_back( g0.id );
                          c1.fanin_list.emplace_back( g1.id );
                          c1.fanin_list.emplace_back( -1 );
                          c1.fanin_list.emplace_back( -2 );
                          c1.root_id = index; 
                          c1.id = _lsg.size();
                          c1.is_comb_supergate = true;
                          _lsg.emplace_back( c1 );

                          for (auto const& g2: _lsg)
                          {
                              if ( g2.num_vars < 1)
                                  continue;
                              comb_supergate<NInputs> c2;
                              c2.num_vars = nfanins + g0.num_vars + g1.num_vars + g2.num_vars - 3;
                              if ( c2.num_vars > NInputs )
                                  break;
                              c2.fanin_list.emplace_back( g0.id );
                              c2.fanin_list.emplace_back( g1.id );
                              c2.fanin_list.emplace_back( g2.id );
                              c2.fanin_list.emplace_back( -1 );
                              c2.root_id = index; 
                              c2.id = _lsg.size();
                              c2.is_comb_supergate = true;
                              _lsg.emplace_back( c2 );

                              for (auto const& g3: _lsg)
                              {
                                  if( g3.num_vars < 1)
                                      continue;
                                  comb_supergate<NInputs> c3;
                                  c3.num_vars = nfanins + g0.num_vars + g1.num_vars + g2.num_vars + g3.num_vars - 4;
                                  if ( c3.num_vars > NInputs )
                                      break;
                                  c3.fanin_list.emplace_back( g0.id );
                                  c3.fanin_list.emplace_back( g1.id );
                                  c3.fanin_list.emplace_back( g2.id );
                                  c3.fanin_list.emplace_back( g3.id );
                                  c3.root_id = index; 
                                  c3.id = _lsg.size();
                                  c3.is_comb_supergate = true;
                                  _lsg.emplace_back( c3 );
                              }
                          }
                      }
                  }
                  break;
              default:
                  //assert( 0 );
                  break;
          }

      }

      if ( _ps.very_verbose )
      {
          std::cout << "size of supergate " << _lsg.size() << std::endl;

          for (auto const& t : _lsg)
          {
              //auto c = std::get<0> ( t );
              auto c =  t ;
              std::cout << "root id " << c.root_id 
                  << "\t id "      << c.id 
                  << "\t is_comb_supergate " << c.is_comb_supergate 
                  << "\t fanins " << c.num_vars 
                  << "\t tt " ;
              kitty::print_binary(c.function);
              std::cout <<std::endl;

              //auto v = std::get<1> ( t );
              auto v = c.fanin_list;
              std::cout << "Printing vector elements" << std::endl; 
              for (auto const& e: v)
              {
                  std::cout << e << " "; 
              }
              std::cout <<std::endl;
          }
          fflush(stdout);
      }
  }

  kitty::dynamic_truth_table compute_tt(comb_supergate<NInputs> const& cg)
  {
      auto root_gate = _gates[cg.root_id];
      if ( !cg.is_comb_supergate )
      {
          return root_gate.function;
      }
      else
      {
          std::vector<kitty::dynamic_truth_table> ttv;
          std::vector<kitty::dynamic_truth_table> a( NInputs, kitty::dynamic_truth_table( NInputs ) );
          for (auto const leaf: cg.fanin_list)
          {
              /* Case 1:  When one of the fanin of the comb_supergate is another gate primitive */
              if ( leaf >= 0 )
              {
                  auto gate = _lsg[leaf];
                  ttv.emplace_back( kitty::extend_to( gate.function, NInputs ) );
              }
              /* Case 2: When it is nth variable */
              else
              {
                  auto val =  std::abs( leaf );
                  kitty::create_nth_var( a[val], val );
                  ttv.emplace_back( kitty::extend_to( a[val], NInputs ) );
              }
          }
          return kitty::compose_truth_table( root_gate.function, ttv );
      }
  }

  bool is_dominated(comb_supergate<NInputs> & cg)
  {
      auto is_new_better = false;
      auto is_old_better = false;
      if ( kitty::is_const0( cg.function ) || kitty::is_const0( ~cg.function ) )
          return false;
      auto count = 0u;
      for (auto g:_lsg)
      {
          if ( g.function == cg.function )
          {
              /* Constant value 0.001 taken from ABC */
              if ( cg.area > g.area + 0.001 )
                  is_old_better = true;
              else 
                  is_new_better = true;

              if ( cg.worst_delay > g.worst_delay + 0.001 )
                  is_old_better = true;
              else 
                  is_new_better = true;
              
              if ( is_new_better )
              {
                  _lsg.erase( _lsg.begin() + count ); 
              }
              else if ( is_old_better )
                  return false;
              else
                  return false;

          }
          else 
              return true;
          ++count;

      }
      return true;
  }

  void generate_library( const std::string& filename = "" )
  {
    for (auto const g: _gates)
    {
        comb_supergate<NInputs> c;
        c.root_id = g.id;
        c.area = g.area;
        c.id = _lsg.size();
        c.is_comb_supergate = false;
        c.num_vars = g.num_vars;
        c.function = kitty::extend_to( g.function, NInputs );
        c.worst_delay = compute_worst_delay( c );
        int fanin_index = -1;
        for (uint8_t i = 0 ; i < g.num_vars; i++)
        {
            c.fanin_list.emplace_back( fanin_index );
            --fanin_index;
            c.tdelay[i] = c.worst_delay;
        }
        _lsg.emplace_back( c ); 
    }


    if( _ps.compute_supergates )
    {
        std::cout << "Computing Combinational supergates " << std::endl; 
        compute_supergates( );
    }

    //if (1)
    //    throw std::runtime_error( " Hello " );

    bool inv = false;
        
    for ( auto const& cg : _lsg )
    {
      auto const& gate = _gates[cg.root_id];

      if ( cg.num_vars > NInputs )
      {
        std::cerr << "WARNING: gate " << gate.name << " IGNORED, too many variables for the library settings" << std::endl;
        continue;
      }

      float worst_delay = cg.worst_delay;
      float gate_area = cg.area;

      if ( cg.num_vars == 1 )
      {
        /* extract inverter delay and area */
        if ( kitty::is_const0( kitty::cofactor1( gate.function, 0 ) ) )
        {
          /* get the smallest area inverter */
          if ( !inv || gate.area < inv_area )
          {
            inv_area = gate.area;
            inv_delay = worst_delay;
            inv_id = gate.id;
            inv = true;
          }
        }
      }

      max_size = std::max( max_size, cg.num_vars);

      uint32_t np_count = 0;

      const auto on_np = [&]( auto const& tt, auto neg, auto const& perm ) {
        supergate<NInputs> sg;
        sg.root = &gate;
        sg.root_cg = &cg;
        //sg.area = gate.area;
        sg.area = gate_area;
        sg.worstDelay = worst_delay;
        sg.polarity = 0;
        sg.permutation = perm;

        for ( auto i = 0u; i < perm.size() && i < NInputs; ++i )
        {
          sg.tdelay[i] = worst_delay;  /* if pin-to-pin delay change to: gate.delay[perm[i]] */
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
          if ( s1.root_cg->num_vars < s2.root_cg->num_vars )
            return true;
          if ( s1.root_cg->num_vars > s2.root_cg->num_vars )
            return true;
          return s1.root_cg->id < s2.root_cg->id;
        } );

        bool to_add = true;
        /* search for duplicated element due to symmetries */
        while ( it != v.end() )
        {
          if ( sg.root_cg->id == it->root_cg->id )
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

      ///* NP enumeration of the function */
      kitty::exact_np_enumeration( cg.function, on_np );

      if ( _ps.verbose )
      {
        std::cout << "Gate " << gate.name << ", num_vars = " << gate.num_vars << ", np entries = " << np_count << std::endl;
      }
    }

    if ( !inv )
    {
      std::cerr << "[i] WARNING: inverter gate has not been detected in the library" << std::endl;
    }

    std::cout << "Size of superlib " << _super_lib.size() << std::endl;
    if ( _ps.very_verbose )
    {
      for ( auto const& entry : _super_lib )
      {
        kitty::print_hex( entry.first );
        std::cout << ": ";
        for ( auto const& gate : entry.second )
        {
          if ( gate.root_cg->is_comb_supergate )
          {
          }
          else
              printf( "%s(d:%.2f, a:%.2f, p:%d) ", gate.root->name.c_str(), gate.worstDelay, gate.area, gate.polarity );
        }
        std::cout << std::endl;
      }
    }
  }

  float compute_area( comb_supergate<NInputs> const& cg)
  {
      auto root_gate = _gates[cg.root_id];

      float total_area = root_gate.area;
      if (cg.is_comb_supergate)
      {
          for (auto const v: cg.fanin_list)
          {
              if (v > 0)
                  total_area += _lsg[v].area;
          }
          return total_area;
      }
      else 
          return root_gate.area; 
  }

  float compute_worst_delay( comb_supergate<NInputs> const& cg )
  {
    float worst_delay = 0.0f;
    auto const& r = _gates[cg.root_id];

		if ( cg.is_comb_supergate )
    {
        /* finding the pin delay for the root gate*/
        float d0 = 0.0f;
        float leaf_delay = 0.0f;
        for ( auto const& pin : r.pins )
        {
            float worst_pin_delay = static_cast<float>( std::max( pin.rise_block_delay, pin.fall_block_delay ) );
            d0 = std::max( d0, worst_pin_delay );
        }

        /* finding the pin delay for the leaf nodes */
        for (auto const v: cg.fanin_list)
        {
            if(v >= 0)
            {
                auto leaf_root_gate = _gates[_lsg[v].root_id];
                leaf_delay = std::max( leaf_delay, _lsg[v].worst_delay);
                
                //auto leaf_gates = _lsg[v];
                if ( leaf_root_gate.num_vars >= 1 )
                {
                    for (auto const& pin: leaf_root_gate.pins )
                    {
                        float worst_pin_delay = static_cast<float>( std::max( pin.rise_block_delay, pin.fall_block_delay ) );
                        leaf_delay = std::max(leaf_delay, worst_pin_delay);
                    }
                }
            }

        }
        return (d0 + leaf_delay);
		}
    else 
    {
        /* consider only block_delay */
        for ( auto const& pin : r.pins )
        {
            float worst_pin_delay = static_cast<float>( std::max( pin.rise_block_delay, pin.fall_block_delay ) );
            worst_delay = std::max( worst_delay, worst_pin_delay );
        }
    }
    return worst_delay;
  }

  void compute_tt( comb_supergate<NInputs> const& cg, std::vector<kitty::dynamic_truth_table>& res_tt )
  {
      auto root_gate = _gates[cg.root_id];

      if ( !cg.is_comb_supergate )
      {
          res_tt.emplace_back( root_gate.function );
      }
      else
      {
          std::vector<kitty::dynamic_truth_table> ttv;
          std::vector<std::vector<kitty::dynamic_truth_table>> multi_ttv;
          std::vector<kitty::dynamic_truth_table> a( NInputs, kitty::dynamic_truth_table( NInputs ) );

          auto pos = 0u;
          for (auto const leaf: cg.fanin_list)
          {
              /* Case 1:  When one of the fanin of the comb_supergate is another gate primitive */
              if ( leaf >= 0 )
              {
                  ++pos;
                  auto gate = _gates[leaf];
                  if( gate.num_vars > root_gate.num_vars )
                      ttv.emplace_back( kitty::shrink_to( gate.function, root_gate.num_vars ) );
                  else 
                      ttv.emplace_back( kitty::extend_to( gate.function, root_gate.num_vars ) );

              }
              /* Case 2: When it is nth variable */
              else
              {
                  auto val =  std::abs( leaf );
                  kitty::create_nth_var( a[val], val );
                  if( a[val].num_vars() > root_gate.num_vars )
                      ttv.emplace_back( kitty::shrink_to( a[val], root_gate.num_vars ) );
                  else 
                      ttv.emplace_back( kitty::extend_to( a[val], root_gate.num_vars ) );
              }
              /* Computing all combinations of intermediary negations */
          }

          generate_all_combinations(pos, multi_ttv, ttv, 0);

          for (auto i: multi_ttv)
          {
              res_tt.emplace_back( kitty::compose_truth_table( root_gate.function, i ) );
          }
      }
  }


  void generate_all_combinations(uint32_t const& pos, std::vector<std::vector<kitty::dynamic_truth_table>>& multi_ttv, std::vector<kitty::dynamic_truth_table>& ttv, uint32_t i)
  {
      std::vector<kitty::dynamic_truth_table> tt1;
      tt1 = ttv;
      if (i >= pos)
      {
          multi_ttv.emplace_back(ttv);
      }
      else
      {
          tt1[i] = ttv[i];
          generate_all_combinations(pos, multi_ttv, tt1, i + 1);

          tt1[i] = ~ttv[i];
          generate_all_combinations(pos, multi_ttv, tt1, i + 1);
      }
  }


private:
  /* inverter info */
  float inv_area{ 0.0 };
  float inv_delay{ 0.0 };
  uint32_t inv_id{ UINT32_MAX };

  unsigned max_size{0}; /* max #fanins of the gates in the library */

  std::vector<gate> const _gates; /* collection of gates */
  std::vector<map_superGate> const _superGates;
  tech_library_params const _ps;
  list_supergate_t _lsg;
  lib_t _super_lib; /* library of enumerated gates */
};


template<typename Ntk, unsigned NInputs>
struct exact_supergate
{
  signal<Ntk> const root;

  /* number of inputs of the supergate */
  uint8_t n_inputs{ 0 };
  /* saved polarities for inputs and/or outputs */
  uint8_t polarity{ 0 };

  /* area */
  float area{ 0 };
  /* worst delay */
  float worstDelay{ 0 };
  /* pin-to-pin delay */
  std::array<float, NInputs> tdelay{ 0 };

  exact_supergate( signal<Ntk> const root )
      : root( root ) {}
};

struct exact_library_params
{
  /* area of a gate */
  float area_gate{ 1.0f };
  /* area of an inverter */
  float area_inverter{ 0.0f };
  /* delay of a gate */
  float delay_gate{ 1.0f };
  /* delay of an inverter */
  float delay_inverter{ 0.0f };

  /* classify in NP instead of NPN */
  bool np_classification{ true };
  /* verbose */
  bool verbose{ false };
};

/*! \brief Library of exact synthesis supergates
 *
 * This class creates a technology library from an exact
 * synthesis database. Each NPN-entry in the database is
 * stored in its NP class by removing the output inverter
 * if present. The class creates supergates from the
 * database computing area and delay information.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      mockturtle::mig_npn_resynthesis mig_resyn{true};
      mockturtle::exact_library<mockturtle::mig_network, mockturtle::mig_npn_resynthesis, 4> lib( mig_resyn );
   \endverbatim
 */
template<typename Ntk, class RewritingFn, unsigned NInputs = 4u>
class exact_library
{
  using supergates_list_t = std::vector<exact_supergate<Ntk, NInputs>>;
  using tt_hash = kitty::hash<kitty::static_truth_table<NInputs>>;
  using lib_t = std::unordered_map<kitty::static_truth_table<NInputs>, supergates_list_t, tt_hash>;

public:
  explicit exact_library( RewritingFn const& rewriting_fn, exact_library_params const& ps = {} )
      : _database(),
        _rewriting_fn( rewriting_fn ),
        _ps( ps ),
        _super_lib()
  {
    generate_library();
  }

  const supergates_list_t* get_supergates( kitty::static_truth_table<NInputs> const& tt ) const
  {
    auto match = _super_lib.find( tt );
    if ( match != _super_lib.end() )
      return &match->second;
    return nullptr;
  }

  const Ntk& get_database() const
  {
    return _database;
  }

  const std::tuple<float, float> get_inverter_info() const
  {
    return std::make_pair( _ps.area_inverter, _ps.delay_inverter );
  }

private:
  void generate_library()
  {
    std::vector<signal<Ntk>> pis;
    for ( auto i = 0u; i < NInputs; ++i )
    {
      pis.push_back( _database.create_pi() );
    }

    /* Compute NPN classes */
    std::unordered_set<kitty::static_truth_table<NInputs>, tt_hash> classes;
    kitty::static_truth_table<NInputs> tt;
    do
    {
      const auto res = kitty::exact_npn_canonization( tt );
      classes.insert( std::get<0>( res ) );
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );

    /* Constuct supergates */
    for ( auto const& entry : classes )
    {
      supergates_list_t supergates_pos;
      supergates_list_t supergates_neg;
      auto const not_entry = ~entry;

      const auto add_supergate = [&]( auto const& f_new ) {
        bool complemented = _database.is_complemented( f_new );
        auto f = f_new;
        if ( _ps.np_classification && complemented )
        {
          f = !f;
        }
        exact_supergate<Ntk, NInputs> sg( f );
        compute_info( sg );
        if ( _ps.np_classification && complemented )
        {
          supergates_neg.push_back( sg );
        }
        else
        {
          supergates_pos.push_back( sg );
        }
        _database.create_po( f );
        return true;
      };

      kitty::dynamic_truth_table function = kitty::extend_to( entry, NInputs );
      _rewriting_fn( _database, function, pis.begin(), pis.end(), add_supergate );
      if ( supergates_pos.size() > 0 )
        _super_lib.insert( { entry, supergates_pos } );
      if ( _ps.np_classification && supergates_neg.size() > 0 )
        _super_lib.insert( { not_entry, supergates_neg } );
    }

    if ( _ps.verbose )
    {
      std::cout << "Classified in " << _super_lib.size() << " entries" << std::endl;
      for ( auto const& pair : _super_lib )
      {
        kitty::print_hex( pair.first );
        std::cout << ": ";

        for ( auto const& gate : pair.second )
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
  void compute_info( exact_supergate<Ntk, NInputs>& sg )
  {
    _database.incr_trav_id();
    /* info does not consider input and output inverters */
    bool compl_root = _database.is_complemented( sg.root );
    auto const root = compl_root ? !sg.root : sg.root;
    sg.area = compute_info_rec( sg, root, 0.0f );

    /* output polarity */
    sg.polarity |= ( unsigned( compl_root ) ) << NInputs;
    /* number of inputs */
    for ( auto i = 0u; i < NInputs; ++i )
    {
      sg.tdelay[i] *= -1; /* invert to positive value */
      if ( sg.tdelay[i] != 0.0f )
        sg.n_inputs++;
    }
    sg.worstDelay *= -1;
  }

  float compute_info_rec( exact_supergate<Ntk, NInputs>& sg, signal<Ntk> const& root, float delay )
  {
    auto n = _database.get_node( root );

    if ( _database.is_constant( n ) )
      return 0.0f;

    float area = 0.0f;
    float tdelay = delay;

    if ( _database.is_pi( n ) )
    {
      sg.tdelay[_database.index_to_node( n ) - 1u] = std::min( sg.tdelay[_database.index_to_node( n ) - 1u], tdelay );
      sg.worstDelay = std::min( sg.worstDelay, tdelay );
      sg.polarity |= ( unsigned( _database.is_complemented( root ) ) ) << ( _database.index_to_node( n ) - 1u );
      return area;
    }

    tdelay -= _ps.delay_gate;

    /* add gate area once */
    if ( _database.visited( n ) != _database.trav_id() )
    {
      area += _ps.area_gate;
      _database.set_value( n, 0u );
      _database.set_visited( n, _database.trav_id() );
    }

    if ( _database.is_complemented( root ) )
    {
      tdelay -= _ps.delay_inverter;
      /* add inverter area only once (shared by fanout) */
      if ( _database.value( n ) == 0u )
      {
        area += _ps.area_inverter;
        _database.set_value( n, 1u );
      }
    }

    _database.foreach_fanin( n, [&]( auto const& child ) {
      area += compute_info_rec( sg, child, tdelay );
    } );

    return area;
  }

private:
  Ntk _database;
  RewritingFn const& _rewriting_fn;
  exact_library_params const _ps;
  lib_t _super_lib;
};

} // namespace mockturtle
