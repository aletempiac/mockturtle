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

  \author Shubham Rai 
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
#include "../io/super_reader.hpp"
#include "../traits.hpp"
#include "../utils/truth_table_cache.hpp"

namespace mockturtle
{

struct sGate_pin
{
  double rise_tdelays;
  double fall_tdelays;
  double max_tdelay;
};

struct sGate
{
    uint32_t id; 
    bool is_super{false};
    struct gate const* root;
    uint32_t num_vars;
    kitty::dynamic_truth_table function;
    double area{0.0f};
    std::vector<sGate_pin> pins{};
    std::vector<sGate> fanins{}; 
};

template<unsigned NInputs = 5u>
class superGate_library
{
public:
    explicit superGate_library( std::vector<mockturtle::gate> const& gates, mockturtle::super_info const& val, std::vector<mockturtle::map_superGate> const& vec_sg )
        : _gates( gates ),
          val( val ),
          vec_sg( vec_sg ),
          sg_lib( )

    {
        generate_library_with_super();
    }

    const std::vector<mockturtle::sGate>* get_sg_library()
    {
        return &sg_lib;
    }

public:
    void  generate_library_with_super()  
    {
        assert( vec_sg.size() != 0 );
        assert( _gates.size() != 0 );

        std::vector<sGate> sg_list;

        /* Creating a hash table of the genlib gates accessible by their names */
        std::unordered_map<std::string, gate> gate_by_name;
        for(auto const g: _gates)
        {
            gate_by_name[g.name] = g;
        }

        /* Creating elementary gates */ 
        for(uint8_t i = 0; i < val.max_num_vars; ++i)
        {
            sGate s;
            s.id        = sg_list.size();
            s.num_vars  = val.max_num_vars;
            s.root      = nullptr;
            kitty::dynamic_truth_table tt{val.max_num_vars};

            kitty::create_nth_var( tt, i );
            s.function = tt;
            std::cout << "s.function ";
            kitty::print_binary( s.function );
            std::cout<<std::endl;

            std::vector<mockturtle::sGate_pin> pp;

            for (uint32_t k = 0; k < val.max_num_vars; ++k)
            {
                pp.emplace_back( sGate_pin{ 0.0f, 0.0f, 0.0f} );
            }
            s.pins = pp;
            sg_list.emplace_back( s ); 
        }

        std::cout << "Size of sg_list " << sg_list.size() << std::endl;


        for(auto const v:vec_sg)
        {
            sGate s;
            s.id = sg_list.size();
            s.is_super = v.is_super; 

            auto g = gate_by_name.find( v.name );
            if ( g != gate_by_name.end() )
                s.root = &g->second;

            s.num_vars = s.root->num_vars;

            if ( s.num_vars != v.fanins_id.size() )
                std::cout << "num_vars != fanins_id.size " << std::endl;

            /* Should not have more than the entries in super file */
            if (sg_list.size( ) > val.num_lines)
            {
                std::cout << "[i] the number of supergates exceed the number of lines in .super file" << std::endl;
            }

            /* Creating Fanins of the gate */ 
            for(auto const f: v.fanins_id)
            {
                if (f < 0)
                {
                    std::cerr << "[i] Supergate entry without any fanins " << std::endl;
                }
                s.fanins.emplace_back( sg_list[f] );
            }

            compute_truth_table( s );
            compute_area( s );

            std::cout << "s.function ";
            kitty::print_binary( s.function );
            std::cout << "\t area " << s.area;
            std::cout<<std::endl;

            sg_list.emplace_back( s );

        }
    }

private:
    void compute_truth_table( mockturtle::sGate& s ) 
    {
        auto root_gate = s.root; 
        if ( !s.is_super )
        {
            s.function = s.root->function; 
        }
        else
        {
            std::vector<kitty::dynamic_truth_table> ttv;
            std::vector<kitty::dynamic_truth_table> a( val.max_num_vars, kitty::dynamic_truth_table( val.max_num_vars ) );
            for (auto const leaf: s.fanins)
            {
                auto tt_leaf = leaf.function;
                ttv.emplace_back( kitty::extend_to( tt_leaf, NInputs ) ); 
            }

            s.function = kitty::compose_truth_table( s.root->function, ttv );
        }
    }

    void compute_delay_parameters( mockturtle::sGate& s)
    {

    }

    void compute_area( mockturtle::sGate& s )
    {
        auto root_area = s.area = s.root->area;
        for( auto const leaf: s.fanins)
        {
            s.area += leaf.area; 
        }

    }

protected:
    std::vector<gate> const _gates; 
    mockturtle::super_info val;
    std::vector<mockturtle::map_superGate> vec_sg;
    std::vector<mockturtle::sGate> sg_lib;

}; /* Class superGate_library */

} /* namespace mockturtle */


