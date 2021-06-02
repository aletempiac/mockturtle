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

const std::vector<sGate>  generate_library_with_super( std::vector<gate> const& gates, std::vector<map_superGate> const& vec_superGates, mockturtle::super_info const& vals)
{
    assert( vec_superGates.size() != 0 );
    assert( gates.size() != 0 );

    std::vector<sGate> sg_list;

    /* Creating a hash table of the genlib gates accessible by their names */
    std::unordered_map<std::string, gate> gate_by_name;
    for(auto const g: gates)
    {
        gate_by_name[g.name] = g;
    }

    /* Creating elementary gates */ 
    for(uint8_t i = 0; i < vals.max_num_vars; ++i)
    {
        sGate s;
        s.id        = sg_list.size();
        s.num_vars  = vals.max_num_vars;
        s.root      = nullptr;
        kitty::dynamic_truth_table tt{vals.max_num_vars};

        kitty::create_nth_var( tt, i );
        s.function = tt;

        std::vector<mockturtle::sGate_pin> pp;

        for (uint32_t k = 0; k < vals.max_num_vars; ++k)
        {
            pp.emplace_back( sGate_pin{ 0.0f, 0.0f, 0.0f} );
        }
        s.pins = pp;
        sg_list.emplace_back( s ); 
    }

    std::cout << "Size of sg_list " << sg_list.size() << std::endl;


    for(auto const s:vec_superGates)
    {
        sGate s1;
        s1.id = sg_list.size();
        s1.is_super = s.is_super; 
        
        auto g = gate_by_name.find( s.name );
        if ( g != gate_by_name.end() )
            s1.root = &g->second;
        else
            return {}; 

        s1.num_vars = s1.root->num_vars;

        if ( s1.num_vars == s.fanins_id.size() )
            std::cout << "num_vars == fanins_id.size " << std::endl;
        else
            std::cout << "num_vars != fanins_id.size " << std::endl;

        /* Creating Fanins of the gate */ 
        for(auto const f: s.fanins_id)
        {
            if (f < 0)
            {
                std::cout << "Gate without any fanins " << std::endl;
                return {};
            }
            s1.fanins.emplace_back( sg_list[f] );
        }
        sg_list.emplace_back( s1 );

    }
    return sg_list;
}

} /* namespace mockturtle */


