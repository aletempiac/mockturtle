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
  double rise_block_delay;
  double fall_block_delay;
};

struct sGate
{
    std::string name;
    uint32_t id; 
    bool is_super{false};
    //gate root{};
    int32_t root_id{-1};
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
    explicit superGate_library( std::vector<mockturtle::gate> const& gates, mockturtle::super_info const& val = {} , std::vector<mockturtle::map_superGate> const& vec_sg = {} )
        : _gates( gates ),
          _val( val ),
          _vec_sg( vec_sg ),
          sg_list( )

    {
        if ( _vec_sg.size() == 0 )
            compute_library_with_genlib();
        else
            generate_library_with_super();
    }

    const std::vector<mockturtle::sGate> get_sg_library() const
    {
        return sg_list;
    }

public:
    void compute_library_with_genlib()
    {
        //assert( _vec_sg.size() == 0 );
        //assert( _gates.size() != 0 );

        for(const auto &g: _gates)
        {
            sGate s;
            s.name      = g.name;
            s.id        = sg_list.size();
            s.num_vars  = g.num_vars;
            s.is_super  = false;
            s.root_id      = g.id;
            s.area      = g.area;
            s.function  = g.function;
            for( auto const &p: g.pins)
            {
                s.pins.emplace_back( sGate_pin{ p.rise_block_delay, 
                        p.fall_block_delay} );
            }
            sg_list.emplace_back( s ); 
        }

    }

    void  generate_library_with_super()  
    {
        assert( _vec_sg.size() != 0 );
        assert( _gates.size() != 0 );

        /* Creating elementary gates */ 
        for(uint8_t i = 0; i < _val.max_num_vars; ++i)
        {
            sGate s;
            s.name      = "elementary_" + std::to_string( i );
            s.id        = sg_list.size();
            s.is_super  = false;
            s.num_vars  = _val.max_num_vars;
            kitty::dynamic_truth_table tt{_val.max_num_vars};

            kitty::create_nth_var( tt, i );
            s.function = tt;

            std::vector<mockturtle::sGate_pin> pp;

            for (uint32_t k = 0; k < _val.max_num_vars; ++k)
            {
                pp.emplace_back( sGate_pin{ 0.0f, 0.0f} );
            }
            s.pins = pp;
            sg_list.emplace_back( s ); 
        }



        for(auto const v:_vec_sg)
        {
            sGate s;
            s.id = sg_list.size();
            s.is_super = v.is_super; 

            bool match_found = false;
            for (auto const& g: _gates )
            {
                if( v.name == g.name )
                {
                    s.root_id = g.id;
                    match_found = true;
                    break;
                }

            }
            if( !match_found )
            {
                std::cout << "Some issue parsing the .super file" << std::endl;
            }

            s.num_vars = _gates[s.root_id].num_vars;
            s.name     = _gates[s.root_id].name + "_super_" + std::to_string( s.id );

            if ( s.num_vars != v.fanins_id.size() )
                std::cout << "num_vars != fanins_id.size " << std::endl;
            if (s.num_vars > _val.max_num_vars )
                std::cout << "num_vars cannot be more than max " << std::endl;

            /* Should not have more than the entries in super file */
            if (sg_list.size( ) > _val.num_lines)
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
            compute_delay_parameters( s );

            sg_list.emplace_back( s );

        }
    }

private:
    void compute_truth_table( mockturtle::sGate& s ) 
    {
        if ( !s.is_super )
        {
            s.function = _gates[s.root_id].function; 
        }
        else
        {
            std::vector<kitty::dynamic_truth_table> ttv;
            for (auto const leaf: s.fanins)
            {
                auto tt_leaf = leaf.function;
                ttv.emplace_back(  tt_leaf ) ; 
            }

            auto func = kitty::compose_truth_table( _gates[s.root_id].function, ttv );
            const auto support = kitty::min_base_inplace( func );
            s.function = kitty::shrink_to( func, static_cast<unsigned int>( support.size() ) );
        }
    }

    void compute_delay_parameters( mockturtle::sGate& s)
    {
        std::vector<mockturtle::sGate_pin> pp;
        const auto& root = _gates[s.root_id];

        /* setting initial delay */
        for (uint32_t k = 0; k < _val.max_num_vars; ++k)
        {
            pp.emplace_back( sGate_pin{ 0.0f, 0.0f} );
        }
        s.pins = pp;

        uint32_t i = 0; // For tracking the fanins of the current supergate s

        /* adding fanin delays */
        for(const auto& p: root.pins)
        {
            auto leaf = s.fanins[i];
            auto rise_block_delay = p.rise_block_delay;
            auto fall_block_delay = p.fall_block_delay;
            
            for(auto k = 0u; k < _val.max_num_vars; k++)
            {
                if ( leaf.pins[k].rise_block_delay >= 0 ) 
                {
                    if ( s.pins[k].rise_block_delay < leaf.pins[k].rise_block_delay + rise_block_delay )
                        s.pins[k].rise_block_delay = leaf.pins[k].rise_block_delay + rise_block_delay;
                }
                if ( leaf.pins[k].fall_block_delay >= 0 ) 
                {
                    if ( s.pins[k].fall_block_delay < leaf.pins[k].fall_block_delay + fall_block_delay )
                        s.pins[k].fall_block_delay = leaf.pins[k].fall_block_delay + fall_block_delay;
                }

            }
            ++i;
        }
    }

    void compute_area( mockturtle::sGate& s )
    {
        s.area = _gates[s.root_id].area;
        for( const auto leaf: s.fanins)
        {
            s.area += leaf.area; 
        }

    }

protected:
    std::vector<gate> const _gates; 
    mockturtle::super_info const _val;
    std::vector<mockturtle::map_superGate> const _vec_sg;
    std::vector<sGate> sg_list;

}; /* Class superGate_library */

} /* namespace mockturtle */

