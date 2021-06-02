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
  \file genlib_reader.hpp
  \brief Reader visitor for GENLIB files

  \author Heinz Riener
*/

#pragma once

#include "../traits.hpp"
#include "./genlib_reader.hpp"

#include <kitty/constructors.hpp>
#include <lorina/genlib.hpp>
#include <fmt/format.h>

namespace mockturtle
{
    
struct map_superGate
{
    unsigned int id;
    std::string name{};
    uint32_t num_vars{0};
    bool is_super{false};
    std::vector<uint32_t> fanins_id;
};

struct super_info
{
    std::string genlib_name{};
    uint32_t max_num_vars{0u};
    uint32_t num_superGates{0u};
    uint32_t num_lines{0};
};


class super_reader: public lorina::super_reader
{
public: 
    explicit super_reader( std::vector<map_superGate>& superGates, super_info& val )
        : superGates( superGates ),
        val( val )
    {}

    virtual void parse_super_info( std::string const& genlib_name, uint32_t const& max_num_vars,
          uint32_t const& max_superGates, uint32_t const& num_lines ) const override
    {
        val.genlib_name       = genlib_name;
        val.max_num_vars      = max_num_vars;
        val.num_superGates    = max_superGates; 
        val.num_lines         = num_lines;

        /* Emplacing elementary gates */ 
        for(uint8_t i = 0; i < val.max_num_vars; ++i)
        {
            map_superGate s;
            s.id       = superGates.size();
            s.num_vars = val.max_num_vars;

            for (uint32_t k = 0; k < val.max_num_vars; ++k)
            {
                s.fanins_id.emplace_back( k );
            }
            superGates.emplace_back( s );
        }
    }

    virtual void on_line ( std::string const& name, bool const& is_super, std::vector<uint32_t> const& fanins_id ) const override
    {
        map_superGate s;

        s.name      = name;
        s.id        = superGates.size(); 
        s.is_super  =  is_super;
        s.fanins_id = fanins_id;
        superGates.emplace_back( s ); 
    }

protected:
    std::vector<map_superGate>& superGates;
    super_info& val;
}; /* super_reader */

} /* namespace mockturtle */ 
