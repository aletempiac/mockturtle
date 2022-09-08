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
  \file write_genlib.hpp
  \brief Write library of gates to GENLIB format

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <fmt/format.h>

#include "genlib_reader.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace mockturtle
{

/*! \brief Writes a library of gates in GENLIB format into output stream
 *
 * \param gates List of gates
 * \param os Output stream
 */
void write_genlib( std::vector<gate> const& gates, std::ostream& os )
{
  for ( auto const& g : gates )
  {
    os << "GATE ";
    std::string name = g.name;
    std::string expression = g.expression;
    const uint32_t name_spacing = std::max( 29u, static_cast<uint32_t>( name.length() + 1 ) );

    os << name.append( std::string( name_spacing - name.length(), ' ' ) );
    os << fmt::format( "{:>5.0f} ", g.area );

    const uint32_t expr_spacing = std::max( 29u, static_cast<uint32_t>( expression.length() + 2 ) );
    os << g.output_name << "=" << expression << ";";

    if ( g.pins.size() != 0 )
    {
      os << std::string( expr_spacing - expression.length(), ' ' );
    }

    for ( auto const& p : g.pins )
    {
      std::string phase;
      if ( p.phase == phase_type::INV )
        phase = "INV";
      else if ( p.phase == phase_type::NONINV )
        phase = "NONINV";
      else
        phase = "UNKNOWN";

      os << fmt::format( "PIN {} {} {} {} {} {} {} {}  ", p.name, phase, p.input_load, p.max_load, p.rise_block_delay, p.rise_fanout_delay, p.fall_block_delay, p.fall_fanout_delay );
    }
    os << "\n";
  }

  os << std::flush;
}

/*! \brief Writes a library of gates in GENLIB format into a file
 *
 * \param gates List of gates
 * \param filename Filename
 */
void write_genlib( std::vector<gate> const& gates, std::string const& filename )
{
  std::ofstream os( filename.c_str(), std::ofstream::out );
  write_genlib( gates, os );
  os.close();
}

} /* namespace mockturtle */