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

  \author Alessandro Tempia Calvino
  \author Shubham Rai
*/

#pragma once

#include "../traits.hpp"

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <lorina/super.hpp>

namespace mockturtle
{

struct supergate_spec
{
  unsigned int id;
  std::string name{};
  bool is_super{ false };
  std::vector<uint32_t> fanin_id;
};

struct super_lib
{
  std::string genlib_name{};
  uint32_t max_num_vars{ 0u };
  uint32_t num_supergates{ 0u };
  uint32_t num_lines{ 0 };
  std::vector<supergate_spec> supergates{};
};

/*! \brief lorina callbacks for SUPER files.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      std::vector<mockturtle::supergates_spec> supergates;
      lorina::read_genlib( "file.super", mockturtle::super_reader( supergates ) );
   \endverbatim
 */
class super_reader : public lorina::super_reader
{
public:
  explicit super_reader( super_lib& lib )
      : lib( lib )
  {
  }

  virtual void on_super_info( std::string const& genlib_name, uint32_t max_num_vars, uint32_t max_superGates, uint32_t num_lines ) const override
  {
    lib.genlib_name = genlib_name;
    lib.max_num_vars = max_num_vars;
    lib.num_supergates = max_superGates; 
    lib.num_lines = num_lines;
  }

  virtual void on_supergate( std::string const& name, bool const& is_super, std::vector<uint32_t> const& fanins_id ) const override
  {

    lib.supergates.emplace_back( supergate_spec{ lib.supergates.size(),
                                                 name,
                                                 is_super,
                                                 fanins_id } );
  }

protected:
  super_lib& lib;
}; /* super_reader */

} /* namespace mockturtle */
