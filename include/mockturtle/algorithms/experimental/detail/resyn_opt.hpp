/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file resyn_opt.hpp
  \brief Optimization scripts for resynthesis

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <iostream>
#include <optional>
#include <vector>

#include "../../../utils/algorithm.hpp"
#include "../../../utils/tech_library.hpp"
#include "../../rewrite.hpp"
#include "../../sim_resub.hpp"
#include "../../node_resynthesis/xag_npn.hpp"

namespace mockturtle
{

namespace detail
{

template<class Ntk>
struct resyn_aig_size
{
  resyn_aig_size()
      : lib()
  {
    xag_npn_resynthesis<Ntk, Ntk, xag_npn_db_kind::aig_complete> resyn;
    lib.add_library( resyn );
  }

  void operator()( Ntk& ntk ) const
  {
    /* run resub */
    resubstitution_params ps;
    resubstitution_stats st;

    ps.max_inserts = 20;
    ps.max_pis = 12;
    ps.max_divisors = std::numeric_limits<uint32_t>::max();
    sim_resubstitution( ntk );
    ntk = cleanup_dangling( ntk );

    /* run rewriting */
    rewrite( ntk, lib );
  }

private:
  exact_library<Ntk> lib;
};

} /* namespace detail */

} /* namespace mockturtle */