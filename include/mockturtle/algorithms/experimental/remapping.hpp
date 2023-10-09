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
  \file rremapping.hpp
  \brief A rremapping engine

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <chrono>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "../../networks/block.hpp"
#include "../../networks/klut.hpp"
#include "../../utils/cuts.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"
#include "../../utils/tech_library.hpp"
#include "../../views/binding_view.hpp"
#include "../../views/cell_view.hpp"
#include "../../views/topo_view.hpp"
#include "../cleanup.hpp"
#include "../cut_enumeration.hpp"
#include "../detail/mffc_utils.hpp"
#include "../detail/switching_activity.hpp"

namespace mockturtle
{

/*! \brief Parameters for remap.
 *
 * The data structure `remap_params` holds configurable parameters
 * with default arguments for `remap`.
 */
struct remap_params
{
  /*! \brief Do area-oriented remapping. */
  bool area_oriented_remapping{ false };

  /*! \brief Required time for delay optimization. */
  double required_time{ 0.0f };

  /*! \brief Maps using multi-output gates */
  bool use_multioutput{ false };

  /*! \brief Window number of PIs */
  uint32_t num_pis{ 12 };

  /*! \brief Window number of POs */
  uint32_t num_pos{ 12 };

  remap_params()
  {
    cut_enumeration_ps.cut_limit = 16;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 16.
   * The maximum cut limit is 15.
   * By default, truth table minimization
   * is performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for remap.
 *
 * The data structure `remap_stats` provides data collected by running
 * `remap`.
 */
struct remap_stats
{
  /*! \brief Area result. */
  double area{ 0 };
  /*! \brief Worst delay result. */
  double delay{ 0 };
  /*! \brief Power result. */
  double power{ 0 };
  /*! \brief Power result. */
  uint32_t inverters{ 0 };

  /*! \brief Mapped multi-output gates. */
  uint32_t multioutput_gates{ 0 };

  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Rempping error. */
  bool remapping_error{ false };

  void report() const
  {
    std::cout << fmt::format( "[i] Area = {:>5.2f}; Delay = {:>5.2f};", area, delay );
    if ( power != 0 )
      std::cout << fmt::format( " Power = {:>5.2f};\n", power );
    else
      std::cout << "\n";
    if ( multioutput_gates )
    {
      std::cout << fmt::format( "[i] Multi-output gates   = {:>5}\n", multioutput_gates );
      std::cout << fmt::format( "[i] Multi-output runtime = {:>5.2f} secs\n", to_seconds( time_multioutput ) );
    }
    std::cout << fmt::format( "[i] Total runtime        = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

template<class Ntk, typename WindowEngine, typename ResynFn, typename MapperFn, unsigned CutSize, unsigned NInputs, classification_type Configuration>
class remap_impl
{
private:
public:
  explicit remap_impl( Ntk& ntk, tech_library<NInputs, Configuration> const& library, remap_params const& ps = {}, remap_stats& pst )
      : ntk( ntk ),
        library( library ),
        ps( ps ),
        pst( pst )
  {}

  void run()
  {
    /* TODO: implement */
  }

private:
  Ntk& ntk;
  // WindowsEngine engine;
  // ResynFn resyn;
  tech_library<NInputs, Configuration> const& library; /* TODO: replace with &&Lib? */
  remap_params& ps;
  remap_stats& st;
};

} /* namespace detail */

/*! \brief Remapping.
 *
 * This function implements a simple technology mapping algorithm.
 * The algorithm maps each node to the first implementation in the technology library.
 *
 * The input must be a binding_view with the gates correctly loaded.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `has_binding`
 *
 * \param ntk Network
 *
 */
template<class Ntk, typename WindowEngine, typename ResynFn, typename MapperFn, unsigned CutSize = 6u, unsigned NInputs, classification_type Configuration>
void remap( Ntk& ntk, tech_library<NInputs, Configuration> const& library, remap_params const& ps = {}, remap_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_has_binding_v<Ntk> || has_has_cell_v<Ntk>, "Ntk does not implement the has_binding or has_cell method" );

  remap_stats st;
  detail::remap_impl<Ntk, CutSize, NInputs, Configuration> p( ntk, library, ps, st );
  auto res = p.run();

  if ( ps.verbose && !st.remapping_error )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */