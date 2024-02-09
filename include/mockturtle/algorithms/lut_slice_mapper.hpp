/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file lut_slice_mapper.hpp
  \brief LUT slice mapper

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <string>

#include <fmt/format.h>

#include "../networks/klut.hpp"
#include "../utils/cost_functions.hpp"
#include "../utils/cuts.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/truth_table_cache.hpp"
#include "../views/mffc_view.hpp"
#include "../views/topo_view.hpp"
#include "cleanup.hpp"
#include "collapse_mapped.hpp"
#include "cut_enumeration.hpp"
#include "exorcism.hpp"
#include "simulation.hpp"

namespace mockturtle
{

/*! \brief Parameters for slice mapper.
 *
 * The data structure `lut_map_slice_params` holds configurable parameters
 * with default arguments for `lut_slice_map`.
 */
struct lut_slice_map_params
{
  /*! \brief Slice mapping policy */
  enum
  {
    labeling,
    cut_based
  } policy = labeling;

  /*! \brief Required delay */
  uint32_t required_delay{ 0u };

  /*! \brief Propagation delay of an LUT. */
  float lut_delay{ 0 };

  /*! \brief Routing delay between LUTs in a slice. */
  float intra_slice_delay{ 0 };

  /*! \brief Routing delay between LUTs in different slices. */
  float inter_slice_delay{ 1 };

  /*! \brief LUTs in slice. */
  uint32_t slice_size{ 8 };

  uint32_t num_priority_cuts{ 8 };

  // /*! \brief Required depth relaxation ratio (%). */
  // uint32_t relax_required{ 0u };

  /*! \brief Number of rounds for area sharing optimization. */
  uint32_t area_share_rounds{ 2u };

  /*! \brief Number of rounds for area flow optimization. */
  uint32_t area_flow_rounds{ 1u };

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t ela_rounds{ 2u };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for slice mapper.
 *
 * The data structure `lut_slice_map_stats` provides data
 * collected by running `lut_slice_map`.
 */
struct lut_slice_map_stats
{
  /*! \brief Area result. */
  uint32_t area{ 0 };
  /*! \brief Worst delay result. */
  float delay{ 0 };
  /*! \brief Slices. */
  uint32_t slices{ 0 };

  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Cut enumeration stats. */
  uint32_t num_cuts{ 0 };

  /*! \brief Delay and area stats for each round. */
  std::vector<std::string> round_stats{};

  void report() const
  {
    for ( auto const& stat : round_stats )
    {
      std::cout << stat;
    }
    std::cout << fmt::format( "[i] Total runtime           = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

#pragma region LUT mapper

class lut_slice_map_impl
{
private:
  template<uint32_t max_size>
  struct slice
  {
    std::array<uint32_t, max_size> luts;
    uint8_t size{ 0 };
    float delay{ 0 };
  };

  template<uint32_t max_size, typename T>
  struct slice_set
  {
    std::array<T, max_size> luts;
    uint8_t size{ 0 };
  };

public:
  static constexpr uint32_t max_slice_num = 32;
  static constexpr uint32_t max_slice_size = 32;
  using slice_t = slice<max_slice_size>;
  using slice_set_t = slice_set<max_slice_num, slice_t>;

public:
  explicit lut_slice_map_impl( klut_network const& ntk, lut_slice_map_params const& ps, lut_slice_map_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    ntk.clear_values();

    /* TODO: implement cut-based approach */
    run_labeling();
  }

private:
  void run_labeling()
  {
    std::vector<uint32_t> label_size;
    label_size.reserve( ntk.size() );

    std::vector<float> node_delay( ntk.size() );

    uint32_t slice_id = 0;

    ntk.foreach_node( [&]( auto const& n ) {
      if ( n < ntk.num_pis() + 2 ) //if ( ntk.is_pi( f ) || ntk.is_constant( f ) )
      {
        node_delay[n] = 0;
        return;
      }

      uint32_t index = 0;
      uint32_t num_critical = 0;
      float d = 0;
      ntk.foreach_fanin( n, [&]( auto f, uint32_t i ) {
        if ( f < ntk.num_pis() + 2 ) //if ( ntk.is_pi( f ) || ntk.is_constant( f ) )
        {
          return;
        }

        if ( node_delay[f] > d )
        {
          d = node_delay[f];
          index = i;
          num_critical = 0;
        }
        else if ( node_delay[f] == d )
        {
          /* multiple critical inputs, slices don't help */
          ++num_critical;
        }
      } );

      uint32_t fanin_id = ntk.value( index );
      if ( num_critical == 0 && d != 0 && label_size[fanin_id] < ps.slice_size )
      {
        /* assign slice to the critical fanin */
        ntk.set_value( n, fanin_id );
        label_size[fanin_id]++;
        node_delay[n] = d + ps.intra_slice_delay;
      }
      else
      {
        /* add a new slice */
        ntk.set_value( n, slice_id++ );
        label_size.push_back( 1 );
        node_delay[n] = 0;
        index = ntk.fanin_size( n );
      }

      /* compute node delay */
      ntk.foreach_fanin( n, [&]( auto f, uint32_t i ) {
        if ( i == index )
          return;
        node_delay[n] = std::max( node_delay[n], node_delay[f] + ps.inter_slice_delay );
      } );
    } );

    /* compute worst delay */
    float worst_delay = 0;
    ntk.foreach_po( [&]( auto const& f ) {
      worst_delay = std::max( worst_delay, node_delay[f] );
    } );

    st.round_stats.push_back( fmt::format( "[i] Labeling : Delay = {:>8.2f}  Area = {:8d}  Slices = {:8d}\n", worst_delay, ntk.num_gates(), slice_id ) );

    st.area = ntk.num_gates();
    st.delay = worst_delay;
    st.slices = slice_id;
  }

private:
  klut_network const& ntk;
  lut_slice_map_params const& ps;
  lut_slice_map_stats& st;

  // uint32_t iteration{ 0 };       /* current mapping iteration */
  // uint32_t area_iteration{ 0 };  /* current area iteration */
  // uint32_t delay{ 0 };           /* current delay of the mapping */
  // uint32_t area{ 0 };            /* current area of the mapping */
  // uint32_t edges{ 0 };           /* current edges of the mapping */
  // uint32_t cuts_total{ 0 };      /* current computed cuts */
  // const float epsilon{ 0.005f }; /* epsilon */
  // LUTCostFn lut_cost{};

  // std::vector<node> topo_order;
  // std::vector<uint32_t> tmp_visited;
  // std::vector<node_lut> node_match;

  // std::vector<cut_set_t> cuts;  /* compressed representation of cuts */
};
#pragma endregion

} /* namespace detail */

/*! \brief Map LUT networks to slices.
 *
 */
void lut_slice_map( klut_network const& ntk, lut_slice_map_params const& ps = {}, lut_slice_map_stats* pst = nullptr )
{
  lut_slice_map_stats st;

  detail::lut_slice_map_impl p( ntk, ps, st );
  p.run();

  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */