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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/experimental/decompose_multioutput.hpp>
#include <mockturtle/algorithms/experimental/emap.hpp>
#include <mockturtle/algorithms/map_adders.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/dont_touch_view.hpp>

#include <experiments.hpp>

using namespace mockturtle;

template<class Ntk, class Library>
bool add_adders_binding_info( Ntk& ntk, Library const& lib )
{
  static_assert( has_is_dont_touch_v<Ntk>, "Ntk does not implement the is_dont_touch method" );
  static_assert( has_add_binding_v<Ntk>, "Ntk does not implement the add_binding method" );

  uint32_t ha_nand_id = 0, ha_xnor_id = 0, fa_min_id = 0, fa_xnor_id = 0;
  auto const& multi_gates = lib.get_multioutput_gates();

  /* get half adder and full adder gates */
  for ( auto const& mg : multi_gates )
  {
    for ( auto const& g : mg )
    {
      if ( g.function._bits[0] == 0x7 )
      {
        ha_nand_id = g.root->id;
      }
      else if ( g.function._bits[0] == 0x9 )
      {
        ha_xnor_id = g.root->id;
      }
      else if ( g.function._bits[0] == 0x17 )
      {
        fa_min_id = g.root->id;
      }
      else if ( g.function._bits[0] == 0x69 )
      {
        fa_xnor_id = g.root->id;
      }
    }
  }

  /* some IDs are not matched */
  if ( !ha_nand_id || !ha_xnor_id || !fa_min_id || !fa_xnor_id )
    return false;

  bool success = true;
  ntk.foreach_node( [&]( auto const& n ) {
    if ( !ntk.is_dont_touch( n ) )
      return;
    
    uint64_t tt = ntk.node_function( n )._bits[0];

    if ( tt == 0x7 )
    {
      ntk.add_binding( n, ha_nand_id );
    }
    else if ( tt == 0x9 )
    {
      ntk.add_binding( n, ha_xnor_id );
    }
    else if ( tt == 0x17 )
    {
      ntk.add_binding( n, fa_min_id );
    }
    else if ( tt == 0x69 )
    {
      ntk.add_binding( n, fa_xnor_id );
    }
    else
    {
      success = false;
    }
  } );

  return success;
}

int main()
{
  using namespace experiments;
  using block_dt_t = dont_touch_view<block_network>;

  experiment<std::string, uint32_t, double, double, uint32_t, double, double, uint32_t, uint32_t, float, float, bool, bool> exp(
      "map_adders", "benchmark", "size", "area_det", "area_emap", "depth", "delay_det", "delay_emap", "adders_det", "adders_emap", "runtime_det", "runtime_emap", "cec_det", "cec_emap" );

  /* library to map to technology */
  fmt::print( "[i] processing technology library\n" );
  std::vector<gate> gates;
  std::ifstream in( "asap7.genlib" );

  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
  {
    return 1;
  }

  tech_library_params tps;
  tps.verbose = true;
  tps.load_multioutput_gates = true;
  tps.load_multioutput_gates_single = true;
  tech_library<6, classification_type::np_configurations> tech_lib( gates, tps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    /* METHOD 1 :map adders in two steps: FA-HA detection followed by emap */
    map_adders_params ps_ma;
    ps_ma.map_inverted = true;
    map_adders_stats st_ma;
    block_network res_det = map_adders( aig, ps_ma, &st_ma );

    block_dt_t block_res = decompose_multioutput<block_network, block_dt_t>( res_det, { true } );
    binding_view<block_dt_t> partial_map_res{ block_res, gates };
    if ( !add_adders_binding_info( partial_map_res, tech_lib ) )
    {
      std::cout << "[e] Failed at adding adder binding info\n";
      return -1;
    }
    double initial_area = partial_map_res.compute_area();

    emap_params ps1;
    emap_stats st1;
    binding_view<klut_network> det_emap = emap<binding_view<block_dt_t>, 6>( partial_map_res, tech_lib, ps1, &st1 );
    bool const cec1 = ( benchmark == "hyp" ) ? true : abc_cec( det_emap, benchmark );
    st1.area -= initial_area / 2; /* area of multioutput gates is counted as double */


    /* METHOD 2: map adders in one step using emap */
    emap_params ps2;
    ps2.map_multioutput = true;
    emap_stats st2;
    binding_view<klut_network> res_emap = emap<aig_network, 6>( aig, tech_lib, ps2, &st2 );

    exp( benchmark, size_before, st1.area, st2.area, depth_before, st1.delay, st2.delay, st_ma.mapped_fa + st_ma.mapped_ha, st2.multioutput_gates, to_seconds( st_ma.time_total ) + to_seconds( st1.time_total ), to_seconds( st2.time_total ), cec1, true );
  }

  exp.save();
  exp.table();

  return 0;
}
