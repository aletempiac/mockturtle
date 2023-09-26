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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/det_randomization.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <mockturtle/utils/stopwatch.hpp>

#include <experiments.hpp>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, float> exp(
      "rewrite_comparison", "benchmark", "size", "size_mig", "depth", "depth_mig", "time_mig" );

  fmt::print( "[i] processing technology library\n" );

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  eps.np_classification = false;
  eps.enable_dont_cares = false;
  exact_library<mig_network, mig_npn_resynthesis> exact_lib( resyn, eps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );
    aig_network aig;
    if ( lorina::read_aiger( "optimized/" + benchmark + ".aig", aiger_reader( aig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = aig.num_gates();
    const uint32_t depth_before = depth_view( aig ).depth();

    map_params ps;
    map_stats st;
    
    mig_network mig;
    if ( lorina::read_aiger( "optimized/" + benchmark + ".aig", aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    float d_time = 0;
    float m_time = 0;

    ps.skip_delay_round = true;
    ps.required_time = std::numeric_limits<float>::max();
    ps.ela_rounds = 2;
    ps.enable_logic_sharing = true;
    ps.use_dont_cares = false;
    ps.window_size = 12;
    ps.logic_sharing_cut_limit = 1;
    ps.cut_enumeration_ps.cut_limit = 8;

    mig_network mig_opt = cleanup_dangling( mig );

    uint32_t i = 3;
    while ( i-- > 0 )
    {
      auto size_before = mig_opt.size();
      // auto depth_before = depth_view( mig_d ).depth();

      auto mig_map = map( mig_opt, exact_lib, ps, &st );

      d_time += to_seconds( st.time_total );
      m_time += to_seconds( st.time_matching );

      // if ( mig_map.size() >= size_before /*|| depth_after >= depth_before*/ )
      //   break;

      mig_opt = mig_map;
    }

    i = 3;
    while ( i-- > 0 )
    {
      auto size_before = mig_opt.size();
      // auto depth_before = depth_view( mig_d ).depth();

      rewrite_params rps;
      rps.use_dont_cares = i == 1;
      rps.allow_zero_gain = true;
      rps.odc_levels = 0;
      rps.window_size = 8;
      rewrite_stats rst;
      rewrite( mig_opt, exact_lib, rps, &rst );

      d_time += to_seconds( st.time_total );

      if ( mig_opt.size() >= size_before /*|| depth_after >= depth_before*/ )
        break;
    }

    while ( true )
    {
      auto size_global_before = mig_opt.size();
      /* resub */
      {
        resubstitution_params ps;
        resubstitution_stats st;
        ps.max_pis = 8u;
        ps.window_size = 8u;
        ps.max_inserts = 3u;
        // ps.progress = false;

        auto mig_resub = cleanup_dangling( mig_opt );

        depth_view depth_mig{mig_resub};
        fanout_view fanout_mig{depth_mig};

        uint32_t const size_before2 = fanout_mig.num_gates();
        mig_resubstitution2( fanout_mig, ps, &st );
        mig_resub = cleanup_dangling( mig_resub );

        d_time += to_seconds( st.time_total );

        if ( mig_resub.num_gates() < size_before2 )
        {
          mig_opt = mig_resub;
        }
      }
      if ( mig_opt.size() >= size_global_before )
        break;

      // mig_opt = det_randomize( mig_opt, i );
    }

    {
      resubstitution_params ps;
      resubstitution_stats st;
      ps.max_pis = 8u;
      ps.window_size = 12u;
      ps.max_inserts = 3u;
      // ps.progress = false;

      auto mig_resub = cleanup_dangling( mig_opt );

      depth_view depth_mig{mig_resub};
      fanout_view fanout_mig{depth_mig};

      uint32_t const size_before2 = fanout_mig.num_gates();
      mig_resubstitution2( fanout_mig, ps, &st );
      mig_resub = cleanup_dangling( mig_resub );

      d_time += to_seconds( st.time_total );

      if ( mig_resub.num_gates() < size_before2 )
      {
        mig_opt = mig_resub;
      }
    }


    // const auto cec = benchmark == "hyp" ? true : abc_cec( mig_opt, benchmark );
    const auto cec = true;
    // std::cout << cec << "\n";
    // if ( benchmark != "hyp" )
    // {
    //   std::cout << abc_cec( mig_a, benchmark ) << std::endl;
    //   std::cout << abc_cec( mig_ae, benchmark ) << std::endl;
    // }

    const uint32_t depth_d = depth_view( mig_opt ).depth();

    exp( benchmark, size_before, mig_opt.num_gates(), depth_before, depth_d, d_time );
  }

  exp.save();
  exp.table();

  return 0;
}
