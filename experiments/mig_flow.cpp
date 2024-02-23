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

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/circuit_validator.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/det_randomization.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/resyn_engines/mig_enumerative.hpp>
#include <mockturtle/algorithms/resyn_engines/mig_resyn.hpp>
#include <mockturtle/algorithms/rewrite.hpp>
#include <mockturtle/algorithms/sim_resub.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
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
      "mig_flow", "benchmark", "size", "size_mig", "depth", "depth_mig", "time_mig" );

  fmt::print( "[i] processing technology library\n" );

  /* library to map to MIGs */
  mig_npn_resynthesis resyn{ true };
  exact_library_params eps;
  eps.np_classification = false;
  eps.compute_dc_classes = true;
  exact_library<mig_network> exact_lib( resyn, eps );

  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    mig_network mig;
    if ( lorina::read_aiger( "optimized/" + benchmark + ".aig", aiger_reader( mig ) ) != lorina::return_code::success )
    {
      continue;
    }

    const uint32_t size_before = mig.num_gates();
    const uint32_t depth_before = depth_view( mig ).depth();

    float d_time = 0;
    mig_network mig_opt = cleanup_dangling( mig );

    uint32_t i = 3;
    while ( i-- > 0 )
    {
      map_params ps;
      map_stats st;
      ps.skip_delay_round = true;
      ps.required_time = std::numeric_limits<float>::max();
      ps.ela_rounds = 2;
      ps.enable_logic_sharing = true;
      ps.use_dont_cares = true;
      ps.window_size = 12;
      ps.logic_sharing_cut_limit = 1;
      ps.cut_enumeration_ps.cut_limit = 8;

      auto size_before = mig_opt.size();
      auto mig_map = map( mig_opt, exact_lib, ps, &st );

      d_time += to_seconds( st.time_total );

      if ( mig_map.size() >= size_before )
        break;

      mig_opt = mig_map;
    }

    i = 3;
    while ( i-- > 0 )
    {
      auto size_before = mig_opt.size();

      rewrite_params rps;
      rps.use_dont_cares = i == 1;
      rps.allow_zero_gain = true;
      rps.window_size = 8;
      rewrite_stats rst;
      rewrite( mig_opt, exact_lib, rps, &rst );

      d_time += to_seconds( rst.time_total );

      if ( mig_opt.size() >= size_before )
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
        // ps.window_size = 8u;
        ps.max_inserts = 2u;

        auto mig_resub = cleanup_dangling( mig_opt );

        depth_view depth_mig{ mig_resub };
        fanout_view fanout_mig{ depth_mig };

        uint32_t const size_before2 = fanout_mig.num_gates();
        mig_resubstitution( fanout_mig, ps, &st );
        mig_resub = cleanup_dangling( mig_resub );

        d_time += to_seconds( st.time_total );

        if ( mig_resub.num_gates() < size_before2 )
        {
          mig_opt = mig_resub;
        }
      }
      if ( mig_opt.size() >= size_global_before )
        break;
    }

    {
      resubstitution_params ps;
      resubstitution_stats st;
      ps.max_pis = 8u;
      ps.max_inserts = std::numeric_limits<uint32_t>::max();

      depth_view depth_mig{ mig_opt };
      fanout_view fanout_mig{ depth_mig };

      using resub_view_t = fanout_view<depth_view<mig_network>>;
      using resyn_engine_t = mig_resyn_topdown<kitty::partial_truth_table, mig_resyn_static_params>;

      using validator_t = circuit_validator<resub_view_t, bill::solvers::bsat2, false, true, false>;
      using resub_impl_t = typename detail::resubstitution_impl<resub_view_t, typename detail::simulation_based_resub_engine<resub_view_t, validator_t, resyn_engine_t>>;

      typename resub_impl_t::engine_st_t engine_st;
      typename resub_impl_t::collector_st_t collector_st;

      resub_impl_t p( fanout_mig, ps, st, engine_st, collector_st );
      p.run();
      d_time += to_seconds( st.time_total );
      mig_opt = cleanup_dangling( mig_opt );
    }

    const auto cec = benchmark == "hyp" ? true : abc_cec( mig_opt, benchmark );
    const uint32_t depth_d = depth_view( mig_opt ).depth();
    std::cout << fmt::format( "Size = {:8d}\t Depth = {:8d}\t cec = {}\n", mig_opt.num_gates(), depth_d, cec );

    exp( benchmark, size_before, mig_opt.num_gates(), depth_before, depth_d, d_time );
  }

  exp.save();
  exp.table();

  return 0;
}