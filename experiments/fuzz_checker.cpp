#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/algorithms/network_fuzz_tester.hpp>
#include <mockturtle/generators/random_logic_generator.hpp>
#include "crypto_experiments.hpp"
#include "experiments.hpp"

#include <fmt/format.h>

#include <string>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/lorina.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/cached.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg3_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/algorithms/tech_mapper.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/choice_utils.hpp>
#include <mockturtle/utils/stopwatch.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/choice_view.hpp>
#include <mockturtle/views/depth_choice_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/topo_view.hpp>

#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/properties/xmgcost.hpp>


#include <lorina/lorina.hpp>
#include <fmt/format.h>
#include <optional>

using namespace mockturtle;

template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
  mockturtle::write_verilog( ntk, "/tmp/xmg_network.v" );
  system( fmt::format( "abc -q \"/tmp/xmg_network.v; &get; &if -a -K {}; &put; write_blif /tmp/xmg_output.blif\"", k ).c_str() );
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/xmg_output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR LUT" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}

int main()
{
    std::vector<mockturtle::gate> gates1, gates2;
    if ( lorina::read_genlib( "smaller.genlib", mockturtle::genlib_reader( gates1 ) ) != lorina::return_code::success )
    {
        std::cout << "ERROR IN" << std::endl;
        std::abort();
        return 0;
    }


    auto gen = mockturtle::default_random_xag_generator( );
    mockturtle::xmg_npn_resynthesis npn_resyn;
    mockturtle::xmg_network xmg;

    auto iteration = 0u;
    auto map_fn = [&] (mockturtle::xag_network xag) -> bool {
        mockturtle::tech_library_params lib_ps;
        lib_ps.very_verbose = false;
        mockturtle::tech_library<5> lib1( gates1, lib_ps );

        ++iteration;    
        std::cout << "Iterations = " << iteration << std::endl;
        auto klut = lut_map( xag, 4 );
        xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn);
        xmg = mockturtle::cleanup_dangling( xmg );
        balancing_params sps;
        balancing_stats st4;
        sop_rebalancing<xmg_network> xmg_balancing;    
        xmg = balancing( xmg, {xmg_balancing}, sps, &st4 );

        sop_rebalancing<xag_network> xag_balancing;    
        xag = balancing( xag, {xag_balancing}, sps, &st4 );

        mockturtle::map_params mps;
        mps.cut_enumeration_ps.cut_size = 4;
        mps.cut_enumeration_ps.cut_limit = 16;
        mps.verbose = true;
        mps.skip_delay_round = true;
        mockturtle::map_stats xag_mst, xmg_mst;

        std::cout << "tech mapping with xag" << std::endl << std::endl;
        fflush( stdout );
        mockturtle::tech_mapping( xag, lib1, mps, &xag_mst );
        std::cout << "tech mapping with XMG" << std::endl << std::endl;
        fflush( stdout );
        mockturtle::tech_mapping( xmg, lib1, mps, &xmg_mst );

        std::cout << "xag area \t " << xag_mst.area << std::endl
                  << "xmg area \t " << xmg_mst.area << std::endl;

        if (xag_mst.area < xmg_mst.area)
            return false;
        else 
            return true;
    };


    mockturtle::fuzz_tester_params ps;
    ps.num_iterations = 50;
    mockturtle::network_fuzz_tester fuzzer( gen, ps );
    fuzzer.run( map_fn );
    return 0;
}
