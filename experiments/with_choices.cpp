#include "crypto_experiments.hpp"
#include "experiments.hpp"

#include <fmt/format.h>

#include <string>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/lorina.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
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
#include <mockturtle/algorithms/xag_resub_withDC.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/balancing/sop_balancing.hpp>
#include <mockturtle/properties/xmgcost.hpp>


using namespace mockturtle;
using namespace experiments;

    template<class Ntk>
bool abc_cec_benchmark( Ntk const& ntk, std::string const& benchmark )
{
    mockturtle::write_bench( ntk, "/tmp/xmg_all_test.bench" );
    std::string command = fmt::format( "abc -q \"cec -n {} /tmp/xmg_all_test.bench\"", benchmark );

    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
    if ( !pipe )
    {
        throw std::runtime_error( "popen() failed" );
    }
    while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
    {
        result += buffer.data();
    }
    std::cout << result << std::endl;

    return result.size() >= 23 && result.substr( 0u, 23u ) == "Networks are equivalent";
}

    template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
    mockturtle::write_verilog( ntk, "/tmp/xmg_all_network.v" );
    system( fmt::format( "abc -q \"read /tmp/xmg_all_network.v; if -K {}; write_blif /tmp/xmg_all_output.blif\"", k ).c_str() );
    mockturtle::klut_network klut;
    if ( lorina::read_blif( "/tmp/xmg_all_output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
    {
        std::cout << "ERROR LUT" << std::endl;
        std::abort();
        return klut;
    }
    return klut;
}

void tech_map()
{
    std::vector<mockturtle::gate> gates1;
    if ( lorina::read_genlib( "smaller.genlib", mockturtle::genlib_reader( gates1 ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR IN" << std::endl;
    std::abort();
    return;
  }

    mockturtle::tech_library_params lib_ps;
    lib_ps.very_verbose = false;
    lib_ps.compute_supergates = true;
    mockturtle::tech_library<6> lib1( gates1, lib_ps );

    mockturtle::xmg_cost_params ps1, ps2;

    experiments::experiment<std::string, std::string, std::string>
        exp2( "RFET_area", "benchmark", "sd_rat", "sd_rat'");

    experiments::experiment<std::string, float, float, float, float, float, float, float, float > exp( "Mapper Comparison", "benchmark", "Area AIG", "Area MIG", "Area XMG ", "Area XAG", "delay AIG", "delay MIG", "delay XMG", "delay XAG" );


    for ( const auto& benchmark : experiments::epfl_benchmarks() )
    {
        //if( benchmark != "bar")
        //    continue;

        fmt::print( "[i] processing {}\n", benchmark );
        fflush( stdout );

        mockturtle::xmg_network xmg;
        mockturtle::aig_network aig;
        mockturtle::mig_network mig;
        mockturtle::xag_network xag;

        /* Option 2 */
        mockturtle::xag_npn_resynthesis<mockturtle::aig_network, mockturtle::aig_network, mockturtle::xag_npn_db_kind::aig_complete> aig_npn_resyn;
        mockturtle::xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete> xag_npn_resyn;
        mockturtle::xmg_npn_resynthesis npn_resyn;
        mockturtle::mig_npn_resynthesis mig_npn_resyn{ true };


        if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
        {
            std::cout << "ERROR IN reading benchmark" << std::endl;
            std::abort();
            return;
        }
        auto klut = lut_map( aig, 4u );

        /* Calling Resynthesis engine */
        xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn);
        xmg = cleanup_dangling( xmg );

        mig = mockturtle::node_resynthesis<mockturtle::mig_network>( klut, mig_npn_resyn );
        mig = cleanup_dangling( mig );

        xag = mockturtle::node_resynthesis<mockturtle::xag_network>( klut, xag_npn_resyn );
        xag = cleanup_dangling( xag );

        /* measuring sd ratios */
        //ps1.reset();
        //num_gate_profile( xmg, ps1 );
        //ps1.report();
        //auto size_before = xmg.num_gates();
        //double sd_rat = ( double( ps1.actual_maj + ps1.actual_xor3 )/  size_before ) * 100;
        //std::string sd_before = fmt::format( "{}/{} = {}", ( ps1.actual_maj + ps1.actual_xor3 ),  size_before, sd_rat );

        //ps2.reset();
        //num_gate_profile( xmg, ps2);
        //ps2.report();
        //auto size_after = xmg.num_gates();
        //sd_rat = ( double( ps2.actual_maj + ps2.actual_xor3 )/  size_after ) * 100;
        //std::string sd_after = fmt::format( "{}/{} = {}", ( ps2.actual_maj + ps2.actual_xor3 ),  size_after, sd_rat );


        mockturtle::depth_view xmg_d{ xmg };
        mockturtle::depth_view mig_d{ mig };
        mockturtle::depth_view aig_d{ aig };
        mockturtle::depth_view xag_d{ xag };

        /* Implementing choices */

        auto best_size = xmg.size();
        mockturtle::depth_view xmg_d_tmp{xmg};
        auto best_depth = xmg_d_tmp.depth();
        mockturtle::functional_reduction_params frp;
        mockturtle::functional_reduction_stats st;
        frp.compute_equivalence_classes = true;
        auto eqpairs = mockturtle::functional_reduction_choices( xmg, frp, &st );
        //mockturtle::functional_reduction( mig, frp, &st );

        mockturtle::choice_view cxmg{xmg};
        mockturtle::reduce_choice_network( cxmg, eqpairs );
        mockturtle::improve_representatives( cxmg );
        xmg = cleanup_dangling( xmg );
        mockturtle::choice_view<mockturtle::xmg_network> cxmg2 = mockturtle::levelize_choice_network( cxmg );

        mockturtle::map_params ps;
        ps.cut_enumeration_ps.cut_size = 6;
        ps.cut_enumeration_ps.cut_limit = 25;
        ps.verbose = false;
        ps.skip_delay_round = true;
        ps.required_time = std::numeric_limits<float>::max();

        mockturtle::map_stats aig_mst, mig_mst, xmg_mst, xag_mst;

        mockturtle::tech_mapping( aig, lib1, ps, &aig_mst );
        fflush( stdout );
        mockturtle::tech_mapping( mig, lib1, ps, &mig_mst );
        fflush( stdout );
        mockturtle::tech_mapping( xmg, lib1, ps, &xmg_mst );
        fflush( stdout );
        mockturtle::tech_mapping( xag, lib1, ps, &xag_mst );
        fflush( stdout );

        exp( benchmark,
                aig_mst.area, mig_mst.area, xmg_mst.area, xag_mst.area,
                aig_mst.delay, mig_mst.delay, xmg_mst.delay, xag_mst.delay );

        //exp2 (benchmark, sd_before, sd_after);
        exp.save();
        exp.table();
        //exp2.save();
        //exp2.table();
    }
    exp.save();
    exp.table(); 
    //exp2.save();
    //exp2.table(); 
}

int main() 
{
    tech_map();
    return 0;
}





