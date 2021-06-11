#include "crypto_experiments.hpp"
#include "sd_experiments.hpp"
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
#include <mockturtle/io/write_bench.hpp>
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
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/algorithms/aig_resub.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/xmg_resub.hpp>
#include <mockturtle/algorithms/xag_resub_withDC.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>

std::vector<std::string> local_benchmarks = {
	"benchmarks_4_4_3_1",
	"benchmarks_4_4_3_2",
	"benchmarks_4_4_3_3",
	"benchmarks_4_4_3_4",
	"benchmarks_4_4_3_5",
	"benchmarks_4_4_3_6",
	"benchmarks_4_4_3_7",
	"benchmarks_4_4_3_8",
	"benchmarks_4_4_3_9", 
	"benchmarks_4_4_3_10",
	//};
	"benchmarks_12_512_131_10",
	"benchmarks_12_512_131_1",
	"benchmarks_12_512_131_2",
	"benchmarks_12_512_131_3",
	"benchmarks_12_512_131_4",
	"benchmarks_12_512_131_5",
	"benchmarks_12_512_131_6",
	"benchmarks_12_512_131_7",
	"benchmarks_12_512_131_8",
	"benchmarks_12_512_131_9",
	"benchmarks_128_231_131_10",
	"benchmarks_128_231_131_1",
	"benchmarks_128_231_131_2",
	"benchmarks_128_231_131_3",
	"benchmarks_128_231_131_4",
	"benchmarks_128_231_131_5",
	"benchmarks_128_231_131_6",
	"benchmarks_128_231_131_7",
	"benchmarks_128_231_131_8",
	"benchmarks_128_231_131_9",
	"benchmarks_255_399_131_10",
	"benchmarks_255_399_131_1",
	"benchmarks_255_399_131_2",
	"benchmarks_255_399_131_3",
	"benchmarks_255_399_131_4",
	"benchmarks_255_399_131_5",
	"benchmarks_255_399_131_6",
	"benchmarks_255_399_131_7",
	"benchmarks_255_399_131_8",
	"benchmarks_255_399_131_9" };

using namespace mockturtle;
using namespace experiments;

	template<class Ntk>
bool abc_cec_benchmark( Ntk const& ntk, std::string const& benchmark )
{
	mockturtle::write_bench( ntk, "/tmp/test.bench" );
	std::string command = fmt::format( "abc -q \"cec -n {} /tmp/test.bench\"", benchmark );

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
	mockturtle::write_verilog( ntk, "/tmp/network.v" );
	system( fmt::format( "abc -q \"/tmp/network.v; &get; &if -a -K {}; &put; write_blif /tmp/output.blif\"", k ).c_str() );
	mockturtle::klut_network klut;
	if ( lorina::read_blif( "/tmp/output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
	{
		std::cout << "ERROR LUT" << std::endl;
		std::abort();
		return klut;
	}
	return klut;
}

	template<typename Ntk>
Ntk ntk_optimization( Ntk const& ntk )
{
	auto des = ntk;
	resubstitution_params ps;
	resubstitution_stats st;

	ps.max_pis = 8u;
	ps.max_inserts = 1u;
	ps.progress = false;
	mockturtle::cut_rewriting_params cr_ps;
	mockturtle::cut_rewriting_stats cr_st;
	cr_ps.cut_enumeration_ps.cut_size = 4;

	float improv = 0;
	float improv_per = 0;
	uint32_t iter = 0;

	while(true)
	{
		auto const size_before = des.size();
		if constexpr( std::is_same<typename Ntk::base_type, mockturtle::aig_network>::value )
		{
			std::cout << "aig" << std::endl;
			mockturtle::xag_npn_resynthesis<mockturtle::aig_network, mockturtle::aig_network, mockturtle::xag_npn_db_kind::aig_complete> aig_npn_resyn;
			mockturtle::cut_rewriting( des, aig_npn_resyn, cr_ps, &cr_st );
			des = mockturtle::cleanup_dangling( des);

			aig_resubstitution( des, ps, &st );
			des = cleanup_dangling( des );
		}
		if constexpr( std::is_same<typename Ntk::base_type, mockturtle::xag_network>::value )
		{
			std::cout << "xag" << std::endl;
			mockturtle::xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete> xag_npn_resyn;
			mockturtle::cut_rewriting( des, xag_npn_resyn, cr_ps, &cr_st );
			des = mockturtle::cleanup_dangling( des);

			using view_t = depth_view<fanout_view<xag_network>>;
			fanout_view<xag_network> fanout_view{des};
			view_t resub_view{fanout_view};
			resubstitution_minmc_withDC( resub_view , ps, &st);

			des = cleanup_dangling( des );
		}
		if constexpr( std::is_same<typename Ntk::base_type, mockturtle::mig_network>::value )
		{
			std::cout << "mig" << std::endl;
			mockturtle::mig_npn_resynthesis mig_npn_resyn{ true };
			mockturtle::cut_rewriting( des, mig_npn_resyn, cr_ps, &cr_st );
			des = mockturtle::cleanup_dangling( des);
			depth_view depth_mig{des};
			fanout_view fanout_mig{depth_mig};

			mig_resubstitution( fanout_mig, ps, &st );
			des = cleanup_dangling( des );
		}
		if constexpr( std::is_same<typename Ntk::base_type, mockturtle::xmg_network>::value )
		{
			std::cout << "xmg" << std::endl;
			mockturtle::xmg_npn_resynthesis xmg_npn_resyn;
			mockturtle::cut_rewriting( des, xmg_npn_resyn, cr_ps, &cr_st );
			des = mockturtle::cleanup_dangling( des);

			xmg_resubstitution( des, ps, &st );
			des = mockturtle::cleanup_dangling( des );
		}

		std::cout << "size after and before  "<< des.size() << " " << size_before << std::endl ;
		improv =  size_before - des.size(); 
		auto diff = std::abs(improv);
		improv_per = 100 * (double(diff/size_before)); //100 * (double((std::abs(improv)) improv/ size_before ));
		std::cout << " improvement " << improv << " improv_per " << improv_per << std::endl;
		std::cout << "Iterations # " << iter++ << std::endl; 
		if (improv_per <= 0.5 )
			break;
	}

	return des;
}

void tech_map()
{

	experiments::experiment<std::string, float, float, float, float, float, float, float, float > exp( "Mapper Comparison", "benchmark", "Area AIG", "Area MIG", "Area XMG ", "Area XAG", "delay AIG", "delay MIG", "delay XMG", "delay XAG" );

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

	for ( const auto& b : local_benchmarks )
	{
		if (b != "benchmarks_4_4_3_8")
		    continue;
		std::string filename{"../experiments/self_dual_benchmarks/"};
		filename = filename + b + ".v";


		fmt::print( "[i] processing {}\n", filename);
		fflush( stdout );

		mockturtle::xmg_network xmg;
		mockturtle::aig_network aig;
		mockturtle::mig_network mig;
		mockturtle::xag_network xag;
    mockturtle::xmg_cost_params ps1, ps2;

		/* Option 2 */
		mockturtle::xag_npn_resynthesis<mockturtle::aig_network, mockturtle::aig_network, mockturtle::xag_npn_db_kind::aig_complete> aig_npn_resyn;
		mockturtle::xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete> xag_npn_resyn;
		mockturtle::xmg_npn_resynthesis npn_resyn;
		mockturtle::mig_npn_resynthesis mig_npn_resyn{ true };

		if ( lorina::read_verilog( filename, mockturtle::verilog_reader( xmg) ) != lorina::return_code::success )
		{
			std::cout << "ERROR IN reading benchmark" << std::endl;
			std::abort();
			return;
		}

		mockturtle::write_verilog( xmg, "resyn_fail.v");

		auto klut = lut_map( xmg, 4u );
		std::cout << "Before Resyn done " << std::endl;

		mockturtle::write_bench( klut, "resyn_bench.v");
		//imig = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, npn_resyn );
		aig = mockturtle::node_resynthesis<mockturtle::aig_network>( klut, aig_npn_resyn );
		aig = cleanup_dangling( aig);

		mig = mockturtle::node_resynthesis<mockturtle::mig_network>( klut, mig_npn_resyn );
		mig = cleanup_dangling( mig );

		xag = mockturtle::node_resynthesis<mockturtle::xag_network>( klut, xag_npn_resyn );
		xag = cleanup_dangling( xag );

		std::cout << "Resyn done " << std::endl;

		aig = ntk_optimization<mockturtle::aig_network> ( aig );
		mig = ntk_optimization<mockturtle::mig_network> ( mig );
		xmg = ntk_optimization<mockturtle::xmg_network> ( xmg );
		xag = ntk_optimization<mockturtle::xag_network> ( xag );

		aig = cleanup_dangling( aig );
		mig = cleanup_dangling( mig );
		xmg = cleanup_dangling( xmg );
		xag = cleanup_dangling( xag );

		mockturtle::depth_view xmg_d{ xmg };
		mockturtle::depth_view mig_d{ mig };
		mockturtle::depth_view aig_d{ aig };
		mockturtle::depth_view xag_d{ xag };
		printf( "###################################################\n" );
		printf( "[i] AIG: n = %d   depth = %d\n",
				mig.size(), mig_d.depth() );
		printf( "[i] MIG: n = %d   depth = %d\n",
				mig.size(), mig_d.depth() );
		printf( "[i] XMG: n = %d   depth = %d\n",
				xmg.size(), xmg_d.depth() );
		printf( "[i] XAG: n = %d   depth = %d\n",
				xag.size(), xag_d.depth() );
		fflush( stdout );

		mockturtle::map_params ps;
		ps.cut_enumeration_ps.cut_size = 6;
		ps.cut_enumeration_ps.cut_limit = 25;
		ps.verbose = true;
		ps.skip_delay_round = true;
		mockturtle::map_stats aig_mst, mig_mst, xmg_mst, xag_mst;

		mockturtle::tech_mapping( aig, lib1, ps, &aig_mst );
		fflush( stdout );
		mockturtle::tech_mapping( mig, lib1, ps, &mig_mst );
		fflush( stdout );
		mockturtle::tech_mapping( xmg, lib1, ps, &xmg_mst );
		fflush( stdout );
		mockturtle::tech_mapping( xag, lib1, ps, &xag_mst );
		fflush( stdout );

		exp( b,
				aig_mst.area, mig_mst.area, xmg_mst.area, xag_mst.area,
				aig_mst.delay, mig_mst.delay, xmg_mst.delay, xag_mst.delay );

		////mockturtle::tech_mapping( xmg, lib2, ps, &mst );
		exp.save();
		exp.table();
	}
	exp.save();
	exp.table();
}

int main()
{
	tech_map();
	return 0;
}