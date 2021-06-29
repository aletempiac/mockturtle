#include "crypto_experiments.hpp"
#include "experiments.hpp"

#include <fmt/format.h>

#include <string>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
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
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
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

std::vector<std::string> local_benchmarks = {
    "adder",
    "bar",
    "div",
    "hyp",
    "log2",
    "max",
    "multiplier",
    "sin",
    "sqrt",
    "square" };
using namespace mockturtle;
using namespace experiments;

template<class Ntk>
bool abc_cec_benchmark( Ntk const& ntk, std::string const& benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/xmg__epfl_test.bench" );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/xmg__epfl_test.bench\"", benchmark );

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
  mockturtle::write_verilog( ntk, "/tmp/xmg__epfl_network.v" );
  system( fmt::format( "abc -q \"read /tmp/xmg__epfl_network.v; if -K {}; write_blif /tmp/xmg__epfl_output.blif\"", k ).c_str() );
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/xmg__epfl_output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
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

  exact_library_params eps;
  map_params ps1;
  ps1.skip_delay_round = true;
  ps1.required_time = std::numeric_limits<float>::max();
  map_stats st1;

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
          exact_library<aig_network, xag_npn_resynthesis<mockturtle::aig_network, mockturtle::aig_network, mockturtle::xag_npn_db_kind::aig_complete> > exact_aig_lib( aig_npn_resyn, eps );
          des = map( des, exact_aig_lib, ps1, &st1 );
          //mockturtle::cut_rewriting( des, aig_npn_resyn, cr_ps, &cr_st );
          des = mockturtle::cleanup_dangling( des);

          aig_resubstitution( des, ps, &st );
          des = cleanup_dangling( des );
      }
      if constexpr( std::is_same<typename Ntk::base_type, mockturtle::xag_network>::value )
      {
          std::cout << "xag" << std::endl;
          mockturtle::xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete> xag_npn_resyn;
          exact_library<xag_network, xag_npn_resynthesis<mockturtle::xag_network, mockturtle::xag_network, mockturtle::xag_npn_db_kind::xag_complete>> exact_xag_lib( xag_npn_resyn, eps );
          des = map( des, exact_xag_lib, ps1, &st1 );
          //mockturtle::cut_rewriting( des, xag_npn_resyn, cr_ps, &cr_st );
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
          exact_library<mig_network, mig_npn_resynthesis> exact_mig_lib( mig_npn_resyn, eps );
          des = map( des, exact_mig_lib, ps1, &st1 );
          //mockturtle::cut_rewriting( des, mig_npn_resyn, cr_ps, &cr_st );
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
          exact_library<xmg_network, xmg_npn_resynthesis> exact_xmg_lib( xmg_npn_resyn, eps );
          des = map( des, exact_xmg_lib, ps1, &st1 );
          //mockturtle::cut_rewriting( des, xmg_npn_resyn, cr_ps, &cr_st );
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

//void tech_map( std::string aig_or_klut, const uint32_t& cut_size, bool delay_round, bool req_time)
void tech_map()
{
    experiments::experiment<std::string, std::string, std::string, bool, bool, bool, bool>
         exp2( "RFET_area", "benchmark", "sd_rat", "sd_rat'", "cec1", "cec2", "cec3", "cec4" );

    experiments::experiment<std::string, float, float, float, float, float, float, float, float > exp( "Mapper Comparison", "benchmark", "Area AIG", "Area MIG", "Area XMG ", "Area XAG", "delay AIG", "delay MIG", "delay XMG", "delay XAG" );

  std::vector<mockturtle::gate> gates;
  if ( lorina::read_genlib( "smaller.genlib", mockturtle::genlib_reader( gates ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR IN" << std::endl;
    std::abort();
    return;
  }

  std::vector<mockturtle::map_superGate> supergates;
  mockturtle::super_info vals;
  if ( lorina::read_super( "orig_smaller.super", mockturtle::super_reader( supergates, vals ) ) != lorina::return_code::success )
  {
      std::cout << "ERROR IN super " << std::endl;
      std::abort();
      return;
  }

  mockturtle::tech_library_params lib_ps;
  lib_ps.very_verbose = false;
  mockturtle::tech_library<5> lib1( gates, lib_ps, supergates, vals);

  /* Option 1 */
  mockturtle::xmg_cost_params ps1, ps2;
  

  /* EPFL benchmarks */
  for ( const auto& benchmark : experiments::all_benchmarks() )
  {
      //if( benchmark != "leon2")
      //continue;

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
    ps1.reset();
    num_gate_profile( xmg, ps1 );
    ps1.report();
    auto size_before = xmg.num_gates();
    double sd_rat = ( double( ps1.actual_maj + ps1.actual_xor3 )/  size_before ) * 100;
    std::string sd_before = fmt::format( "{:>12.2f}", sd_rat );

    aig = cleanup_dangling( aig );
    mig = cleanup_dangling( mig );
    xmg = cleanup_dangling( xmg );
    xag = cleanup_dangling( xag );

    aig = ntk_optimization<mockturtle::aig_network> ( aig );
    mig = ntk_optimization<mockturtle::mig_network> ( mig );
    xmg = ntk_optimization<mockturtle::xmg_network> ( xmg );
		xag = ntk_optimization<mockturtle::xag_network> ( xag );

    aig = cleanup_dangling( aig );
    mig = cleanup_dangling( mig );
    xmg = cleanup_dangling( xmg );
		xag = cleanup_dangling( xag );

    ps2.reset();
    num_gate_profile( xmg, ps2);
    ps2.report();
    auto size_after = xmg.num_gates();
    sd_rat = ( double( ps2.actual_maj + ps2.actual_xor3 )/  size_after ) * 100;
    std::string sd_after = fmt::format( "{:>12.2f}", sd_rat );


    mockturtle::depth_view xmg_d{ xmg };
    mockturtle::depth_view mig_d{ mig };
    mockturtle::depth_view aig_d{ aig };
		mockturtle::depth_view xag_d{ xag };
    //printf( "###################################################\n" );
    //printf( "[i] AIG: n = %d   depth = %d\n",
    //        mig.size(), mig_d.depth() );
    //printf( "[i] MIG: n = %d   depth = %d\n",
    //        mig.size(), mig_d.depth() );
    //printf( "[i] XMG: n = %d   depth = %d\n",
    //        xmg.size(), xmg_d.depth() );
    fflush( stdout );

    mockturtle::map_params ps;
    //ps.verbose = true;
    ps.skip_delay_round = true;
    ps.required_time = std::numeric_limits<float>::max();
   
    mockturtle::map_stats aig_mst, mig_mst, xmg_mst, xag_mst;

    auto res1 = mockturtle::map( aig, lib1, ps, &aig_mst );
    auto res2 = mockturtle::map( mig, lib1, ps, &mig_mst );
    auto res3 = mockturtle::map( xmg, lib1, ps, &xmg_mst );
    auto res4 = mockturtle::map( xag, lib1, ps, &xag_mst );

    const auto cec1 =  abc_cec( res1, benchmark );
    const auto cec2 =  abc_cec( res2, benchmark );
    const auto cec3 =  abc_cec( res3, benchmark );
    const auto cec4 =  abc_cec( res4, benchmark );

    exp( benchmark,
            aig_mst.area, mig_mst.area, xmg_mst.area, xag_mst.area,
				aig_mst.delay, mig_mst.delay, xmg_mst.delay, xag_mst.delay );

    exp2 ( benchmark, sd_before, sd_after, cec1, cec2, cec3, cec4 );
    //exp2 (benchmark, sd_before, sd_after);
    //mockturtle::tech_mapping( xmg, lib2, ps, &mst );
    exp.save();
    exp.table();
    exp2.save();
    exp2.table();
  }
  exp.save();
  exp.table();
  exp2.save();
  exp2.table();
}

//int main( int argc, char* argv[])
int main() 
{

  tech_map( );//argv[1], std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
  return 0;
}
