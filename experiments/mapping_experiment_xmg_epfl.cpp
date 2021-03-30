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
//
//
//std::vector<std::string> local_benchmarks_iwls = {
//  "aes_core",
//  "mem_ctrl",
//  "voter"
//};
//
//
//std::vector<std::string> benchmarks_aqfp_v = {
//    //"5xp1",
//    "C1908_orig",
//    "C432_orig",
//    "C880_orig",
//    "C5315_orig",
//    "count_orig",
//    //"dist_orig",
//    "i5_orig",
//    "i6_orig",
//    "k2_orig",
//    "majority_orig",
//    "x1_orig"
//};
//
using namespace mockturtle;
using namespace experiments;

template<class Ntk>
bool abc_cec_benchmark( Ntk const& ntk, std::string const& benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/xmg_test.bench" );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/xmg_test.bench\"", benchmark );

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

template<typename Ntk>
Ntk ntk_optimization( Ntk const& ntk )
{
  auto des = ntk;
  resubstitution_params ps;
  resubstitution_stats st;

  ps.max_pis = 8u;
  ps.max_inserts = 1u;
  ps.progress = true;
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
          //std::cout << "aig" << std::endl;
          exact_resynthesis_params eps;
          //eps.cache = std::make_shared<exact_resynthesis_params::cache_map_t>();
          exact_aig_resynthesis<aig_network> aig_exact( false, eps );
          mockturtle::cached_resynthesis<mockturtle::aig_network, decltype( aig_exact )> cached_aig_exact( aig_exact, 4, "exact_aig_cache4_cr.v" );
          mockturtle::cut_rewriting( des, cached_aig_exact, cr_ps, &cr_st );
          des = mockturtle::cleanup_dangling( des);

          aig_resubstitution( des, ps, &st );
          des = cleanup_dangling( des );
      }
      if constexpr( std::is_same<typename Ntk::base_type, mockturtle::mig_network>::value )
      {
          //std::cout << "mig" << std::endl;
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
          //std::cout << "xmg" << std::endl;
          mockturtle::xmg3_npn_resynthesis<mockturtle::xmg_network> xmg_npn_resyn;
          mockturtle::cut_rewriting( des, xmg_npn_resyn, cr_ps, &cr_st );
          des = mockturtle::cleanup_dangling( des);
          
          xmg_resubstitution( des, ps, &st );
          des = mockturtle::cleanup_dangling( des );
      }

      //std::cout << "size after and before  "<< des.size() << " " << size_before << std::endl ;
      improv =  size_before - des.size(); 
      auto diff = std::abs(improv);
      improv_per = 100 * (double(diff/size_before)); //100 * (double((std::abs(improv)) improv/ size_before ));
      //std::cout << " improvement " << improv << " improv_per " << improv_per << std::endl;
      //std::cout << "Iterations # " << iter++ << std::endl; 
      if (improv_per <= 0.5 )
          break;
  }

  return des;
}

void tech_map( std::string aig_or_klut, const uint32_t& cut_size, bool delay_round, bool req_time)
{

    std::string filename = "epfl";
    filename = filename + aig_or_klut + std::to_string(cut_size) + (delay_round == 0 ? "_false" : "_true") + (req_time == 0 ? "_def": "_max") + ".txt" ;
    //std::cout << "log file" << filename;
    std::ofstream outs;
    outs.open(filename.c_str());

    outs << "aig(0) or klut(1)   "      << aig_or_klut << std::endl;
    outs << "cut size = "               << cut_size    << std::endl;
    outs << "delay round (0/1)=  "      << (delay_round ? "true" : "false") << std::endl;
    outs << "required time (def/max)= " << (req_time ? "true" : "false")  << std::endl;

    experiments::experiment<std::string, std::string, std::string>
         exp2( "RFET_area", "benchmark", "sd_rat", "sd_rat'");

  experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, float, float, float, float, float, float, float, float, float> exp( "Mapper Comparison", "benchmark", "size AIG", "size MIG", "Size XMG", "depth AIG", "depth MIG", "depth XMG", "Area AIG", "Area MIG", "Area XMG ", "delay AIG", "delay MIG", "delay XMG" );

  std::vector<mockturtle::gate> gates1, gates2;
  if ( lorina::read_genlib( "smaller.genlib", mockturtle::genlib_reader( gates1 ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR IN" << std::endl;
    std::abort();
    return;
  }
  //if ( lorina::read_genlib( "mcnc_smaller.genlib", mockturtle::genlib_reader( gates2 ) ) != lorina::return_code::success )
  //{
  //    std::cout << "ERROR IN" << std::endl;
  //    std::abort();
  //    return;
  //}

  //for ( auto const& g : gates1 )
  //{
  //  std::cout << g.name << std::endl;
  //}
  mockturtle::tech_library_params lib_ps;
  lib_ps.very_verbose = false;
  mockturtle::tech_library<5> lib1( gates1, lib_ps );
  //mockturtle::tech_library<5> lib2( gates2, lib_ps );

  /* Option 1 */
  mockturtle::exact_xmg_resynthesis_params xmg3_exact_ps;
  xmg3_exact_ps.use_xor3 = true;
  xmg3_exact_ps.num_candidates = 10u;
  mockturtle::exact_xmg_resynthesis<mockturtle::xmg_network> xmg3_exact( xmg3_exact_ps );
  mockturtle::cached_resynthesis<mockturtle::xmg_network, decltype( xmg3_exact )> cached_xmg3_exact( xmg3_exact, 4, "exact_xmg3_cache4.v" );
  mockturtle::xmg_cost_params ps1, ps2;
  

  exact_resynthesis_params eps;
  //eps.cache = std::make_shared<exact_resynthesis_params::cache_map_t>();
  exact_aig_resynthesis<aig_network> aig_exact( false, eps );
  mockturtle::cached_resynthesis<mockturtle::aig_network, decltype( aig_exact )> cached_aig_exact( aig_exact, 4, "exact_aig_cache4_cr.v" );

  /* EPFL benchmarks */
  for ( const auto& benchmark : experiments::epfl_benchmarks() )
  {
      //if( benchmark != "adder")
      //continue;

    /* Crypto Benchmarks */
    //for ( auto const& benchmark : crypto_experiments::crypto_benchmarks( ))
    //{
    //  //if (benchmark != "sd_test")
    //  //    continue;
    fmt::print( "[i] processing {}\n", benchmark );
    fflush( stdout );

    mockturtle::xmg_network xmg;
    mockturtle::aig_network aig;
    mockturtle::mig_network mig;

    /* Option 2 */
    mockturtle::xmg3_npn_resynthesis<mockturtle::xmg_network> npn_resyn;
    mockturtle::mig_npn_resynthesis mig_npn_resyn{ true };

    //if ( lorina::read_verilog( crypto_experiments::benchmark_path( benchmark ), mockturtle::verilog_reader( aig ) ) != lorina::return_code::success )
    //{
    //    std::cout << "ERROR IN reading benchmark" << std::endl;
    //    std::abort();
    //    return;
    //}
    

    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), mockturtle::aiger_reader( aig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR IN reading benchmark" << std::endl;
      std::abort();
      return;
    }
    balancing_params sps;
    balancing_stats st4;
    sop_rebalancing<aig_network> sop_balancing;    
    aig = balancing( aig, {sop_balancing}, sps, &st4 );

    auto klut = lut_map( aig, 4u );


    /* Calling Resynthesis engine */
    if(aig_or_klut == "aig")
    {
        xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( aig, cached_xmg3_exact );
        xmg = cleanup_dangling( xmg );

        aig = mockturtle::node_resynthesis<mockturtle::aig_network>( aig, cached_aig_exact );
        aig = cleanup_dangling( aig );

        mig = mockturtle::node_resynthesis<mockturtle::mig_network>( aig, mig_npn_resyn );
        mig = cleanup_dangling( mig );
    }
    else
    {
        xmg = mockturtle::node_resynthesis<mockturtle::xmg_network>( klut, cached_xmg3_exact );
        xmg = cleanup_dangling( xmg );

        aig = mockturtle::node_resynthesis<mockturtle::aig_network>( klut, cached_aig_exact );
        aig = cleanup_dangling( aig );

        mig = mockturtle::node_resynthesis<mockturtle::mig_network>( klut, mig_npn_resyn );
        mig = cleanup_dangling( mig );

    }

    /* measuring sd ratios */
    ps1.reset();
    num_gate_profile( xmg, ps1 );
    ps1.report();
    auto size_before = xmg.num_gates();
    double sd_rat = ( double( ps1.actual_maj + ps1.actual_xor3 )/  size_before ) * 100;
    std::string sd_before = fmt::format( "{}/{} = {}", ( ps1.actual_maj + ps1.actual_xor3 ),  size_before, sd_rat );


    aig = ntk_optimization<mockturtle::aig_network> ( aig );
    mig = ntk_optimization<mockturtle::mig_network> ( mig );
    xmg = ntk_optimization<mockturtle::xmg_network> ( xmg );

    aig = cleanup_dangling( aig );
    mig = cleanup_dangling( mig );
    xmg = cleanup_dangling( xmg );

    ps2.reset();
    num_gate_profile( xmg, ps2);
    ps2.report();
    auto size_after = xmg.num_gates();
    sd_rat = ( double( ps2.actual_maj + ps2.actual_xor3 )/  size_after ) * 100;
    std::string sd_after = fmt::format( "{}/{} = {}", ( ps2.actual_maj + ps2.actual_xor3 ),  size_after, sd_rat );


    mockturtle::depth_view xmg_d{ xmg };
    mockturtle::depth_view mig_d{ mig };
    mockturtle::depth_view aig_d{ aig };
    //printf( "###################################################\n" );
    //printf( "[i] AIG: n = %d   depth = %d\n",
    //        mig.size(), mig_d.depth() );
    //printf( "[i] MIG: n = %d   depth = %d\n",
    //        mig.size(), mig_d.depth() );
    //printf( "[i] XMG: n = %d   depth = %d\n",
    //        xmg.size(), xmg_d.depth() );
    fflush( stdout );

    mockturtle::map_params ps;
    ps.cut_enumeration_ps.cut_size = cut_size;
    ps.cut_enumeration_ps.cut_limit = 25;
    ps.verbose = true;
    if (delay_round)
        ps.skip_delay_round = true;
    else
        ps.skip_delay_round = false;
    if (req_time)
        ps.required_time = std::numeric_limits<float>::max();
   
    mockturtle::map_stats aig_mst, mig_mst, xmg_mst;

    mockturtle::tech_mapping( aig, lib1, ps, &aig_mst );
    fflush( stdout );
    mockturtle::tech_mapping( mig, lib1, ps, &mig_mst );
    fflush( stdout );
    mockturtle::tech_mapping( xmg, lib1, ps, &xmg_mst );
    fflush( stdout );

    exp( benchmark, aig.size(), mig.size(), xmg.size(),
         aig_d.depth(), mig_d.depth(), xmg_d.depth(),
         aig_mst.area, mig_mst.area, xmg_mst.area,
         aig_mst.delay, mig_mst.delay, xmg_mst.delay );

    exp2 (benchmark, sd_before, sd_after);
    //mockturtle::tech_mapping( xmg, lib2, ps, &mst );
    exp.save();
    exp.table();
    exp2.save();
    exp2.table();
  }
  outs.close();
  outs.open(filename.c_str(), std::ios::app);
  exp.save("1");
  exp.table("1", outs); 
  exp2.save("1");
  exp2.table("1",outs); 
  outs.close();
}

int main( int argc, char* argv[])
{
  std::cout << "aig(0) or klut(1)   "      << argv[1] << std::endl;
  std::cout << "cut size = "               << argv[2] << std::endl;
  std::cout << "delay round (0/1)=  "      << argv[3] << std::endl;
  std::cout << "required time (def/max)= " << argv[4] << std::endl;

  tech_map( argv[1], std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
  return 0;
}
