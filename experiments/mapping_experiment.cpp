#include "experiments.hpp"

#include <fmt/format.h>

#include <string>
#include <vector>

#include <lorina/lorina.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
// #include <mockturtle/algorithms/node_resynthesis/xmg_npn.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/views/topo_view.hpp>
#include <mockturtle/views/choice_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/depth_choice_view.hpp>
#include <mockturtle/utils/choice_utils.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/tech_mapper.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/utils/stopwatch.hpp>



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
  "square"
};


std::vector<std::string> local_benchmarks_iwls = {
  "aes_core",
  "mem_ctrl",
  "voter"
};


std::vector<std::string> benchmarks_aqfp_v = {
    //"5xp1",
    "C1908_orig",
    "C432_orig",
    "C880_orig",
    "C5315_orig",
    "count_orig",
    //"dist_orig",
    "i5_orig",
    "i6_orig",
    "k2_orig",
    "majority_orig",
    "x1_orig"
};


template<class Ntk>
bool abc_cec_benchmark( Ntk const& ntk, std::string const& benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "../../abc -q \"cec -n {} /tmp/test.bench\"", benchmark );

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
  system( fmt::format( "../../abc -q \"/tmp/network.v; &get; &if -a -K {}; &put; write_blif /tmp/output.blif\"", k ).c_str() );
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR LUT" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}


template<class Ntk>
void map_core( Ntk const& imig, mockturtle::exact_library<mockturtle::mig_network, mockturtle::mig_npn_resynthesis, 4>& lib, std::string const& name, experiments::experiment<std::string, uint32_t, uint32_t, float, uint32_t, uint32_t, float, float>& exp, float& size_avg, float& depth_avg )
{
  // mockturtle::mig_npn_resynthesis mig_resyn{true};
  // mockturtle::xag_npn_resynthesis<mockturtle::xag_network> xag_resyn;

  mockturtle::depth_view imig_d{imig};
  printf( "###################################################\n");
  printf( "[i] read_benchmark %s\n", name.c_str() );

  printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
          imig.num_pis(), imig.num_pos(), imig.num_gates(), imig.size(), imig_d.depth() );

  Ntk mig;
  // mockturtle::xag_network res;
  mig = cleanup_dangling( imig );
  float time_i = 0;

  // auto best_size = mig.size();
  mockturtle::depth_view mig_d_tmp{mig};
  // auto best_depth = mig_d_tmp.depth();
  mockturtle::functional_reduction_params frp;
  // mockturtle::functional_reduction_stats st;
  frp.compute_equivalence_classes = true;
  // auto eqpairs = mockturtle::functional_reduction_choices( mig, frp, &st );
  // mockturtle::functional_reduction( mig, frp, &st );

  // mockturtle::choice_view cmig{mig};
  // mockturtle::reduce_choice_network( cmig, eqpairs );
  // mockturtle::improve_representatives( cmig );
  // mig = cleanup_dangling( mig );
  // mockturtle::choice_view<mockturtle::mig_network> cmig2 = mockturtle::levelize_choice_network( cmig );

  mockturtle::map_params ps;
  ps.verbose = true;
  ps.skip_delay_round = false;
  mockturtle::map_stats mst;
  // auto res = mockturtle::map_choices<mockturtle::choice_view<mockturtle::mig_network>, mockturtle::mig_network, mockturtle::mig_npn_resynthesis>( cmig2, mig_resyn, ps, &mst );
  // auto res = mockturtle::map( mig, mig_resyn, ps, &mst );
  auto res = mockturtle::tech_map( mig, lib, ps, &mst );
  // auto res = mockturtle::map<mockturtle::mig_network, mockturtle::xag_network, mockturtle::xag_npn_resynthesis<mockturtle::xag_network>>( mig, xag_resyn, ps, &mst );
  time_i += mockturtle::to_seconds( mst.time_total );

  // mockturtle::functional_reduction( res, frp, &st );
  // res= cleanup_dangling( res );

  // mockturtle::depth_view res_d_tmp{res};

  // if ( res.size() >= best_size )
  // if ( res.size() >= best_size && res_d_tmp.depth() >= best_depth )
  // mig = res;

  mockturtle::depth_view res_d{res};
  // mockturtle::depth_view migr_d{migr};
  printf( "[i] RES: i/o = %d / %d n = %d / %d depth = %d\n",
          res.num_pis(), res.num_pos(), res.num_gates(), res.size(), res_d.depth() );

  float size_impr = ( ( ( (float) imig.num_gates() ) - res.num_gates() ) ) / ( (float) imig.num_gates() ) * 100;
  float depth_impr = ( ( ( (float) imig_d.depth() ) - res_d.depth() ) ) / ( (float) imig_d.depth() ) * 100;
  // uint32_t time = static_cast<uint32_t>( time_i );
  // uint32_t time_cut = static_cast<uint32_t>( mockturtle::to_seconds( pst.time_total ) );

  size_avg += size_impr;
  depth_avg += depth_impr;

  // mockturtle::write_verilog( mig, "itest.v" );
  // mockturtle::write_verilog( res, "test.v" );

  // auto result = abc_cec_benchmark( res, name );
  // assert( result );
  // std::cout << result << std::endl;

  exp( name, imig.num_gates(), res.num_gates(), size_impr, imig_d.depth(), res_d.depth(), depth_impr, time_i );
}


void map()
{
  experiments::experiment<std::string, uint32_t, uint32_t, float, uint32_t, uint32_t, float, float> exp( "Mapper Comparison", "benchmark", "size MIG", "Size Map MIG", "Impr. Size",
                                                                                                            "depth MIG", "depth Map MIG", "Impr. depth", "Map Time (s)" );

  float size_avg = 0.0f, depth_avg = 0.0f;
  auto i = 0u;

  mockturtle::mig_npn_resynthesis mig_resyn{true};

  mockturtle::exact_library<mockturtle::mig_network, mockturtle::mig_npn_resynthesis, 4> lib( mig_resyn );

  for ( const auto& b : local_benchmarks )
  {
    std::string filename{"../test/assets/"};
    filename = filename + b + ".v";
    mockturtle::mig_network imig;
    if ( lorina::read_verilog( filename, mockturtle::verilog_reader( imig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR IN" << std::endl;
      std::abort();
      return;
    }
    map_core( imig, lib, b, exp, size_avg, depth_avg );
    i++;
  }
  // for ( const auto& b : local_benchmarks_iwls )
  // {
  //   std::string filename{"../test/assets/"};
  //   filename = filename + b + ".aig";
  //   mockturtle::mig_network imig;
  //   if ( lorina::read_aiger( filename, mockturtle::aiger_reader( imig ) ) != lorina::return_code::success )
  //   {
  //     std::cout << "ERROR IN" << std::endl;
  //     std::abort();
  //     return;
  //   }
  //   map_core( imig, lib, b, exp, size_avg, depth_avg );
  //   i++;
  // }
  exp.save();
  exp.table();
  printf( "Size avg: %.2f; Depth avg: %.2f\n", size_avg / i, depth_avg / i ); 
}

void tech_map()
{
  std::vector<mockturtle::gate> gates;
  // std::string const file {
  //   "GATE zero 0 O=0;\n"
  //   "GATE one 0 O=1;\n"
  //   "GATE inverter 1 O=!a; PIN * INV 1 999 1.0 1.0 1.0 1.0\n"
  //   "GATE buffer 2 O=a; PIN * NONINV 1 999 1.0 1.0 1.0 1.0\n"
  //   "GATE nand2 1.5 O=!(ab); PIN * NONINV 1 999 1.0 1.0 1.0 1.0\n"
  //   // "GATE or 4 O={ab}; PIN * NONINV 1 999 1.0 1.0 1.0 1.0\n"
  //   "GATE maj3 2.5 O=<abc>; PIN * NONINV 1 999 1.0 1.0 1.0 1.0 1.0\n"
  //   "GATE xor2 4 O=[ab]; PIN * NONINV 1 999 1.0 1.0 1.0 1.0\n"
  //   "GATE xor3 5 O=[abc]; PIN * NONINV 1 999 1.0 1.0 1.0 1.0\n"
  // };

  std::ifstream in( "../../smaller.genlib" );
  // std::istringstream in( file );
  if ( lorina::read_genlib( in, mockturtle::genlib_reader( gates ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR IN" << std::endl;
    std::abort();
    return;
  }
  mockturtle::tech_library<5> lib( gates );

  /* map to library */
  for ( const auto& b : local_benchmarks )
  {
    std::string filename{"../test/assets/"};
    filename = filename + b + ".v";
    mockturtle::aig_network inet;
    if ( lorina::read_verilog( filename, mockturtle::verilog_reader( inet ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR IN" << std::endl;
      std::abort();
      return;
    }
  // for ( const auto& b : local_benchmarks_iwls )
  // {
  //   std::string filename{"../test/assets/"};
  //   filename = filename + b + ".aig";
  //   mockturtle::aig_network inet;
  //   if ( lorina::read_aiger( filename, mockturtle::aiger_reader( inet ) ) != lorina::return_code::success )
  //   {
  //     std::cout << "ERROR IN" << std::endl;
  //     std::abort();
  //     return;
  //   }
    mockturtle::depth_view inet_d{inet};
    printf( "###################################################\n");
    printf( "[i] read_benchmark %s\n", b.c_str() );
    printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
            inet.num_pis(), inet.num_pos(), inet.num_gates(), inet.size(), inet_d.depth() );

    mockturtle::aig_network net;
    net = cleanup_dangling( inet );
    mockturtle::map_params ps;
    ps.cut_enumeration_ps.cut_size = lib.max_gate_size();
    ps.cut_enumeration_ps.cut_limit = 15;
    ps.verbose = true;
    ps.skip_delay_round = false;
    ps.area_flow_rounds = 1;
    // ps.required_time = 12000;
    ps.ela_rounds = 2;
    mockturtle::map_stats mst;

    auto res = mockturtle::tech_mapping( net, lib, ps, &mst );

    mockturtle::depth_view res_d{res};
    printf( "[i] KLUT: i/o = %d / %d n = %d / %d depth = %d\n",
            res.num_pis(), res.num_pos(), res.num_gates(), res.size(), res_d.depth() );

    // auto result = abc_cec_benchmark( res, filename );
    // assert( result );
    // std::cout << result << std::endl;
  }
}


int main()
{
  // map();
  tech_map();
  return 0;
}
