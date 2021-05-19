#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>

#include <mockturtle/networks/aqfp.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>

#include <mockturtle/algorithms/aqfp_resynthesis.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_db.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_fanout_resyn.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_node_resyn.hpp>

#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>

#include <mockturtle/algorithms/aqfp_mapper.hpp>
#include <mockturtle/algorithms/mapper.hpp>

#include <mockturtle/properties/aqfpcost.hpp>

#include <mockturtle/utils/tech_library.hpp>

#include "../experiments/experiments.hpp"

#include <fmt/format.h>
#include <lorina/lorina.hpp>

std::vector<std::string> mcnc = {
    "5xp1",
    "c1908",
    "c432",
    "c5315",
    "c880",
    "chkn",
    "count",
    "dist",
    "in5",
    "in6",
    "k2",
    "m3",
    "max512",
    "misex3",
    "mlp4",
    "prom2",
    "sqr6",
    "x1dn",
};

std::string mcnc_path( std::string const& benchmark_name )
{
  #ifndef EXPERIMENTS_PATH
    return fmt::format( "{}.v", benchmark_name );
  #else
    return fmt::format( "{}benchmarks_aqfp/{}.v", EXPERIMENTS_PATH, benchmark_name );
  #endif
}

template<typename Result>
bool has_better_cost( Result& current, Result& previous )
{
  if ( current.first < previous.first )
    return true;

  if ( current.first > previous.first )
    return false;

  return current.second < previous.second;
}

template<typename Result>
bool has_better_level( Result& current, Result& previous )
{
  if ( current.second < previous.second )
    return true;

  if ( current.second > previous.second )
    return false;

  return current.first < previous.first;
}

std::map<mockturtle::aqfp_node_resyn_strategy, std::string> strategy_name = {
    { mockturtle::aqfp_node_resyn_strategy::cost_based, "cost" },
    { mockturtle::aqfp_node_resyn_strategy::level_based, "level" },
};

template<class Ntk>
bool abc_cec_with_path( const Ntk& ntk, std::string benchmark_path, std::string benchmark_name )
{
  mockturtle::write_bench( ntk, fmt::format("/tmp/test_{}.bench", benchmark_name ) );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/test_{}.bench\"", benchmark_path, benchmark_name );

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

  if ( result.size() >= 23 && result.substr( 0u, 23u ) == "Networks are equivalent" )
  {
    return true;
  }
  return false;
}

template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, std::string type, uint32_t k = 4, std::string name = {} )
{
  std::string tempfile1 = "temp1_" + name + ".blif";
  std::string tempfile2 = "temp2_" + name + ".blif";

  mockturtle::write_blif( ntk, tempfile1 );

  if ( type == "new" )
  {
    system( fmt::format( "abc -q \"{}; &get; &if -K {}; &put; write_blif {}\" >> /dev/null 2>&1", tempfile1, k, tempfile2 ).c_str() );
  }
  else if (type == "new-a")
  {
    system( fmt::format( "abc -q \"{}; &get; &if -a -K {}; &put; write_blif {}\" >> /dev/null 2>&1", tempfile1, k, tempfile2 ).c_str() );
  }
  else if (type == "old")
  {
    system( fmt::format( "abc -q \"{}; if -K {}; write_blif {}\" >> /dev/null 2>&1", tempfile1, k, tempfile2 ).c_str() );
  }
  else if (type == "old-a")
  {
    system( fmt::format( "abc -q \"{}; if -a -K {}; write_blif {}\" >> /dev/null 2>&1", tempfile1, k, tempfile2 ).c_str() );
  }

  mockturtle::klut_network klut;
  if ( lorina::read_blif( tempfile2, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "FATAL NEW LUT MAP - Reading mapped network failed! " << tempfile1 << " " << tempfile2 << std::endl;
    std::abort();
    return klut;
  }

  system( fmt::format( "rm {}", tempfile1 ).c_str() );
  system( fmt::format( "rm {}", tempfile2 ).c_str() );
  //klut = mockturtle::cleanup_dangling(klut);
  return klut;
}

template<typename T>
auto count_majorities( T& ntk )
{
  std::unordered_map<uint32_t, uint32_t> counts;
  ntk.foreach_gate( [&]( auto n ) { counts[ntk.fanin_size( n )]++; } );
  return counts;
}

template<typename ExpT>
void do_experiment(
    ExpT& exp,
    std::string benchmark_path,
    std::unordered_map<uint32_t, double>& gate_costs,
    std::unordered_map<uint32_t, double>& splitters,
    mockturtle::aqfp_db<>& db,
    mockturtle::aqfp_node_resyn_strategy strategy,
    uint32_t iterations,
    std::string lutmap,
    bool pi_buffers,
    bool pi_splitters,
    bool po_buffers,
    mockturtle::aqfp_exact_library<mockturtle::aqfp_network, 4u>& lib,
    double& avg_jjs_size_impr,
    double& avg_jjs_lev_impr )
{
  mockturtle::aqfp_network_cost cost_fn( gate_costs, splitters, pi_buffers, pi_splitters, po_buffers );
  mockturtle::aqfp_node_resyn node_resyn( db, { splitters, strategy, pi_splitters } );

  uint32_t max_branching_factor = std::max_element( splitters.begin(), splitters.end(), [&]( auto s1, auto s2 ) { return s1.first < s2.first; } )->first;
  mockturtle::aqfp_fanout_resyn fanout_resyn( max_branching_factor, pi_splitters );

  auto last_dir_delim_ind = benchmark_path.find_last_of( "/" );
  std::string benchmark_name = ( last_dir_delim_ind == std::string::npos ) ? benchmark_path : benchmark_path.substr( last_dir_delim_ind + 1 );

  auto last_ext_delim_ind = benchmark_name.find_last_of( "." );
  assert( last_ext_delim_ind != std::string::npos );

  bool is_verilog = ( benchmark_name.substr( last_ext_delim_ind + 1 ) == "v" );
  if ( !is_verilog )
  {
    assert( benchmark_name.substr( last_ext_delim_ind + 1 ) == "aig" );
  }
  benchmark_name = benchmark_name.substr( 0, last_ext_delim_ind );

  mockturtle::mig_network mig;
  if ( is_verilog )
  {
    lorina::read_verilog( benchmark_path, mockturtle::verilog_reader( mig ) );
  }
  else
  {
    lorina::read_aiger( benchmark_path, mockturtle::aiger_reader( mig ) );
  }

  mockturtle::map_params ps;
  // ps.skip_delay_round = true;
  // ps.required_time = std::numeric_limits<float>::max();
  mockturtle::map_stats st;

  fmt::print( "processing benchmark {} type {}\n", benchmark_name, ( is_verilog ? "verilog" : "aiger" ) );
  fmt::print( "\tpi: {:4d} po: {:4d} size: {:6d}\n", mig.num_pis(), mig.num_pos(), mig.num_gates() );

  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "iter: " << std::setw(2) << 1 << " " << std::flush;

  mockturtle::klut_network klut_orig_lut = lut_map( mig, lutmap, 4, benchmark_name );
  mockturtle::klut_network klut_orig_map = mockturtle::map_aqfp( mig, lib, ps, &st );

  mockturtle::aqfp_network opt_aqfp_lut, opt_aqfp_map;
  auto res_lut = mockturtle::aqfp_resynthesis( opt_aqfp_lut, klut_orig_lut, node_resyn, fanout_resyn );
  auto res_map = mockturtle::aqfp_resynthesis( opt_aqfp_map, klut_orig_map, node_resyn, fanout_resyn );

  auto maj_counts_lut = count_majorities( opt_aqfp_lut );
  std::pair<double, uint32_t> res_orig_lut = { cost_fn( opt_aqfp_lut, res_lut.node_level, res_lut.po_level ), res_lut.critical_po_level() };
  auto maj_counts_map = count_majorities( opt_aqfp_map );
  std::pair<double, uint32_t> res_orig_map = { cost_fn( opt_aqfp_map, res_map.node_level, res_map.po_level ), res_map.critical_po_level() };

  auto res_opt_lut = res_orig_lut;
  auto res_opt_map = res_orig_map;

  for ( auto i = 2u; i <= iterations; i++ )
  {
    std::cout << "\b\b\b" << std::setw(2) << i << " " << std::flush;

    auto klut_opt_lut = lut_map( opt_aqfp_lut, lutmap, 4, benchmark_name );
    auto klut_opt_map = mockturtle::map_aqfp( opt_aqfp_map, lib, ps, &st );

    opt_aqfp_lut = mockturtle::aqfp_network();
    res_lut = mockturtle::aqfp_resynthesis( opt_aqfp_lut, klut_opt_lut, node_resyn, fanout_resyn );
    std::pair<double, uint32_t> res_temp_lut = { cost_fn( opt_aqfp_lut, res_lut.node_level, res_lut.po_level ), res_lut.critical_po_level() };

    opt_aqfp_map = mockturtle::aqfp_network();
    res_map = mockturtle::aqfp_resynthesis( opt_aqfp_map, klut_opt_map, node_resyn, fanout_resyn );
    std::pair<double, uint32_t> res_temp_map = { cost_fn( opt_aqfp_map, res_map.node_level, res_map.po_level ), res_map.critical_po_level() };


    if ( strategy == mockturtle::aqfp_node_resyn_strategy::cost_based )
    {
      if ( has_better_cost( res_temp_lut, res_opt_lut ) )
      {
        res_opt_lut = res_temp_lut;
        maj_counts_lut = count_majorities( opt_aqfp_lut );
      }
      if ( has_better_cost( res_temp_map, res_opt_map ) )
      {
        res_opt_map = res_temp_map;
        maj_counts_map = count_majorities( opt_aqfp_map );
      }
    }
    else
    {
      assert( strategy == mockturtle::aqfp_node_resyn_strategy::level_based );
      if ( has_better_level( res_temp_lut, res_opt_lut ) )
      {
        res_opt_lut = res_temp_lut;
        maj_counts_lut = count_majorities( opt_aqfp_lut );
      }
      if ( has_better_level( res_temp_map, res_opt_map ) )
      {
        res_opt_map = res_temp_map;
        maj_counts_map = count_majorities( opt_aqfp_map );
      }
    }
  }
  std::cout << "\n";

  // auto t2 = std::chrono::high_resolution_clock::now();

  bool cec1 = abc_cec_with_path( opt_aqfp_lut, benchmark_path, benchmark_name );
  bool cec2 = abc_cec_with_path( opt_aqfp_map, benchmark_path, benchmark_name );

  // auto t3 = std::chrono::high_resolution_clock::now();

  // auto exp_time = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  // auto ver_time = std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 ).count();

  double jj_size_impr = ( (double)res_opt_lut.first - (double)res_opt_map.first ) / (double)res_opt_lut.first * 100;
  double jj_level_impr = ( (double)res_opt_lut.second - (double)res_opt_map.second ) / (double)res_opt_lut.second * 100;

  avg_jjs_size_impr += jj_size_impr;
  avg_jjs_lev_impr += jj_level_impr;

  exp( benchmark_name, (uint32_t)res_opt_lut.first, (uint32_t)res_opt_map.first, jj_size_impr, res_opt_lut.second, res_opt_map.second, jj_level_impr, cec1, cec2 );
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  using namespace experiments;
  using namespace mockturtle;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  aqfp_db<> db3( gate_costs, splitters );

  std::ifstream db_file3( std::string( "../experiments/db3.txt" ) );

  assert( db_file3.is_open() );
  db3.load_db_from_file( db_file3 );
  db_file3.close();

  /* create a tech library with the AQFP maj3 database */
  aqfp_exact_library<aqfp_network, 4u> lib( db3 );

  experiment<std::string, uint32_t, uint32_t, double, uint32_t, uint32_t, double, bool, bool> exp(
      "aqfp_mapper", "benchmark", "lut map JJs", "aqfp map JJs", "JJs improvement", "lut map JJ level", "aqfp map JJ level", "JJ levels impr", "cec1", "cec2" );

  /* params */
  unsigned iter = 10u;
  std::string lm( "new" );
  bool pib = false;
  bool pis = true;
  bool pob = true;

  double avg_jjs_size_impr = 0;
  double avg_jjs_lev_impr = 0;

  for ( auto const& benchmark : mcnc )
  {
    do_experiment( exp, mcnc_path( benchmark ), gate_costs, splitters, db3, aqfp_node_resyn_strategy::cost_based, iter, lm, pib, pis, pob, lib, avg_jjs_size_impr, avg_jjs_lev_impr );
  }

  exp.save();
  exp.table();

  std::cout << "AVG JJ size improvement :" << std::setw(5) << avg_jjs_size_impr / mcnc.size() << "\n";
  std::cout << "AVG JJ level improvement:" << std::setw(5) << avg_jjs_lev_impr / mcnc.size() << "\n";

  return 0;
}
