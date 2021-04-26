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
#include <mockturtle/algorithms/cleanup.hpp>

#include <mockturtle/properties/aqfpcost.hpp>

#include "../experiments/experiments.hpp"

#include <fmt/format.h>
#include <lorina/lorina.hpp>

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

template<class Ntk>
bool abc_cec_with_path( const Ntk& ntk, std::string benchmark_path )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/test.bench\"", benchmark_path );

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
    mockturtle::aqfp_db<>& db3,
    mockturtle::aqfp_db<>& db5,
    mockturtle::aqfp_node_resyn_strategy strategy,
    uint32_t iterations,
    std::string lutmap,
    bool pi_buffers,
    bool pi_splitters,
    bool po_buffers )
{
  mockturtle::aqfp_network_cost cost_fn( gate_costs, splitters, pi_buffers, pi_splitters, po_buffers );
  mockturtle::aqfp_node_resyn node_resyn_3( db3, { splitters, strategy, pi_splitters } );
  mockturtle::aqfp_node_resyn node_resyn_5( db5, { splitters, strategy, pi_splitters } );

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

  fmt::print( "processing benchmark {} type {}\n", benchmark_name, ( is_verilog ? "verilog" : "aiger" ) );
  fmt::print( "\tpi: {:4d} po: {:4d} size: {:6d}\n", mig.num_pis(), mig.num_pos(), mig.num_gates() );

  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "iter: " << std::setw(2) << 1 << " " << std::flush;

  mockturtle::klut_network klut_orig = lut_map( mig, lutmap, 4, benchmark_name );

  mockturtle::aqfp_network opt_aqfp;
  mockturtle::aqfp_network opt_aqfp5;
  auto res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_orig, node_resyn_3, fanout_resyn );
  auto res5 = mockturtle::aqfp_resynthesis( opt_aqfp5, klut_orig, node_resyn_5, fanout_resyn );

  auto maj_counts = count_majorities( opt_aqfp5 );
  std::pair<double, uint32_t> res_orig = { cost_fn( opt_aqfp5, res5.node_level, res5.po_level ), res5.critical_po_level() };

  auto res_opt = res_orig;

  for ( auto i = 2u; i <= iterations; i++ )
  {
    std::cout << "\b\b\b" << std::setw(2) << i << " " << std::flush;

    auto klut_opt = lut_map( opt_aqfp, lutmap, 4, benchmark_name );

    opt_aqfp = mockturtle::aqfp_network();
    opt_aqfp5 = mockturtle::aqfp_network();
    res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_opt, node_resyn_3, fanout_resyn );
    res5 = mockturtle::aqfp_resynthesis( opt_aqfp5, klut_opt, node_resyn_5, fanout_resyn );
    std::pair<double, uint32_t> res_temp = { cost_fn( opt_aqfp5, res5.node_level, res5.po_level ), res5.critical_po_level() };

    if ( strategy == mockturtle::aqfp_node_resyn_strategy::cost_based )
    {
      if ( has_better_cost( res_temp, res_opt ) )
      {
        res_opt = res_temp;
        maj_counts = count_majorities( opt_aqfp5 );
      }
    }
    else
    {
      assert( strategy == mockturtle::aqfp_node_resyn_strategy::level_based );
      if ( has_better_level( res_temp, res_opt ) )
      {
        res_opt = res_temp;
        maj_counts = count_majorities( opt_aqfp5 );
      }
    }
  }
  std::cout << "\n";

  auto t2 = std::chrono::high_resolution_clock::now();

  bool cec = abc_cec_with_path( opt_aqfp5, benchmark_path );
  //bool cec = experiments::abc_cec_with_path( opt_aqfp5, benchmark_path );

  auto t3 = std::chrono::high_resolution_clock::now();

  auto exp_time = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
  auto ver_time = std::chrono::duration_cast<std::chrono::milliseconds>( t3 - t2 ).count();

  exp( benchmark_name, (uint32_t)res_opt.first, res_opt.second, maj_counts[3], maj_counts[5], exp_time / 1000.0, ver_time / 1000.0, cec );
}

std::map<mockturtle::aqfp_node_resyn_strategy, std::string> strategy_name = {
    { mockturtle::aqfp_node_resyn_strategy::cost_based, "cost" },
    { mockturtle::aqfp_node_resyn_strategy::level_based, "level" },
};

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  mockturtle::aqfp_db<> db3( gate_costs, splitters );
  mockturtle::aqfp_db<> db5( gate_costs, splitters );

  if ( argc < 4 )
  {
    std::cerr << "Not enough arguments: Need db3 db5 benchmark_name\n";
  }

  std::ifstream db_file3( ( argc < 2 ) ? std::string( "db1.txt" ) : std::string( argv[1] ) );
  std::ifstream db_file5( ( argc < 3 ) ? std::string( "db12.txt" ) : std::string( argv[2] ) );

  assert( db_file3.is_open() );
  assert( db_file5.is_open() );
  db3.load_db_from_file( db_file3 );
  db5.load_db_from_file( db_file5 );
  db_file3.close();
  db_file5.close();

  std::vector<mockturtle::aqfp_node_resyn_strategy> strategies = {
      mockturtle::aqfp_node_resyn_strategy::cost_based,
      mockturtle::aqfp_node_resyn_strategy::level_based,
  };

  std::vector<size_t> iterations = {
      1,
      10,
  };

  std::vector<std::string> lutmaps = {
      "new",
      "new-a",
      "old",
      "old-a",
  };

  std::vector<std::tuple<bool, bool, bool>> configs = {
      { false, false, true },
//      { false, true, true },
//      { false, true, false },
  };

  std::vector<std::string> benchmarks;

  if ( argc > 3 )
  {
    std::string b( argv[3] );
    bool is_verilog = std::string( argv[4] ) == "verilog";
    if ( is_verilog )
    {
      benchmarks.push_back( fmt::format( "./benchmarks/{}.v", b ) );
    }
    else
    {
      benchmarks.push_back( experiments::benchmark_path( b ) );
    }
  }
  else
  {
    for ( auto b : mcnc )
    {
      benchmarks.push_back( fmt::format( "./benchmarks/{}.v", b ) );
    }
  }

  for ( auto s : strategies )
  {
    for ( auto iter : iterations )
    {
      for ( auto lm : lutmaps )
      {
        for ( auto [pib, pis, pob] : configs )
        {
          std::string exp_name = fmt::format( "aqfp_resyn strategy={} iter={} lutmap={} pi_buffers={} pi_splitters={} po_buffers={}", strategy_name[s], iter, lm, pib, pis, pob );
          experiments::experiment<std::string, uint32_t, uint32_t, uint32_t, uint32_t, double, double, bool> exp( exp_name, "benchmark", "JJ count", "JJ level", "maj 3 count", "maj 5 count", "resyn time", "verify time", "cec" );
          fmt::print( "\n\n\nexperiment: {}\n", exp_name );
          for ( auto path : benchmarks )
          {
            do_experiment( exp, path, gate_costs, splitters, db3, db5, s, iter, lm, pib, pis, pob );
            // exp.save();
            // exp.table();
          }
          exp.save();
          exp.table();
        }
      }
    }
  }

  // std::ofstream db_stats( "stats.csv" );
  // db3.print_usage_state( db_stats );
  // db5.print_usage_state( db_stats );
  // db_stats.close();

  return 0;
}
