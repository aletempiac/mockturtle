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

#include <mockturtle/properties/aqfpcost.hpp>

#include <mockturtle/utils/tech_library.hpp>

#include "../experiments/experiments.hpp"

#include <fmt/format.h>
#include <lorina/lorina.hpp>

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  mockturtle::aqfp_db<> db3( gate_costs, splitters );

  std::ifstream db_file3( std::string( "db1.txt" ) );

  assert( db_file3.is_open() );
  db3.load_db_from_file( db_file3 );
  db_file3.close();

  mockturtle::aqfp_exact_library<mockturtle::aqfp_network, 4u> lib( db3 );

  return 0;
}
