#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <thread>
#include <vector>

#include <mockturtle/algorithms/aqfp_resynthesis/detail/db_builder.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/detail/dag_gen.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/detail/db_utils.hpp>

#include <fmt/format.h>

std::vector<uint32_t> string_to_uint_vec( std::string str )
{
  std::stringstream ss( str );
  std::vector<uint32_t> res;
  uint32_t temp;
  while ( ss >> temp )
  {
    res.push_back( temp );
  }
  return res;
}

std::unordered_map<uint32_t, uint32_t> string_to_uint_uint_map( std::string str )
{
  std::stringstream ss( str );
  std::unordered_map<uint32_t, uint32_t> res;
  uint32_t t1, t2;
  while ( ss >> t1 >> t2 )
  {
    res[t1] = t2;
  }
  return res;
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  if ( argc < 2 )
  {
    std::cerr << fmt::format( "Not enough arguments. Usage: {} cmd [opt]\n", std::string( argv[0] ) );
    return 0;
  }

  auto num_threads = std::thread::hardware_concurrency();
  if ( num_threads == 0u )
    num_threads = 1u;
  std::cerr << fmt::format( "Will be using {} threads\n", num_threads );

  mockturtle::dag_generator_params params;
  params.allowed_num_fanins = { 3u };
  params.max_gates_of_fanin = { { 3u, 7u } };
  params.max_gates = 7u;
  params.max_levels = 7u;
  params.max_num_in = 4u;

  params.verbose = 1u;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  std::string cmd( argv[1] );
  if ( cmd == "generate-dags" )
  {
    if ( argc < 3 )
    {
      std::cerr << fmt::format( "Usage: {} generate-dags dag_file\n", std::string( argv[0] ) );
      return 0;
    }

    std::string dag_file_prefix( argv[2] );

    if ( argc > 3 )
    {
      if ( argc != 8 )
      {
        std::cerr << fmt::format( "Usage: {} generate-dags dag_file allowed_num_fanins max_gates_of_fanin max_gates max_level max_num_in\n", std::string( argv[0] ) );
        return 0;
      }

      params.allowed_num_fanins = string_to_uint_vec( std::string( argv[3] ) );
      params.max_gates_of_fanin = string_to_uint_uint_map( std::string( argv[4] ) );
      params.max_gates = std::stoul( std::string( argv[5] ) );
      params.max_levels = std::stoul( std::string( argv[6] ) );
      params.max_num_in = std::stoul( std::string( argv[7] ) );
    }

    mockturtle::generate_aqfp_dags( params, dag_file_prefix, num_threads );
  }
  else if ( cmd == "compute-costs" )
  {
    if ( argc != 4 )
    {
      std::cerr << fmt::format( "Usage: {} compute-costs dag_file cost_file\n", std::string( argv[0] ) );
      return 0;
    }

    std::string dag_file_prefix( argv[2] );
    std::string cost_file_prefix( argv[3] );

    mockturtle::compute_aqfp_dag_costs( gate_costs, splitters, dag_file_prefix, cost_file_prefix, num_threads );
  }
  else if ( cmd == "generate-db" )
  {
    if ( argc != 5 )
    {
      std::cerr << fmt::format( "Usage: {} generate-db dag_file cost_file db_file\n", std::string( argv[0] ) );
      return 0;
    }

    std::string dag_file_prefix( argv[2] );
    std::string cost_file_prefix( argv[3] );
    std::string db_file_prefix( argv[4] );

    mockturtle::generate_aqfp_db( gate_costs, splitters, dag_file_prefix, cost_file_prefix, db_file_prefix, num_threads );
  }
  else if ( cmd == "db-from-scratch" )
  {
    if ( argc < 3 )
    {
      std::cerr << fmt::format( "Usage: {} db-from-scratch file_prefix\n", std::string( argv[0] ) );
      return 0;
    }

    std::string file_prefix( argv[2] );

    if ( argc > 3 )
    {
      if ( argc != 8 )
      {
        std::cerr << fmt::format( "Usage: {} db-from-scratchs file_prefix allowed_num_fanins max_gates_of_fanin max_gates max_level max_num_in\n", std::string( argv[0] ) );
        return 0;
      }

      params.allowed_num_fanins = string_to_uint_vec( std::string( argv[3] ) );
      params.max_gates_of_fanin = string_to_uint_uint_map( std::string( argv[4] ) );
      params.max_gates = std::stoul( std::string( argv[5] ) );
      params.max_levels = std::stoul( std::string( argv[6] ) );
      params.max_num_in = std::stoul( std::string( argv[7] ) );
    }

    mockturtle::generate_aqfp_db( params, gate_costs, splitters, file_prefix, num_threads );
  }
  else if ( cmd == "db-merge" )
  {
    if ( argc != 5 )
    {
      std::cerr << fmt::format( "Usage: {} db-merge input_file_1 input_file_2 output_file\n", std::string( argv[0] ) );
      return 0;
    }

    std::string f1( argv[2] );
    std::string f2( argv[3] );
    std::string fo( argv[4] );

    mockturtle::aqfp_db_builder<> builder;
    std::ifstream is1( f1 );
    std::ifstream is2( f2 );
    assert( is1.is_open() && is2.is_open() );
    builder.load_db_from_file( is1 );
    builder.load_db_from_file( is2 );
    is1.close();
    is2.close();

    builder.remove_redundant();

    std::ofstream os( fo );
    assert( os.is_open() );
    builder.save_db_to_file( os );
    os.close();
  }
  else
  {
    std::cerr << fmt::format( "Invalid command {}. Must be one of the following:\n"
                              "\tgenerate-dags   -- for generating DAGs\n"
                              "\tcompute-costs   -- for costing DAGs\n"
                              "\tgenerate-db     -- for generating the AQFP database\n"
                              "\tdb-from-scratch -- for generating the AQFP database from scratch\n"
                              "\tdb-merge        -- for merging two databases\n",
                              cmd );
  }

  return 0;
}
