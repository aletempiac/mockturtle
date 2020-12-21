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
#include <mockturtle/algorithms/detail/database_generator.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig4_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/functional_reduction.hpp>
// #include <mockturtle/views/fanout_limit_view.hpp>
#include <mockturtle/views/topo_view.hpp>
#include <mockturtle/views/choice_view.hpp>
#include <mockturtle/utils/choice_utils.hpp>


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
  "square",
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
  std::string command = fmt::format( "../../abc/abc -q \"cec -n {} /tmp/test.bench\"", benchmark );

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

  return result.size() >= 23 && result.substr( 0u, 23u ) == "Networks are equivalent";
}


template<typename Ntk>
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
  mockturtle::write_verilog( ntk, "/tmp/network.v" );
  system( fmt::format( "../../abc/abc -q \"/tmp/network.v; &get; &if -a -K {}; &put; write_blif /tmp/output.blif\"", k ).c_str() );
  mockturtle::klut_network klut;
  if ( lorina::read_blif( "/tmp/output.blif", mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "ERROR LUT" << std::endl;
    std::abort();
    return klut;
  }
  return klut;
}


void create_database()
{
  /* enumerate NPN representatives */
  std::unordered_set<kitty::dynamic_truth_table, kitty::hash<kitty::dynamic_truth_table>> classes;
  kitty::dynamic_truth_table tt( 4u );
  auto i = 0;
  do
  {
    i++;
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  std::cout << "[i] enumerated "
            << ( 1 << ( 1 << tt.num_vars() ) ) << " functions into "
            << classes.size() << " classes." << std::endl;

  /* generate database with exact MIG synthesis */
  mockturtle::mig_network mig;
  mockturtle::exact_mig_resynthesis_params ps;
  ps.num_candidates = 4u;
  mockturtle::detail::database_generator_params ps_db;
  mockturtle::exact_mig_resynthesis<mockturtle::mig_network> exact( ps );
  ps_db.verbose = true;
  ps_db.multiple_candidates = true;
  mockturtle::detail::database_generator dbgen( mig, exact, ps_db );
  for ( const auto& f : classes )
  {
    dbgen.add_function( f );

    std::cout << ".";
    std::cout.flush();
  }
  mockturtle::write_verilog( mig, "db.v" );
}


void synthesis()
{
  mockturtle::mig_network mig_db;
  read_verilog( "db.v", mockturtle::verilog_reader( mig_db ) );

  mockturtle::mig4_npn_resynthesis_params ps;
  //ps.multiple_depth = true;
  mockturtle::mig4_npn_resynthesis<mockturtle::mig_network> mig_resyn( mockturtle::detail::to_index_list( mig_db ), ps );

  auto counter = 0;
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

    //write_verilog( imig, "out_" + b + ".v" );

    mockturtle::depth_view imig_d{imig};
    printf( "###################################################\n");
    printf( "[i] read_benchmark %s\n", filename.c_str() );

    printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
            imig.num_pis(), imig.num_pos(), imig.size() - imig.num_pis() - 1, imig.size(), imig_d.depth() );

    // auto klut = lut_map( imig, 4u );

    mockturtle::mig_network mig;
    //mig = mockturtle::node_resynthesis<mockturtle::mig_network>( klut, mig_resyn );
    mig = cleanup_dangling( imig );

    mockturtle::functional_reduction_params frp;
    mockturtle::functional_reduction_stats st;
    frp.compute_equivalence_classes = true;

    mockturtle::choice_view_params cps;

    /* compute the equivalence classes and substitute representatives in the net */
    // auto eqclasses = mockturtle::functional_reduction_eqclasses( mig, frp, &st );
    // mockturtle::choice_view cmig{mig};
    //eqclasses.print_eqclasses();
    //mockturtle::functional_reduction( mig, frp, &st );
    //mig = cleanup_dangling( mig );
    //std::cout << "size: " << mig.size() << " ";
    //mockturtle::choice_view<mockturtle::mig_network> cv{mig};
    // mockturtle::choice_network<mockturtle::mig_network>( cmig, eqclasses );
    //std::cout << mig.size() << std::endl;
    //eqclasses.print_eqclasses();

    mockturtle::cut_rewriting_params psc;
    psc.cut_enumeration_ps.cut_size = 4;

    int i = 0;
    while ( i++ < 10 )
    {
      mockturtle::mig_network new_mig;
      auto const mig_gates_before = mig.num_gates();

      auto eqpairs = mockturtle::functional_reduction_eqclasses( mig, frp, &st );
      // mig = cleanup_dangling( mig );
      mockturtle::choice_view cmig( mig, cps );
      mockturtle::reduce_choice_network( cmig, eqpairs );
      mockturtle::choice_view<mockturtle::mig_network> cmig2 = mockturtle::cleanup_choice_network( cmig );
      // cmig.print_choice_classes();

      //mockturtle::topo_view mig_topo{ mig };
      //new_mig = mockturtle::cut_rewriting_eqclasses<mockturtle::mig_network>( mig, eqclasses, mig_resyn, psc );
      //new_mig = mockturtle::cut_rewriting<mockturtle::mig_network>( mig_topo, mig_resyn, psc );
      // std::cout << "to cut-rewriting" << std::endl;
      new_mig = mockturtle::cut_rewriting_area_flow<mockturtle::mig_network>( cmig2, mig_resyn, psc );
      new_mig = cleanup_dangling( new_mig );

      if ( new_mig.num_gates() > mig.num_gates() ) {
        new_mig = cleanup_dangling( mig );
      }
      std::cout << "i: " << i << "; gates size " << new_mig.num_gates() << "/" << mig.num_gates() << std::endl;
      //eqclasses.print_eqclasses();
      //mockturtle::functional_reduction( new_mig, frp, &st );

      //auto clean_new_mig = cleanup_dangling( new_mig );
      //new_mig = cleanup_dangling( new_mig );

      //mockturtle::depth_view depth_mig_new{new_mig};

      if ( new_mig.num_gates() >= mig_gates_before ) {
        break;
      }
      mig = new_mig;
      std::cout << i;
    }
    std::cout << std::endl;
    mig = cleanup_dangling( mig );
    //mockturtle::functional_reduction( mig, frp, &st );
    //mig = cleanup_dangling( mig );

    mockturtle::depth_view mig_d{mig};

    printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
            mig.num_pis(), mig.num_pos(), mig.size() - mig.num_pis() - 1, mig.size(), mig_d.depth() );


    //klut = lut_map( mig, 6u );

    //std::string const out_filename = "out_" + b + ".blif";
    //write_blif( klut, out_filename );

    /* abc ceq */
    //if ( counter > 3 ) {
      //auto result = abc_cec_benchmark( mig, filename );
      //assert( result );
    //}

    /* run ABC print_stats for double checking and to see depth information for the mapped network */
    //std::string const shell_command = "../../abc/abc -c \"" + out_filename + "; print_stats;\"";
    //system( shell_command.c_str() );
    counter++;
  }
}


void synthesis_iwls()
{
  mockturtle::mig_network mig_db;
  read_verilog( "db.v", mockturtle::verilog_reader( mig_db ) );

  mockturtle::mig4_npn_resynthesis_params ps;
  //ps.multiple_depth = true;
  mockturtle::mig4_npn_resynthesis<mockturtle::mig_network> mig_resyn( mockturtle::detail::to_index_list( mig_db ), ps );

  auto counter = 0;
  for ( const auto& b : local_benchmarks_iwls )
  {
    std::string filename{"../test/assets/"};
    filename = filename + b + ".aig";

    mockturtle::mig_network imig;
    if ( lorina::read_aiger( filename, mockturtle::aiger_reader( imig ) ) != lorina::return_code::success )
    {
      std::cout << "ERROR IN" << std::endl;
      std::abort();
      return;
    }

    //write_verilog( imig, "out_" + b + ".v" );

    mockturtle::depth_view imig_d{imig};
    printf( "###################################################\n");
    printf( "[i] read_benchmark %s\n", filename.c_str() );

    printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
            imig.num_pis(), imig.num_pos(), imig.size() - imig.num_pis() - 1, imig.size(), imig_d.depth() );

    mockturtle::mig_network mig;
    mig = cleanup_dangling( imig );

    mockturtle::functional_reduction_params frp;
    mockturtle::functional_reduction_stats st;
    frp.compute_equivalence_classes = true;
    //eqclasses_init.print_eqclasses();
    //mig = cleanup_dangling( mig );

    mockturtle::cut_rewriting_params psc;
    psc.cut_enumeration_ps.cut_size = 4;

    // auto eqclasses = mockturtle::functional_reduction_eqclasses( mig, frp, &st );
    //eqclasses.print_eqclasses();
    // mockturtle::choice_network( mig, eqclasses );
    //eqclasses.print_eqclasses();

    int i = 0;
    while ( i++ < 10 )
    {
      mockturtle::mig_network new_mig;
      auto const mig_gates_before = mig.num_gates();

      auto eqpairs = mockturtle::functional_reduction_eqclasses( mig, frp, &st );
      //mig = cleanup_dangling( mig );
      mockturtle::choice_view cmig{mig};
      mockturtle::reduce_choice_network( cmig, eqpairs );
      mockturtle::choice_view<mockturtle::mig_network> cmig2 = mockturtle::cleanup_choice_network( cmig );


      // cmig.print_choice_classes();
      //mockturtle::topo_view mig_topo{ mig };
      //new_mig = mockturtle::cut_rewriting_eqclasses<mockturtle::mig_network>( mig_topo, eqclasses, mig_resyn, psc );
      //new_mig = mockturtle::cut_rewriting<mockturtle::mig_network>( mig_topo, mig_resyn, psc );
      new_mig = mockturtle::cut_rewriting_area_flow<mockturtle::mig_network>( cmig2, mig_resyn, psc );

      if ( new_mig.num_gates() > mig.num_gates() ) {
        new_mig = cleanup_dangling( mig );
      }

      std::cout << "i: " << i << "; gates size " << new_mig.num_gates() << "/" << mig.num_gates() << std::endl;


      // eqclasses = mockturtle::functional_reduction_eqclasses( new_mig, frp, &st );
      // mockturtle::choice_network( new_mig, eqclasses );

      //mockturtle::depth_view depth_mig_new{new_mig};

      if ( new_mig.num_gates() >= mig_gates_before ) {
        break;
      }
      mig = new_mig;
    }

    mockturtle::depth_view mig_d{mig};

    printf( "[i] MIG: i/o = %d / %d n = %d / %d depth = %d\n",
            mig.num_pis(), mig.num_pos(), mig.size() - mig.num_pis() - 1, mig.size(), mig_d.depth() );


    //klut = lut_map( mig, 6u );

    //std::string const out_filename = "out_" + b + ".blif";
    //write_blif( klut, out_filename );

    /* abc ceq */
    //if ( counter > 3 ) {
      //auto result = abc_cec_benchmark( mig, filename );
      //assert( result );
    //}

    /* run ABC print_stats for double checking and to see depth information for the mapped network */
    //std::string const shell_command = "../../abc/abc -c \"" + out_filename + "; print_stats;\"";
    //system( shell_command.c_str() );
    counter++;
  }
}


int main()
{
  //create_database();
  // synthesis();
  synthesis_iwls();
  return 0;
}
