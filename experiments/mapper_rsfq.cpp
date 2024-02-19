/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <lorina/genlib.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_network_conversion.hpp>
#include <mockturtle/algorithms/rsfq/rsfq_path_balancing.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/retiming.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/networks/xmg.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/rsfq_view.hpp>

#include <experiments.hpp>

constexpr int fDFF   = 0;
constexpr int fNOT   = 1;
constexpr int fMERGE = 2;
constexpr int fOR    = 3;
constexpr int fAND   = 4;
constexpr int fXOR   = 5;
constexpr int fOR3   = 6;
constexpr int fAND3  = 7;
constexpr int fMAJ3  = 8;
constexpr int fCB    = 9;
constexpr int fSPL   = 10;
// constexpr int fPI    = 11;
// constexpr int fNOFUNC= 99;

// Removed input buffers in AND/OR gates
constexpr std::array<int,12> COSTS_CONNECT = {6, 9, 7, 3, 3, 11, 11, 11, 11, 7, 3, 0};

typedef mockturtle::klut_network klut;
typedef mockturtle::xag_network xag;
typedef mockturtle::xmg_network xmg;
typedef mockturtle::aig_network aig;

template <typename Ntk>
std::tuple<mockturtle::binding_view<klut>, mockturtle::map_stats, int, int, bool> map_with_pb 
( 
  const std::string & benchmark, 
  const Ntk & tech_indep_ntk, 
  const mockturtle::tech_library<4u, mockturtle::classification_type::p_configurations> & tech_lib, 
  std::unordered_map<std::string, int> & nDFF_global, 
  bool area_oriented = false 
)
{
  mockturtle::map_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = true;
  ps.cut_enumeration_ps.cut_limit = 24;
  ps.buffer_pis = false;
  if (area_oriented)
  {
      ps.skip_delay_round = true;
      ps.required_time = std::numeric_limits<float>::max();
  }
  mockturtle::map_stats st;
  mockturtle::binding_view<klut> res = map( tech_indep_ntk, tech_lib, ps, &st );

  // mockturtle::emap_params ps;
  // ps.cut_enumeration_ps.cut_limit = 24;
  // ps.buffer_pis = false;
  // ps.area_oriented_mapping = area_oriented;
  // mockturtle::emap_stats st;
  // mockturtle::binding_view<klut> res = emap_klut( tech_indep_ntk, tech_lib, ps, &st );
  mockturtle::depth_view<mockturtle::binding_view<klut>> dv { res };

  std::map<klut::node, int> dff_count;
  std::map<klut::node, int> fanout_count;

  /* RSFQ path balancing */
  auto balanced_res = mockturtle::rsfq_path_balancing( res );

  mockturtle::retime_params rps;
  mockturtle::retime_stats rst;
  auto net = mockturtle::rsfq_generic_network_create_from_mapped( balanced_res );
  mockturtle::retime( net, rps, &rst );
  auto retime_res = mockturtle::rsfq_mapped_create_from_generic_network( net );

  uint32_t num_ext_dffs = retime_res.num_dffs();
  
  uint32_t num_int_dffs = 0;
  retime_res.foreach_node( [&]( auto const& n ) {
    if ( !retime_res.has_binding( n ) )
      return;
    auto ndff = nDFF_global[retime_res.get_binding( n ).name];
    // fmt::print("Adding {} internal DFFs from {}\n", ndff,  retime_res.get_binding( n ).name);
    num_int_dffs += ndff;
  } );

  /* RSFQ splitter insertion */
  uint32_t num_splitters = 0;
  retime_res.foreach_node( [&]( auto const& n ) {
    if ( !retime_res.is_constant( n ) )
      num_splitters += retime_res.fanout_size( n ) - 1;
  } );

  bool cec = rsfq_check_buffering( retime_res );
  cec &= benchmark == "hyp" ? true : experiments::abc_cec( retime_res, benchmark );

  // Internal DFF area is already counted in the library
  int total_area = st.area + COSTS_CONNECT[fSPL] * num_splitters +  COSTS_CONNECT[fDFF] * num_ext_dffs;
  // fmt::print("\t{} : Int: {}, Ext: {}, ratio: {}\n", benchmark, num_int_dffs, num_ext_dffs, (float)num_int_dffs / (num_int_dffs + num_ext_dffs) );
  return std::make_tuple( res, st, num_int_dffs + num_ext_dffs, total_area, cec );
}

// Function to read unordered_map from CSV file
std::unordered_map<std::string, int> readCSV(const std::string& filename) 
{
    std::ifstream infile(filename);             // Open the input file stream
    std::unordered_map<std::string, int> map;   // Create the unordered_map
    
    std::string line;
    std::getline(infile, line);                 // Ignore the header row

    fmt::print("READING CSV : {}\n", filename);
    // Read each subsequent row and add the key-value pair to the unordered_map
    while (std::getline(infile, line)) 
    {
        std::stringstream ss(line);
        std::string key;
        int value;
        std::getline(ss, key, ',');
        ss >> value;
        map[key] = value;
        // fmt::print("READ {} : {}\n", key, value);
    }
    infile.close(); // Close the input file stream
    return map;
}

#if false
  // some old libraries

  /* Gate costs are based on CONNECT library (from Japan) */
  // const std::string DATABASE_PATH { "LIBRARY_2023_05_19_CONNECT.genlib" } ;
  /* CONNECT library (from Japan) */
  // const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_VANILLA_CONNECT.genlib" } ;
  // LIBRARY_2023_05_19_CONNECT.genlib
  // const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_06_26_CONNECT_1111.genlib" } ;
  // const std::string DATABASE_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/LIBRARY_2023_05_19_CONNECT.genlib" } ;

  /*The number of internal DFFs within each cell. 
  Some of them are necessary not only for path balancing but also 
  for synchronizing the pulses for AND gates. I include them 
  in total DFF count */
  // const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle_alessandro/build/nDFF_2023_05_08_CONNECT.csv" } ; 
  // const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/nDFF_2023_06_27_CONNECT_CONSERVATIVE.csv" } ; 
  // const std::string NDFF_PATH { "/Users/brainkz/Documents/GitHub/mockturtle/build/NDFF_PARSED_2023_06_27_CONNECT_CONSERVATIVE.csv" } ; //updated costs here

#endif


const std::string DATABASE_2_INPUT_PATH { "rsfq_databases/LIBRARY_2IN.genlib" } ;
const std::string DATABASE_3_INPUT_PATH { "rsfq_databases/LIBRARY_3IN.genlib" } ;

const std::string NDFF_2_INPUT_PATH { "rsfq_databases/NDFF_2IN.csv" } ;
const std::string NDFF_3_INPUT_PATH { "rsfq_databases/NDFF_3IN.csv" } ; 

int main()
{
  using namespace mockturtle;

  // Import benchmarks
  auto benchmarks1 = experiments::epfl_benchmarks( experiments::adder | experiments::sin | experiments::cavlc | experiments::int2float | experiments::priority | experiments::i2c | experiments::voter | experiments::dec );
  auto benchmarks2 = experiments::iscas_benchmarks( experiments::c432 | experiments::c499 | experiments::c880 | experiments::c1355 | experiments::c1908 | experiments::c3540 | experiments::c5315 | experiments::c7552 );
  benchmarks1.insert(benchmarks1.end(), benchmarks2.begin(), benchmarks2.end());
  
  experiments::experiment<std::string, 
            double, double, double, 
            double, double, double, 
            double, double, double, 
            float, float, float> exp(
      "mapper", "benchmark", 
      "#DFF(2)", "#DFF(3)", "#DFF(ratio)", 
      "area(2)", "area(3)", "area(ratio)", 
      "delay(2)", "delay(3)", "delay(ratio)", 
      "time(2)",  "time(3)", "time(ratio)");

  fmt::print( "[i] processing technology library\n" );

  /* library to map to technology */
  std::ifstream inputFile_2_IN( DATABASE_2_INPUT_PATH );
  std::ifstream inputFile_3_IN( DATABASE_3_INPUT_PATH );
  
  std::vector<gate> gates_2_IN;
  std::vector<gate> gates_3_IN;

  std::unordered_map<std::string, int> NDFF_2_IN = readCSV( NDFF_2_INPUT_PATH );
  std::unordered_map<std::string, int> NDFF_3_IN = readCSV( NDFF_3_INPUT_PATH );
  
  if ( lorina::read_genlib( inputFile_2_IN, genlib_reader( gates_2_IN ) ) != lorina::return_code::success ) { return 1; }
  if ( lorina::read_genlib( inputFile_3_IN, genlib_reader( gates_3_IN ) ) != lorina::return_code::success ) { return 1; }

  mockturtle::tech_library_params tps; // tps.verbose = true;
  tps.load_minimum_size_only = false;
  tps.remove_dominated_gates = false;
  tech_library<4, mockturtle::classification_type::p_configurations> tech_lib_2_in( gates_2_IN, tps );
  tech_library<4, mockturtle::classification_type::p_configurations> tech_lib_3_in( gates_3_IN, tps );


  for ( auto const& benchmark : benchmarks1 )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    mockturtle::aig_network tech_indep_ntk;
    if ( lorina::read_aiger( experiments::benchmark_path( benchmark ), aiger_reader( tech_indep_ntk ) ) != lorina::return_code::success ) { continue; }
  
    auto [res_2_in, st_2_in, PB_DFF_2_in, nJJ_2_in, cec_2_in] = map_with_pb(benchmark, tech_indep_ntk, tech_lib_2_in, NDFF_2_IN, false);
    auto [res_3_in, st_3_in, PB_DFF_3_in, nJJ_3_in, cec_3_in] = map_with_pb(benchmark, tech_indep_ntk, tech_lib_3_in, NDFF_3_IN, false);

    exp( benchmark,   PB_DFF_2_in,   PB_DFF_3_in, static_cast<double>(PB_DFF_3_in)  /PB_DFF_2_in  , 
                         nJJ_2_in,      nJJ_3_in, static_cast<double>(nJJ_3_in)     /nJJ_2_in     , 
                    st_2_in.delay, st_3_in.delay, static_cast<double>(st_3_in.delay)/st_2_in.delay, 
                    to_seconds( st_2_in.time_total * 1000 ), to_seconds( st_3_in.time_total * 1000 ), to_seconds( st_3_in.time_total )/to_seconds( st_2_in.time_total ) ); 
    
    exp.save();
    exp.table(); // print continuously
  }

  return 0;
}

/*

  
  // PBMap baseline
  // std::unordered_map<std::string, std::tuple<double,double,double>> PBMAP;
  // PBMAP["int2float"] = std::make_tuple(  270,   6432,  16);
  // PBMAP["priority"]  = std::make_tuple( 9064, 102085, 127);
  // PBMAP["sin"]       = std::make_tuple(13666, 215318, 182);
  // PBMAP["cavlc"]     = std::make_tuple(  522,  16339,  17);
  // PBMAP["dec"]       = std::make_tuple(    8,   5469,   4);
  // PBMAP["c499"]      = std::make_tuple(  476,   7758,  13);
  // PBMAP["c880"]      = std::make_tuple(  774,  12909,  22);
  // PBMAP["c1908"]     = std::make_tuple(  696,  12013,  20);
  // PBMAP["c3540"]     = std::make_tuple( 1159,  28300,  31);
  // PBMAP["c5315"]     = std::make_tuple( 2908,  52033,  23);
  // PBMAP["c7552"]     = std::make_tuple( 2429,  48482,  19);

with XMG:
| benchmark |  #DFF(2) |  #DFF(3) | #DFF(ratio) |   area(2) |   area(3) | area(ratio) | delay(2) | delay(3) | delay(ratio) | time(2) |  time(3) | time(ratio) |
|     adder | 48768.00 | 32388.00 |        0.66 | 156369.00 | 109561.00 |        0.70 |   128.00 |    86.00 |         0.67 |  441.37 |   885.58 |        2.01 |
|       sin | 17627.00 | 11988.00 |        0.68 | 136606.00 | 105377.00 |        0.77 |    92.00 |    69.00 |         0.75 | 4941.55 | 13128.23 |        2.66 |
|     cavlc |   987.00 |   433.00 |        0.44 |  15624.00 |  13311.00 |        0.85 |    15.00 |     7.00 |         0.47 |  148.93 |   494.74 |        3.32 |
|       dec |    16.00 |     2.00 |        0.12 |   6324.00 |   8848.00 |        1.40 |     4.00 |     2.00 |         0.50 |   86.58 |   318.61 |        3.68 |
|       i2c |  4830.00 |  2184.00 |        0.45 |  38335.00 |  27148.00 |        0.71 |    16.00 |     7.00 |         0.44 |  252.00 |   687.65 |        2.73 |
| int2float |   443.00 |   182.00 |        0.41 |   5973.00 |   4521.00 |        0.76 |    15.00 |     6.00 |         0.40 |   59.88 |   188.39 |        3.15 |
|  priority | 14754.00 | 10084.00 |        0.68 |  68177.00 |  44016.00 |        0.65 |    84.00 |    62.00 |         0.74 |  283.31 |   996.63 |        3.52 |
|     voter |  8446.00 |  3348.00 |        0.40 | 189622.00 | 185785.00 |        0.98 |    41.00 |    25.00 |         0.61 | 8846.34 | 27816.97 |        3.14 |
|      c432 |  1180.00 |   535.00 |        0.45 |   6905.00 |   4948.00 |        0.72 |    21.00 |    11.00 |         0.52 |   61.66 |   271.14 |        4.40 |
|      c499 |   512.00 |   430.00 |        0.84 |   5607.00 |   4847.00 |        0.86 |     8.00 |     7.00 |         0.88 |  354.31 |   886.68 |        2.50 |
|      c880 |  1179.00 |   780.00 |        0.66 |   8650.00 |   7283.00 |        0.84 |    15.00 |    10.00 |         0.67 |  123.91 |   356.36 |        2.88 |
|     c1355 |   448.00 |   398.00 |        0.89 |   5703.00 |   5055.00 |        0.89 |     8.00 |     7.00 |         0.88 |  474.98 |  1281.17 |        2.70 |
|     c1908 |   799.00 |   435.00 |        0.54 |   5497.00 |   5667.00 |        1.03 |    13.00 |     8.00 |         0.62 |  239.18 |   619.66 |        2.59 |
|     c3540 |  1556.00 |   936.00 |        0.60 |  20820.00 |  19131.00 |        0.92 |    23.00 |    13.00 |         0.57 |  458.06 |  1358.27 |        2.97 |
|     c5315 |  3727.00 |  2798.00 |        0.75 |  33263.00 |  31139.00 |        0.94 |    15.00 |    11.00 |         0.73 | 1039.38 |  2798.78 |        2.69 |
|     c7552 |  4744.00 |  2201.00 |        0.46 |  31792.00 |  25084.00 |        0.79 |    19.00 |     9.00 |         0.47 | 1211.71 |  3125.85 |        2.58 |

with XAG:
| benchmark |  #DFF(2) |  #DFF(3) | #DFF(ratio) |   area(2) |   area(3) | area(ratio) | delay(2) | delay(3) | delay(ratio) | time(2) |  time(3) | time(ratio) |
|     adder | 48768.00 | 32548.00 |        0.67 | 156369.00 | 110241.00 |        0.71 |   128.00 |    86.00 |         0.67 |  403.67 |   819.39 |        2.03 |
|       sin | 18796.00 | 12271.00 |        0.65 | 140548.00 | 107157.00 |        0.76 |    92.00 |    69.00 |         0.75 | 4352.81 | 12406.18 |        2.85 |
|     cavlc |   999.00 |   435.00 |        0.44 |  15656.00 |  13326.00 |        0.85 |    15.00 |     7.00 |         0.47 |  135.48 |   463.75 |        3.42 |
|       dec |    16.00 |     2.00 |        0.12 |   6324.00 |   8848.00 |        1.40 |     4.00 |     2.00 |         0.50 |   84.17 |   291.84 |        3.47 |
|       i2c |  4791.00 |  2218.00 |        0.46 |  38205.00 |  27206.00 |        0.71 |    16.00 |     7.00 |         0.44 |  226.60 |   605.00 |        2.67 |
| int2float |   433.00 |   182.00 |        0.42 |   5949.00 |   4521.00 |        0.76 |    15.00 |     6.00 |         0.40 |   53.70 |   163.92 |        3.05 |
|  priority | 14724.00 | 10205.00 |        0.69 |  68123.00 |  44337.00 |        0.65 |    84.00 |    62.00 |         0.74 |  247.44 |   969.58 |        3.92 |
|     voter |  7819.00 |  3633.00 |        0.46 | 186541.00 | 186363.00 |        1.00 |    41.00 |    25.00 |         0.61 | 7393.02 | 25590.47 |        3.46 |
|      c432 |  1192.00 |   535.00 |        0.45 |   6965.00 |   4948.00 |        0.71 |    21.00 |    11.00 |         0.52 |   53.66 |   251.43 |        4.69 |
|      c499 |   518.00 |   430.00 |        0.83 |   5601.00 |   4847.00 |        0.87 |     8.00 |     7.00 |         0.88 |  317.20 |   793.30 |        2.50 |
|      c880 |  1183.00 |   780.00 |        0.66 |   8641.00 |   7283.00 |        0.84 |    15.00 |    10.00 |         0.67 |  113.18 |   323.12 |        2.85 |
|     c1355 |   454.00 |   398.00 |        0.88 |   5697.00 |   5055.00 |        0.89 |     8.00 |     7.00 |         0.88 |  424.24 |  1180.26 |        2.78 |
|     c1908 |   804.00 |   450.00 |        0.56 |   5547.00 |   5655.00 |        1.02 |    13.00 |     8.00 |         0.62 |  207.90 |   582.31 |        2.80 |
|     c3540 |  1562.00 |   923.00 |        0.59 |  20883.00 |  19113.00 |        0.92 |    23.00 |    13.00 |         0.57 |  400.10 |  1271.54 |        3.18 |
|     c5315 |  3807.00 |  2805.00 |        0.74 |  33282.00 |  31145.00 |        0.94 |    15.00 |    11.00 |         0.73 |  939.55 |  2562.27 |        2.73 |
|     c7552 |  4566.00 |  2188.00 |        0.48 |  31096.00 |  24864.00 |        0.80 |    19.00 |     9.00 |         0.47 | 1055.31 |  2963.02 |        2.81 |

*/