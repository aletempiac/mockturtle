#include <catch.hpp>

#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/algorithms/acd66.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/klut.hpp>

using namespace mockturtle;

// TEST_CASE( "ACD function 5 vars FS 2", "[acd66]" )
// {
//   kitty::static_truth_table<5> tt;
//   tt._bits = 0x01221002;

//   acd66_params ps;
//   ps.lut_size = 4;
//   detail::acd66_impl acd( tt, 5, ps );

//   std::vector<uint32_t> late_arriving = { 0, 1 };

//   CHECK( acd.run( late_arriving ) == 3 );

//   auto res = acd.get_result_ntk();
//   CHECK( res.has_value() );

//   klut_network klut = *res;

//   const auto tt_res = simulate<kitty::static_truth_table<5>>( klut )[0];

//   CHECK( tt_res == tt );
// }

// TEST_CASE( "ACD function 6 vars FS 2", "[acd66]" )
// {
//   kitty::static_truth_table<6> tt;
//   tt._bits = 0x8804800184148111;

//   acd66_params ps;
//   ps.lut_size = 4;
//   detail::acd66_impl acd( tt, 6, ps );

//   CHECK( acd.run( 2 ) == 4 );
// }

// TEST_CASE( "ACD function 6 vars FS 1", "[acd66]" )
// {
//   kitty::static_truth_table<6> tt;
//   tt._bits = 0x8804800184148111;

//   acd66_params ps;
//   ps.lut_size = 4;
//   detail::acd66_impl acd( tt, 6, ps );

//   CHECK( acd.run( 1 ) == UINT32_MAX );
// }

// TEST_CASE( "ACD function 6 vars FS 3", "[acd66]" )
// {
//   kitty::static_truth_table<6> tt;
//   tt._bits = 0x8804800184148111;
  
//   acd66_params ps;
//   ps.lut_size = 4;
//   detail::acd66_impl acd( tt, 6, ps );

//   CHECK( acd.run( 3 ) == UINT32_MAX );
// }

// TEST_CASE( "ACD function 6 vars multiple FS", "[acd66]" )
// {
//   kitty::static_truth_table<6> tt;
//   tt._bits = 0x8804800184148111;

//   acd66_params ps;
//   ps.lut_size = 4;
//   detail::acd66_impl acd( tt, 6, ps );

//   CHECK( acd.run() == 4 );
// }

TEST_CASE( "ACD66 8 vars failing", "[acd66]" )
{
  uint64_t tt_c[4] = { 0xf1e3c78f1f3e7cf8, 0xe1c3870f1e3c78f0, 0xf0e0c3830f0e3c38, 0xe0c083030e0c3830 };

  acd66_params ps;
  acd66_stats st;
  acd66_impl acd( 8, ps, &st );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == false );
  CHECK( acd.compute_decomposition() == -1 );
}

TEST_CASE( "ACD66 7 vars 1 BS and 0 SS", "[acd66]" )
{
  uint64_t tt_c[2] = { 0x7f807f807f80807f, 0x7f807f807f807f80 };

  acd66_params ps;
  acd66_stats st;
  acd66_impl acd( 7, ps, &st );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD66 8 vars 1 BS and 0 SS", "[acd66]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "1000000000000000000000000000000000000000000000000000000000000000" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd66_params ps;
  acd66_stats st;
  acd66_impl acd( 8, ps, &st );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD66 8 vars 1 BS and 1 SS", "[acd66]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "FF30CE44FDB8FDB8FF30CE44FDB8FDB8FF30FF30CCCCCCCCFF30CE44CCCCCCCC" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd66_params ps;
  acd66_stats st;
  acd66_impl acd( 8, ps, &st );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD66 8 vars 1 BS and 1 SS and 3 U", "[acd66]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "100000000000000000000000000000000000000000000000000000000000FFFF" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd66_params ps;
  acd66_stats st;
  acd66_impl acd( 8, ps, &st );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  CHECK( acd.compute_decomposition() == 0 );
}

// TEST_CASE( "ACD function 6 vars late arriving", "[acd66]" )
// {
//   kitty::static_truth_table<6> tt;
//   tt._bits = 0x8804800184148111;

//   acd66_params ps;
//   ps.lut_size = 4;

//   std::vector<uint32_t> late_arriving = { 1 };

//   {
//     detail::acd66_impl acd( tt, 6, ps );
//     CHECK( acd.run( late_arriving ) == 4 );

//     auto res = acd.get_result_ntk();
//     CHECK( res.has_value() );

//     klut_network klut = *res;

//     const auto tt_res = simulate<kitty::static_truth_table<6>>( klut )[0];

//     CHECK( tt_res == tt );
//   }

//   {
//     detail::acd66_impl acd( tt, 6, ps );

//     late_arriving.push_back( 2 );
//     CHECK( acd.run( late_arriving ) == UINT32_MAX );
//   }

//   {
//     detail::acd66_impl acd( tt, 6, ps );

//     late_arriving.pop_back();
//     late_arriving.push_back( 0 );
//     CHECK( acd.run( late_arriving ) == 4 );
//   }
// }

// TEST_CASE( "ACD function 8 vars FS 2", "[acd66]" )
// {
//   kitty::static_truth_table<8> tt;
//   kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

//   acd66_params ps;
//   ps.lut_size = 6;
//   detail::acd66_impl acd( tt, 8, ps );

//   CHECK( acd.run( 2 ) == 4 );
// }

// TEST_CASE( "ACD function 8 vars multiple FS", "[acd66]" )
// {
//   kitty::static_truth_table<8> tt;
//   kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

//   acd66_params ps;
//   ps.lut_size = 6;
//   detail::acd66_impl acd( tt, 8, ps );

//   CHECK( acd.run() == 4 );
// }

// TEST_CASE( "ACD function 8 vars late arriving", "[acd66]" )
// {
//   kitty::static_truth_table<8> tt;
//   kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

//   acd66_params ps;
//   ps.lut_size = 6;

//   std::vector<uint32_t> late_arriving = { 2 };

//   {
//     detail::acd66_impl acd( tt, 8, ps );
//     CHECK( acd.run( late_arriving ) == 4 );

//     auto res = acd.get_result_ntk();
//     CHECK( res.has_value() );

//     klut_network klut = *res;

//     const auto tt_res = simulate<kitty::static_truth_table<8>>( klut )[0];

//     CHECK( tt_res == tt );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.push_back( 3 );
//     CHECK( acd.run( late_arriving ) == 4 );

//     auto res = acd.get_result_ntk();
//     CHECK( res.has_value() );

//     klut_network klut = *res;

//     const auto tt_res = simulate<kitty::static_truth_table<8>>( klut )[0];

//     CHECK( tt_res == tt );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.push_back( 6 );
//     CHECK( acd.run( late_arriving ) == 7 );

//     auto res = acd.get_result_ntk();
//     CHECK( res.has_value() );

//     klut_network klut = *res;

//     const auto tt_res = simulate<kitty::static_truth_table<8>>( klut )[0];

//     CHECK( tt_res == tt );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.pop_back();
//     late_arriving.pop_back();
//     late_arriving.push_back( 7 );
//     CHECK( acd.run( late_arriving ) == 5 );

//     auto res = acd.get_result_ntk();
//     CHECK( res.has_value() );

//     klut_network klut = *res;

//     const auto tt_res = simulate<kitty::static_truth_table<8>>( klut )[0];

//     CHECK( tt_res == tt );
//   }
// }

// TEST_CASE( "ACD function 8 vars DSD late arriving", "[acd66]" )
// {
//   kitty::static_truth_table<8> tt;
//   kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

//   acd66_params ps;
//   ps.lut_size = 6;

//   std::vector<uint32_t> late_arriving = { 2 };

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     CHECK( acd.run_dsd( late_arriving ) == 4 );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.push_back( 3 );
//     CHECK( acd.run_dsd( late_arriving ) == 4 );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.push_back( 6 );
//     CHECK( acd.run_dsd( late_arriving ) == 7 );
//   }

//   {
//     detail::acd66_impl acd( tt, 8, ps );

//     late_arriving.pop_back();
//     late_arriving.pop_back();
//     late_arriving.push_back( 7 );
//     CHECK( acd.run_dsd( late_arriving ) == 6 );
//   }
// }

// TEST_CASE( "ACD function 7 vars with don't cares", "[acd66]" )
// {
//   kitty::static_truth_table<7> tt;
//   kitty::create_from_hex_string( tt, "02020202020200020202020202020202" );

//   acd66_params ps;
//   ps.lut_size = 6;

//   std::vector<uint32_t> late_arriving = { 3, 6 };

//   detail::acd66_impl acd( tt, 7, ps );
//   CHECK( acd.run( late_arriving ) == 3 );
//   CHECK( acd.verify_equivalence() );
// }

// TEST_CASE( "ACD for shannon", "[acd66]" )
// {
//   kitty::static_truth_table<5> tt;
//   kitty::create_from_hex_string( tt, "1241fa11" );

//   acd66_params ps;
//   ps.lut_size = 4;

//   std::vector<uint32_t> late_arriving = { 1 };

//   detail::acd66_impl acd( tt, 5, ps );
//   CHECK( acd.run( late_arriving ) == 4 );
//   CHECK( acd.verify_equivalence() );
// }

// TEST_CASE( "ACD on hard function", "[acd66]" )
// {
//   kitty::static_truth_table<8> tt;
//   kitty::create_from_hex_string( tt, "71f3ffff71f3333371f3777771f3111130717777307111113071333330710000" );

//   acd66_params ps;
//   ps.lut_size = 6;

//   detail::acd66_impl acd( tt, 8, ps );
//   std::vector<uint32_t> late_arriving = { 4, 7 };
//   CHECK( acd.run( late_arriving ) == 5 );
// }
