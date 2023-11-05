#include <catch.hpp>

#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/algorithms/ac_decomposition.hpp>

using namespace mockturtle;

TEST_CASE( "ACD function 6 vars FS 2", "[ac_decomposition]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  ac_decomposition_params ps;
  ps.lut_size = 4;
  detail::ac_decomposition_impl acd( tt, 6, ps );

  CHECK( acd.run_no_permutations( 2 ) == 4 );
  CHECK( acd.run( 2 ) == 4 );
  CHECK( acd.run_offset( 2, 1 ) == 4 );
}

TEST_CASE( "ACD function 6 vars FS 1", "[ac_decomposition]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  ac_decomposition_params ps;
  ps.lut_size = 4;
  detail::ac_decomposition_impl acd( tt, 6, ps );

  CHECK( acd.run_no_permutations( 1 ) == 3 );
  CHECK( acd.run( 1 ) == 3 );
}

TEST_CASE( "ACD function 6 vars FS 3", "[ac_decomposition]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;
  
  ac_decomposition_params ps;
  ps.lut_size = 4;
  detail::ac_decomposition_impl acd( tt, 6, ps );

  CHECK( acd.run_no_permutations( 3 ) == 8 );
  CHECK( acd.run( 3 ) == 5 );
  CHECK( acd.run_offset( 3, 2 ) == 5 );
}

TEST_CASE( "ACD function 6 vars multiple FS", "[ac_decomposition]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  ac_decomposition_params ps;
  ps.lut_size = 4;
  detail::ac_decomposition_impl acd( tt, 6, ps );

  CHECK( acd.run() == 4 );
}

TEST_CASE( "ACD function 6 vars late arriving", "[ac_decomposition]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  ac_decomposition_params ps;
  ps.lut_size = 4;

  std::vector<uint32_t> late_arriving = { 1 };

  {
    detail::ac_decomposition_impl acd( tt, 6, ps );
    CHECK( acd.run( late_arriving ) == 4 );
  }

  {
    detail::ac_decomposition_impl acd( tt, 6, ps );

    late_arriving.push_back( 2 );
    CHECK( acd.run( late_arriving ) == UINT32_MAX );
  }

  {
    detail::ac_decomposition_impl acd( tt, 6, ps );

    late_arriving.pop_back();
    late_arriving.push_back( 0 );
    CHECK( acd.run( late_arriving ) == 4 );
  }
}

TEST_CASE( "ACD function 8 vars FS 2", "[ac_decomposition]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

  ac_decomposition_params ps;
  ps.lut_size = 6;
  detail::ac_decomposition_impl acd( tt, 8, ps );

  CHECK( acd.run_no_permutations( 2 ) == 9 );
  CHECK( acd.run( 2 ) == 4 );
  CHECK( acd.run_offset( 2, 1 ) == 7 );
}

TEST_CASE( "ACD function 8 vars multiple FS", "[ac_decomposition]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

  ac_decomposition_params ps;
  ps.lut_size = 6;
  detail::ac_decomposition_impl acd( tt, 8, ps );

  CHECK( acd.run() == 4 );
}

TEST_CASE( "ACD function 8 vars late arriving", "[ac_decomposition]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

  ac_decomposition_params ps;
  ps.lut_size = 6;

  std::vector<uint32_t> late_arriving = { 2 };

  {
    detail::ac_decomposition_impl acd( tt, 8, ps );

    CHECK( acd.run( late_arriving ) == 4 );
  }

  {
    detail::ac_decomposition_impl acd( tt, 8, ps );

    late_arriving.push_back( 3 );
    CHECK( acd.run( late_arriving ) == 4 );
  }

  {
    detail::ac_decomposition_impl acd( tt, 8, ps );

    late_arriving.push_back( 6 );
    CHECK( acd.run( late_arriving ) == 7 );
  }

  {
    detail::ac_decomposition_impl acd( tt, 8, ps );

    late_arriving.pop_back();
    late_arriving.pop_back();
    late_arriving.push_back( 7 );
    CHECK( acd.run( late_arriving ) == 6 );
  }
}
