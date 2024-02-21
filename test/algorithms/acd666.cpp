#include <catch.hpp>

#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/algorithms/acd666.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/klut.hpp>

using namespace mockturtle;

TEST_CASE( "ACD666 failing 2", "[acd666]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "fffffffffffffffffffff88888888fff00000ffffffff0000000088888888000" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd666_impl acd( 8 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD666 8 vars failing", "[acd666]" )
{
  uint64_t tt_c[4] = { 0xf1e3c78f1f3e7cf8, 0xe1c3870f1e3c78f0, 0xf0e0c3830f0e3c38, 0xe0c083030e0c3830 };

  acd666_impl acd( 8 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == false );
  // CHECK( acd.compute_decomposition() == -1 );
}

TEST_CASE( "ACD666 7 vars 1 BS and 0 SS", "[acd666]" )
{
  uint64_t tt_c[2] = { 0x7f807f807f80807f, 0x7f807f807f807f80 };

  acd666_impl acd( 7 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  // CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD666 8 vars 1 BS and 0 SS", "[acd666]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "1000000000000000000000000000000000000000000000000000000000000000" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd666_impl acd( 8 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  // CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD666 8 vars 1 BS and 1 SS", "[acd666]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "FF30CE44FDB8FDB8FF30CE44FDB8FDB8FF30FF30CCCCCCCCFF30CE44CCCCCCCC" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd666_impl acd( 8 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  // CHECK( acd.compute_decomposition() == 0 );
}

TEST_CASE( "ACD666 8 vars 1 BS and 1 SS and 3 U", "[acd666]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "100000000000000000000000000000000000000000000000000000000000FFFF" );

  uint64_t tt_c[4] = { tt._bits[0], tt._bits[1], tt._bits[2], tt._bits[3] };

  acd666_impl acd( 8 );

  // unsigned char decompArray[92];

  CHECK( acd.run( tt_c ) == true );
  // CHECK( acd.compute_decomposition() == 0 );
}