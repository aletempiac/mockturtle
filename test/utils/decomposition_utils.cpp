#include <catch.hpp>

#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/static_truth_table.hpp>
#include <mockturtle/utils/decomposition_utils.hpp>

using namespace mockturtle;

TEST_CASE( "Function 6 vars FS 2", "[decomposition_utils]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  CHECK( acd_column_multiplicity( tt, 2 ) == 4 );
  CHECK( acd_enumerate_combinations( tt, 2 ) == 4 );
}

TEST_CASE( "Function 6 vars FS 1", "[decomposition_utils]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  CHECK( acd_column_multiplicity( tt, 1 ) == 3 );

  CHECK( acd_enumerate_combinations( tt, 1 ) == 3 );
}

TEST_CASE( "Function 6 vars FS 3", "[decomposition_utils]" )
{
  kitty::static_truth_table<6> tt;
  tt._bits = 0x8804800184148111;

  CHECK( acd_column_multiplicity( tt, 3 ) == 8 );
  CHECK( acd_enumerate_combinations( tt, 3 ) == 5 );
}

TEST_CASE( "Function 8 vars FS 2", "[decomposition_utils]" )
{
  kitty::static_truth_table<8> tt;
  kitty::create_from_hex_string( tt, "000000001000200000000000000000020000100100001000C009800BC00D800F" );

  CHECK( acd_column_multiplicity( tt, 2 ) == 9 );
  CHECK( acd_enumerate_combinations( tt, 2 ) == 5 );
}
