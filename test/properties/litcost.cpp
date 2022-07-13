#include <catch.hpp>

#include <kitty/cube.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <mockturtle/properties/litcost.hpp>

using namespace mockturtle;

TEST_CASE( "count literals constant'", "[litcost]" )
{
  std::vector<kitty::cube> sop = { {0, 0} };

  CHECK( literal_cost( sop, 2 ) == 0 );
}

TEST_CASE( "count literals ab'", "[litcost]" )
{
  std::vector<kitty::cube> sop = { {2,3} };

  CHECK( literal_cost( sop, 2 ) == 2 );
}

TEST_CASE( "count literals ab + ac", "[litcost]" )
{
  std::vector<kitty::cube> sop = { {6, 6}, {5, 5} };

  CHECK( literal_cost( sop, 3 ) == 3 );
}

TEST_CASE( "count literals tricky 4-input", "[litcost]" )
{
  std::vector<kitty::cube> sop = { {8, 12}, {8, 11}, {11, 11}, {5, 15}, {9, 15} };

  CHECK( literal_cost( sop, 4 ) == 12 );
}

TEST_CASE( "count literals from truth table", "[litcost]" )
{
  kitty::dynamic_truth_table tt( 3 );
  tt._bits[0] = 0xe8;

  CHECK( literal_cost( tt ) == 5 );
}

TEST_CASE( "count literals from truth table with don't cares", "[litcost]" )
{
  kitty::dynamic_truth_table tt( 3 );
  kitty::dynamic_truth_table dc( 3 );
  tt._bits[0] = 0xe8;
  dc._bits[0] = 0x88;

  CHECK( literal_cost( tt, dc ) == 3 );
}
