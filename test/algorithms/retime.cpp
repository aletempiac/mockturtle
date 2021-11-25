#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/constructors.hpp>

#include <mockturtle/networks/generic.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/algorithms/retime.hpp>

using namespace mockturtle;

TEST_CASE( "Retime forward 1", "[retime]" )
{
  generic_network ntk;
  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto b1 = ntk.create_latch( a );
  const auto b2 = ntk.create_latch( b );
  const auto f = ntk.create_and( b1, b2 );

  ntk.create_po( f );

  retime( ntk );
  // CHECK( luts.size() == 6u );
}

TEST_CASE( "Retime backward 1", "[retime]" )
{
  generic_network ntk;
  const auto a = ntk.create_pi();
  const auto x1 = ntk.create_not( a );
  const auto x2 = ntk.create_buf( a );

  const auto b1 = ntk.create_latch( x1 );
  const auto b2 = ntk.create_latch( x2 );

  ntk.create_po( b1 );
  ntk.create_po( b2 );

  retime( ntk );
  // CHECK( luts.size() == 6u );
}

TEST_CASE( "Zero retime forward", "[retime]" )
{
  generic_network ntk;
  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto b1 = ntk.create_latch( a );
  const auto b2 = ntk.create_latch( b );

  const auto x1 = ntk.create_and( b1, b2 );
  const auto x2 = ntk.create_or( b1, b2 );

  ntk.create_po( x1 );
  ntk.create_po( x2 );

  retime( ntk );
  // CHECK( luts.size() == 6u );
}

TEST_CASE( "Retime forward 2", "[retime]" )
{
  generic_network ntk;
  const auto a = ntk.create_pi();
  const auto b = ntk.create_pi();
  const auto c = ntk.create_pi();
  const auto b1 = ntk.create_latch( a );
  const auto b2 = ntk.create_latch( b );
  const auto b3 = ntk.create_latch( c );

  const auto x1 = ntk.create_and( b1, b2 );
  const auto x2 = ntk.create_or( x1, b3 );

  const auto b4 = ntk.create_latch( x2 );

  ntk.create_po( b4 );

  retime( ntk );
  // CHECK( luts.size() == 6u );
}