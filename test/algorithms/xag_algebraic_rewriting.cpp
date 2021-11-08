#include <catch.hpp>

#include <mockturtle/traits.hpp>
#include <mockturtle/algorithms/xag_algebraic_rewriting.hpp>
#include <mockturtle/networks/xag.hpp>
#include <mockturtle/views/depth_view.hpp>

using namespace mockturtle;

TEST_CASE( "AIG depth optimization with associativity", "[xag_algebraic_rewriting]" )
{
  xag_network xag;

  const auto a = xag.create_pi();
  const auto b = xag.create_pi();
  const auto c = xag.create_pi();
  const auto d = xag.create_pi();

  const auto f1 = xag.create_and( a, b );
  const auto f2 = xag.create_and( f1, c );
  const auto f3 = xag.create_and( f2, d );

  xag.create_po( f3 );

  depth_view depth_xag{xag};

  CHECK( depth_xag.depth() == 3 );

  xag_algebraic_depth_rewriting( depth_xag );

  CHECK( depth_xag.depth() == 2 );
}

TEST_CASE( "XOR depth optimization with associativity", "[xag_algebraic_rewriting]" )
{
  xag_network xag;

  const auto a = xag.create_pi();
  const auto b = xag.create_pi();
  const auto c = xag.create_pi();
  const auto d = xag.create_pi();

  const auto f1 = xag.create_xor( a, b );
  const auto f2 = xag.create_xor( f1, c );
  const auto f3 = xag.create_xor( f2, d );

  xag.create_po( f3 );

  depth_view depth_xag{xag};

  CHECK( depth_xag.depth() == 3 );

  xag_algebraic_depth_rewriting( depth_xag );

  CHECK( depth_xag.depth() == 2 );
}

TEST_CASE( "AND-XOR depth optimization with distributivity", "[xag_algebraic_rewriting]" )
{
  xag_network xag;

  const auto a = xag.create_pi();
  const auto b = xag.create_pi();
  const auto c = xag.create_pi();
  const auto d = xag.create_pi();
  const auto e = xag.create_pi();

  const auto f1 = xag.create_nand( a, b );
  const auto f2 = xag.create_and( f1, c );
  const auto f3 = xag.create_xor( f2, d );
  const auto f4 = xag.create_and( f3, e );

  xag.create_po( f4 );

  depth_view depth_xag{xag};

  CHECK( depth_xag.depth() == 4 );

  xag_algebraic_depth_rewriting( depth_xag );

  CHECK( depth_xag.depth() == 3 );
}

TEST_CASE( "AND-OR depth optimization with distributivity negation", "[xag_algebraic_rewriting]" )
{
  xag_network xag;

  const auto a = xag.create_pi();
  const auto b = xag.create_pi();
  const auto c = xag.create_pi();
  const auto d = xag.create_pi();
  const auto e = xag.create_pi();

  const auto f1 = xag.create_nand( a, b );
  const auto f2 = xag.create_and( f1, c );
  const auto f3 = xag.create_or( f2, d );
  const auto f4 = xag.create_and( f3, e );

  xag.create_po( f4 );

  depth_view depth_xag{xag};

  CHECK( depth_xag.depth() == 4 );

  xag_algebraic_depth_rewriting( depth_xag );

  CHECK( depth_xag.depth() == 3 );
}

TEST_CASE( "AND-XOR depth optimization with distributivity negation", "[xag_algebraic_rewriting]" )
{
  xag_network xag;

  const auto a = xag.create_pi();
  const auto b = xag.create_pi();
  const auto c = xag.create_pi();
  const auto d = xag.create_pi();
  const auto e = xag.create_pi();

  const auto f1 = xag.create_nand( a, b );
  const auto f2 = xag.create_nand( f1, c );
  const auto f3 = xag.create_xor( f2, d );
  const auto f4 = xag.create_and( f3, e );

  xag.create_po( f4 );

  depth_view depth_xag{xag};

  CHECK( depth_xag.depth() == 4 );

  xag_algebraic_depth_rewriting( depth_xag );

  CHECK( depth_xag.depth() == 3 );
}