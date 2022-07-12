#include <catch.hpp>
#include <unordered_set>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>
#include <kitty/static_truth_table.hpp>
#include <kitty/npn.hpp>

#include <mockturtle/algorithms/node_resynthesis/sop_factoring.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/node_map.hpp>

using namespace mockturtle;

TEST_CASE( "SOP factoring for 4-NPN functions", "[sop_factoring]" )
{
  using tt_hash = kitty::hash<kitty::dynamic_truth_table>;

  std::unordered_set<kitty::dynamic_truth_table, tt_hash> classes;
  kitty::dynamic_truth_table tt( 4 );
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();

  std::vector<aig_network::signal> pis = {a, b, c, d};

  sop_factoring<aig_network> resyn;

  default_simulator<kitty::dynamic_truth_table> sim( aig.num_pis() );

  for ( auto const& t : classes )
  {
    resyn( aig, t, pis.begin(), pis.end(), [&]( auto const& f ) {
      aig.create_po( f );
    } );
  }

  auto tts = simulate<kitty::dynamic_truth_table, aig_network>( aig, sim );

  auto it = classes.begin();
  for ( auto i = 0u; i < classes.size(); ++i, ++it )
  {
    CHECK( *it == tts[i] );
  }
}

TEST_CASE( "SOP factoring for 4-NPN functions with don't cares", "[sop_factoring]" )
{
  using tt_hash = kitty::hash<kitty::dynamic_truth_table>;

  std::unordered_set<kitty::dynamic_truth_table, tt_hash> classes;
  kitty::dynamic_truth_table tt( 4 );
  do
  {
    const auto res = kitty::exact_npn_canonization( tt );
    classes.insert( std::get<0>( res ) );
    kitty::next_inplace( tt );
  } while ( !kitty::is_const0( tt ) );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();

  std::vector<aig_network::signal> pis = {a, b, c, d};

  sop_factoring<aig_network> resyn;

  default_simulator<kitty::dynamic_truth_table> sim( aig.num_pis() );

  kitty::dynamic_truth_table dc( 4 );
  dc._bits[0] = 0x6;

  for ( auto const& t : classes )
  {
    resyn( aig, t, dc, pis.begin(), pis.end(), [&]( auto const& f ) {
      aig.create_po( f );
    } );
  }

  auto tts = simulate<kitty::dynamic_truth_table, aig_network>( aig, sim );

  auto it = classes.begin();
  for ( auto i = 0u; i < classes.size(); ++i, ++it )
  {
    CHECK( ( *it & ~dc ) == ( tts[i] & ~dc ) );
  }
}
