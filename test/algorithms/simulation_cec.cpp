#include <catch.hpp>

#include <mockturtle/algorithms/mapper.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/simulation_cec.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/xag.hpp>

#include <fmt/format.h>
#include <lorina/aiger.hpp>

using namespace mockturtle;

TEST_CASE( "CEC simple AIG", "[cec]" )
{
  aig_network aig1;

  const auto a1 = aig1.create_pi();
  const auto b1 = aig1.create_pi();
  const auto c1 = aig1.create_pi();
  const auto f11 = aig1.create_and( a1, b1 );
  const auto f21 = aig1.create_and( c1, f11 );
  aig1.create_po( f21 );

  aig_network aig2;

  const auto a2 = aig2.create_pi();
  const auto b2 = aig2.create_pi();
  const auto c2 = aig2.create_pi();
  const auto f12 = aig2.create_and( b2, c2 );
  const auto f22 = aig2.create_and( a2, f12 );
  aig2.create_po( f22 );

  simulation_cec_stats st;
  CHECK( *simulation_cec( aig1, aig2, &st ) );
  CHECK( st.split_var == 3 );
  CHECK( st.rounds == 1 );
}


TEST_CASE( "CEC different #PIs", "[cec]" )
{
  xag_network xag1;

  const auto a1 = xag1.create_pi();
  const auto b1 = xag1.create_pi();
  const auto c1 = xag1.create_pi();
  const auto d1 = xag1.create_pi();
  const auto f11 = xag1.create_and( a1, b1 );
  const auto f21 = xag1.create_and( c1, f11 );
  const auto f31 = xag1.create_and( d1, f21 );
  const auto f41 = xag1.create_xor( f31, f21 );
  xag1.create_po( f41 );

  xag_network xag2;

  const auto a2 = xag2.create_pi();
  const auto b2 = xag2.create_pi();
  const auto c2 = xag2.create_pi();
  const auto f12 = xag2.create_and( b2, c2 );
  const auto f22 = xag2.create_and( a2, f12 );
  xag2.create_po( f22 );

  simulation_cec_stats st;
  CHECK( !*simulation_cec( xag1, xag2, &st ) );
}

TEST_CASE( "CEC different #POs", "[cec]" )
{
  aig_network aig1;

  const auto a1 = aig1.create_pi();
  const auto b1 = aig1.create_pi();
  const auto c1 = aig1.create_pi();
  const auto d1 = aig1.create_pi();
  const auto f11 = aig1.create_and( a1, b1 );
  const auto f21 = aig1.create_and( c1, f11 );
  const auto f31 = aig1.create_and( d1, f21 );
  aig1.create_po( f21 );
  aig1.create_po( f31 );

  aig_network aig2;

  const auto a2 = aig2.create_pi();
  const auto b2 = aig2.create_pi();
  const auto c2 = aig2.create_pi();
  const auto f12 = aig2.create_and( b2, c2 );
  const auto f22 = aig2.create_and( a2, f12 );
  aig2.create_po( f22 );

  simulation_cec_stats st;
  CHECK( !*simulation_cec( aig1, aig2, &st ) );
}

TEST_CASE( "CEC too many PIs", "[cec]" )
{
  aig_network aig1;

  auto const result = lorina::read_aiger( fmt::format( "{}/adder.aig", BENCHMARKS_PATH ), aiger_reader( aig1 ) );
  CHECK( result == lorina::return_code::success );

  simulation_cec_stats st;
  CHECK( simulation_cec( aig1, aig1, &st ) == std::nullopt );
}

TEST_CASE( "CEC different AIGs", "[cec]" )
{
  aig_network aig1;

  const auto a1 = aig1.create_pi();
  const auto b1 = aig1.create_pi();
  const auto f1 = aig1.create_or( a1, b1 );
  aig1.create_po( f1 );

  aig_network aig2;

  const auto a2 = aig2.create_pi();
  const auto b2 = aig2.create_pi();
  const auto f2 = aig2.create_xor( a2, b2 );
  aig2.create_po( f2 );

  simulation_cec_stats st;
  CHECK( !*simulation_cec( aig1, aig2, &st ) );
  CHECK( st.split_var == 2 );
  CHECK( st.rounds == 1 );
}

TEST_CASE( "CEC on optimized design 2", "[cec]" )
{
  xag_network xag1;

  auto const result = lorina::read_aiger( fmt::format( "{}/int2float.aig", BENCHMARKS_PATH ), aiger_reader( xag1 ) );
  CHECK( result == lorina::return_code::success );

  xag_npn_resynthesis<xag_network> resyn;

  exact_library<xag_network, xag_npn_resynthesis<xag_network>> lib( resyn );

  xag_network xag2 = map( xag1, lib );

  simulation_cec_stats st;
  CHECK( *simulation_cec( xag1, xag2, &st ) );
  CHECK( st.split_var == 11 );
  CHECK( st.rounds == 1 );
}

/*
TEST_CASE( "CEC on optimized design 1", "[cec]" )
{
  mig_network mig1;

  auto const result = lorina::read_aiger( fmt::format( "{}/sin.aig", BENCHMARKS_PATH ), aiger_reader( mig1 ) );
  CHECK( result == lorina::return_code::success );

  mig_npn_resynthesis resyn{ true };

  exact_library<mig_network, mig_npn_resynthesis> lib( resyn );

  mig_network mig2 = map( mig1, lib );

  simulation_cec_stats st;
  CHECK( *simulation_cec( mig1, mig2, &st ) );
  CHECK( st.split_var == 18 );
  CHECK( st.rounds == 64 );
}

TEST_CASE( "CEC on badly optimized design 1", "[cec]" )
{
  mig_network mig1;

  auto const result = lorina::read_aiger( fmt::format( "{}/sin.aig", BENCHMARKS_PATH ), aiger_reader( mig1 ) );
  CHECK( result == lorina::return_code::success );

  mig_npn_resynthesis resyn{ true };

  exact_library<mig_network, mig_npn_resynthesis> lib( resyn );

  mig_network mig2 = map( mig1, lib );

  auto a = mig2.pi_at( 22 );
  auto b = mig2.pi_at( 23 );
  auto f = mig2.create_and( mig2.make_signal( a ), mig2.make_signal( b ) );
  auto g = mig2.create_or( f, mig2.make_signal( mig2.index_to_node( 4331 ) ) );
  mig2.replace_in_outputs( mig2.index_to_node( 4331 ), g );

  simulation_cec_stats st;
  CHECK( !*simulation_cec( mig1, mig2, &st ) );
  CHECK( st.split_var == 18 );
  CHECK( st.rounds == 64 );
}
*/


