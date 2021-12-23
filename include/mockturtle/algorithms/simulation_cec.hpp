/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
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

/*!
  \file simulation_cec.hpp
  \brief Simulation-based CEC

  EPFL CS-472 2021 Final Project Option 2
*/

#pragma once

#include <math.h>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"

namespace mockturtle
{

/* Statistics to be reported */
struct simulation_cec_stats
{
  /*! \brief Split variable (simulation size). */
  uint32_t split_var{ 0 };

  /*! \brief Number of simulation rounds. */
  uint32_t rounds{ 0 };
};

namespace detail
{

template<class Ntk>
class simulation_cec_impl
{
public:
  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ),
        _st( st )
  {
  }

  bool run()
  {
    /* TODO: write your implementation here */

    // Step 1: update simulation_cec_stats _st :
    uint32_t num_nodes = _ntk.size();     // V
    uint32_t num_inputs = _ntk.num_pis(); // n

    std::cout << "\[simulation cec\] start new run ------------------------" << std::endl;
    std::cout << "\[simulation cec\] num_inputs is " << num_inputs << std::endl;

    assert( num_inputs <= 40u );

    if ( num_inputs <= 6u )
    {
      _st.split_var = num_inputs;
    }
    else
    {
      uint32_t tmp = floor( log2( pow( 2, 29 ) / num_nodes - 32 ) + 3 );

      // std::cout << "\[simulation cec\] tmp is " << tmp << std::endl;

      _st.split_var = ( num_inputs <= tmp ? num_inputs : tmp );
    }
    _st.rounds = ( 1u << ( num_inputs - _st.split_var ) );

    // For test only
    // _st.split_var = 7u;
    // _st.rounds = ( 1u << ( num_inputs - _st.split_var ) ); 

    std::cout << "\[simulation cec\] split_var is " << _st.split_var << std::endl;
    std::cout << "\[simulation cec\] rounds is " << _st.rounds << std::endl;

    // Step 2: initialize patters (pattern_t), simulator :
    pattern_t patterns( _ntk );
    default_simulator<kitty::dynamic_truth_table> sim( num_inputs );
    bool result = true;

    // Step 3: iterate
    if ( _st.rounds == 1u )
    {
      patterns.reset();
      simulate_nodes<kitty::dynamic_truth_table>( _ntk, patterns, sim );
      _ntk.foreach_po( [&]( signal const& f, uint32_t j )
                       {
                         std::vector<uint64_t> tt_bits = ( _ntk.is_complemented( f ) ? ~patterns[f] : patterns[f] )._bits;
                         // std::cout << "\t\t j = " << j << "; size of tt =  " << tt_bits.size() << std::endl;

                         for ( std::vector<uint64_t>::iterator it = tt_bits.begin(); it != tt_bits.end(); ++it )
                         {
                           // std::cout << "\[simulation cec\] result is " << *it << std::endl;
                           if ( *it != 0x0 )
                           {
                             result = result && false;
                             std::cout << "\[simulation cec\] inconsistency detected." << std::endl;
                           }
                         }
                       } );
    }
    else
    {
      for ( uint32_t round = 0u; round < _st.rounds; round = round + 1u )
      {
        std::cout << "\[simulation cec\] round-1 is " << round << std::endl;

        patterns.reset();

        uint32_t split_var = _st.split_var;

        _ntk.foreach_pi( [&]( node const& n, uint32_t i )
                         {
                           if ( i >= split_var )
                           {
                             uint32_t index = i - split_var;
                             // std::cout << "\[simulation cec\] - update_patterns, index is " << index << std::endl;

                             std::cout << "\[simulation cec\] - update_patterns, update to " << (( round >> index ) % 2u) << std::endl;
                             if ( ( round >> index ) % 2u )
                             { // if this input should be set as 1:
                               patterns[n] = ~kitty::dynamic_truth_table( _ntk.num_pis() );
                             }
                             else
                             { // if this input should be set as 0:
                               patterns[n] = kitty::dynamic_truth_table( _ntk.num_pis() );
                             }
                           }
                         } );
        std::cout << "\[simulation cec\] end of update_patterns" << std::endl;
        
        simulate_nodes<kitty::dynamic_truth_table>( _ntk, patterns, sim );

        _ntk.foreach_po( [&]( signal const& f, uint32_t j )
                         {
                           std::vector<uint64_t> tt_bits = ( _ntk.is_complemented( f ) ? ~patterns[f] : patterns[f] )._bits;
                           for ( std::vector<uint64_t>::iterator it = tt_bits.begin(); it != tt_bits.end(); ++it )
                           {
                             if ( *it != 0x0 )
                             {
                               result = result && false;
                               std::cout << "\[simulation cec\] inconsistency detected." << std::endl;
                             }
                           }
                         } );
        if ( !result )
        {
          return result;
        }
      }
    }
    return result;
  }

private:
  /* you can add additional methods here */

private:
  Ntk& _ntk;
  simulation_cec_stats& _st;
  /* you can add other attributes here */
};

} // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
 *
 * This function implements a simulation-based combinational equivalence checker.
 * The implementation creates a miter network and run several rounds of simulation
 * to verify the functional equivalence. For memory and speed reasons this approach
 * is limited up to 40 input networks. It returns an optional which is `nullopt`,
 * if the network has more than 40 inputs.
 */
template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}

} // namespace mockturtle
