/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
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
  \file sop_utils.hpp
  \brief Utilities for sum-of-products

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <cstdint>
#include <vector>
#include <cassert>
extern "C"
{
#include <espresso.h>
}


namespace mockturtle
{

enum class sop_simplify_type : uint32_t
{
  simp_espresso = 0,
  simp_exact = 1
};

namespace detail
{

inline bool cube_has_lit( uint64_t const cube, uint64_t const lit )
{
  return ( cube & ( static_cast<uint64_t>( 1 ) << lit ) ) > 0;
}

inline uint32_t cube_count_literals( uint64_t cube )
{
  uint32_t count;
  for ( count = 0; cube; ++count ) { cube &= cube - 1u; };
  return count;
}

inline int64_t sop_literals_occurrences( std::vector<uint64_t> const& sop, uint32_t const num_lit )
{
  /* find the first literal which occurs more than once */
  for ( uint64_t lit = 0; lit < num_lit; ++lit )
  {
    unsigned occurrences = 0;
    for ( auto const& cube : sop )
    {
      if ( cube_has_lit( cube, lit ) )
        ++occurrences;

      if ( occurrences > 1 )
        return lit;
    }
  }

  /* each literal appears once */
  return -1;
}

inline int64_t sop_least_occurrent_literal( std::vector<uint64_t> const& sop, uint32_t const num_lit )
{
  int64_t min_lit = -1;
  uint32_t min_occurrences = UINT32_MAX;

  /* find the first literal which occurs more than once */
  for ( uint64_t lit = 0; lit < num_lit; ++lit )
  {
    uint32_t occurrences = 0;
    for ( auto const& cube : sop )
    {
      if ( cube_has_lit( cube, lit ) )
        ++occurrences;
    }

    if ( occurrences > 1 && occurrences < min_occurrences )
    {
      min_lit = static_cast<int64_t>( lit );
      min_occurrences = occurrences;
    }
  }

  return min_lit;
}

inline int64_t sop_most_occurrent_literal_masked( std::vector<uint64_t> const& sop, uint64_t const cube, uint32_t const num_lit )
{
  int64_t max_lit = -1;
  uint32_t max_occurrences = 1;

  for ( uint64_t lit = 0; lit < num_lit; ++lit )
  {
    if ( !cube_has_lit( cube, lit ) )
      continue;

    uint32_t occurrences = 0;
    for ( auto const& c : sop )
    {
      if ( cube_has_lit( c, lit ) )
        ++occurrences;
    }

    if ( occurrences > max_occurrences )
    {
      max_lit = static_cast<int64_t>( lit );
      max_occurrences = occurrences;
    }
  }

  return max_lit;
}

inline bool sop_maximal_cube_literal( std::vector<uint64_t> const& sop, uint64_t& cube, uint32_t const lit )
{
  uint64_t max_cube = UINT64_MAX;
  uint32_t occurrences = 0;

  for ( auto const& c : sop )
  {
    if ( cube_has_lit( c, lit ) )
    {
      ++occurrences;
      max_cube &= c;
    }
  }

  cube = max_cube;
  return occurrences > 1;
}

inline void sop_best_literal( std::vector<uint64_t> const& sop, std::vector<uint64_t>& result, uint64_t const cube, uint32_t const num_lit )
{
  int64_t max_lit = sop_most_occurrent_literal_masked( sop, cube, num_lit );
  assert( max_lit >= 0 );

  result.push_back( static_cast<uint64_t>( 1 ) << max_lit );
}

inline pset_family sop_to_espresso( std::vector<uint64_t> const& sop, uint32_t const num_vars )
{
  pset_family cover;
  pset set;

  assert( cube.size == 2 * num_vars );

  if ( sop.size() == 0 )
  {
    /* const 0 */
    cover = sf_new( 0, cube.size );
    return cover;
  }
  else if ( sop.size() == 1 && sop[0] == 0 )
  {
    /* const 1 */
    cover = sf_new( 1, cube.size );
    set = GETSET( cover, cover->count++ );
    set_copy( set, cube.fullset );
    return cover;
  }

  cover = sf_new( sop.size(), cube.size );
  for ( auto const& c : sop )
  {
    set = GETSET( cover, cover->count++ );
    set_copy( set, cube.fullset );
    for ( int i = 0; i < num_vars; ++i )
    {
      if ( ( ( c >> ( 2 * i ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i + 1 );
      }
      else if ( ( ( c >> ( 2 * i + 1 ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i );
      }
    }
  }

  return cover;
}

inline std::vector<uint64_t> espresso_to_sop( pset_family cover )
{
  std::vector<uint64_t> sop;
  pset set;

  uint32_t num_vars = cover->sf_size >> 1;

  for ( int i = 0; i < cover->count; ++i )
  {
    uint64_t c = 0;
    set = GETSET( cover, i );
    for ( int j = 0; j < num_vars; ++j )
    {
      int lit = GETINPUT( set, j );

      if ( lit == ESPRESSO_ZERO )
      {
        c |= 1 << ( 2 * j );
      }
      else if ( lit == ESPRESSO_ONE )
      {
        c |= 1 << ( 2 * j + 1 );
      }
    }
    sop.push_back( c );
  }

  return sop;
}

inline pset_family espresso_minimize( pset_family onset, pset_family dcset )
{
  pset_family offset;
  bool dcset_tofree;

  /* create the dcset */
  dcset_tofree = ( dcset == NULL );
  if ( dcset == NULL )
      dcset = sf_new( 1, cube.size );

  dcset->wsize = onset->wsize;
  dcset->sf_size = onset->sf_size;

  // derive the offset
  if ( dcset->sf_size == 0 || dcset->count == 0 )
      offset = complement( cube1list( onset ) );
  else
      offset = complement( cube2list( onset, dcset ) ); 

  // perform minimization
  onset = espresso( onset, dcset, offset );

  // free covers
  sf_free( offset );

  if ( dcset_tofree )
    sf_free( dcset );

  return onset;
}

inline pset_family exact_minimize( pset_family onset, pset_family dcset )
{
  pset_family offset;
  bool dcset_tofree;

  /* create the dcset */
  dcset_tofree = ( dcset == NULL );
  if ( dcset == NULL )
      dcset = sf_new( 1, cube.size );

  dcset->wsize = onset->wsize;
  dcset->sf_size = onset->sf_size;

  // derive the offset
  if ( dcset->sf_size == 0 || dcset->count == 0 )
      offset = complement( cube1list( onset ) );
  else
      offset = complement( cube2list( onset, dcset ) ); 

  // perform minimization
  onset = minimize_exact( onset, dcset, offset, 1 );

  // free covers
  sf_free( offset );

  if ( dcset_tofree )
    sf_free( dcset );

  return onset;
}

} // namespace detail

inline uint32_t sop_count_literals( std::vector<uint64_t> const& sop )
{
  uint32_t lit_count = 0;

  for ( auto const& c : sop )
  {
    lit_count += detail::cube_count_literals( c );
  }

  return lit_count;
}

/*! \brief Makes a SOP cube free
*
* This method checks for a common cube divisor in
* the SOP. If found, the SOP is divided by that cube.
* 
* \param sop
*/
inline void sop_make_cube_free( std::vector<uint64_t>& sop )
{
  /* find common cube */
  uint64_t mask = UINT64_MAX;
  for ( auto const& c : sop )
  {
    mask &= c;
  }

  if ( mask == 0 )
    return;
  
  /* make cube free */
  for ( auto& c : sop )
  {
    c &= ~mask;
  }
}

/*! \brief Checks if a SOP is cube free
*
* This method checks for a common cube divisor in
* the SOP.
* 
* \param sop
*/
inline bool sop_is_cube_free( std::vector<uint64_t> const& sop )
{
  /* find common cube */
  uint64_t mask = UINT64_MAX;
  for ( auto const& c : sop )
  {
    mask &= c;
  }

  return mask == 0;
}

/*! \brief Algebraic division by literal
*
* This method divides a SOP inplace by a literal and
* stores the resulting quotient in the original SOP.
* 
* \param sop
* \param lit
*/
inline void sop_divide_by_literal_inplace( std::vector<uint64_t>& sop, uint64_t const lit )
{
  uint32_t p = 0;
  for ( auto i = 0; i < sop.size(); ++i )
  {
    if ( detail::cube_has_lit( sop[i], lit ) )
    {
      sop[p++] = sop[i] & ( ~( static_cast<uint64_t>( 1 ) << lit ) );
    }
  }

  sop.resize( p );
}

/*! \brief Algebraic division by a cube
*
* This method divides a SOP (divident) by the divisor
* and stores the resulting quotient and reminder.
* 
* \param Divident
* \param Divisor
* \param Quotient
* \param Reminder
*/
inline void sop_divide_by_cube( std::vector<uint64_t> const& divident, std::vector<uint64_t> const& divisor, std::vector<uint64_t>& quotient, std::vector<uint64_t>& reminder )
{
  assert( divisor.size() == 1 );

  quotient.clear();
  reminder.clear();

  uint32_t p = 0;
  for ( auto const& c : divident )
  {
    if ( ( c & divisor[0] ) == divisor[0] )
    {
      quotient.push_back( c & ( ~divisor[0] ) );
    }
    else
    {
      reminder.push_back( c );
    }
  }
}

/*! \brief Algebraic division by a cube
*
* This method divides a SOP (divident) by the divisor
* and stores the resulting quotient.
* 
* \param Divident
* \param Divisor
* \param Quotient
*/
inline void sop_divide_by_cube_no_reminder( std::vector<uint64_t> const& divident, uint64_t const& divisor, std::vector<uint64_t>& quotient )
{
  quotient.clear();

  uint32_t p = 0;
  for ( auto const& c : divident )
  {
    if ( ( c & divisor ) == divisor )
    {
      quotient.push_back( c & ~divisor );
    }
  }
}

/*! \brief Algebraic division
*
* This method divides a SOP (divident) by the divisor
* and stores the resulting quotient and reminder.
* 
* \param Divident
* \param Divisor
* \param Quotient
* \param Reminder
*/
inline void sop_divide( std::vector<uint64_t>& divident, std::vector<uint64_t> const& divisor, std::vector<uint64_t>& quotient, std::vector<uint64_t>& reminder )
{
  /* divisor contains a single cube */
  if ( divisor.size() == 1 )
  {
    sop_divide_by_cube( divident, divisor, quotient, reminder );
    return;
  }

  quotient.clear();
  reminder.clear();

  /* perform division */
  for ( auto i = 0; i < divident.size(); ++i )
  {
    auto const c = divident[i];

    /* cube has been already covered */
    if ( detail::cube_has_lit( c, 63 ) )
      continue;

    uint32_t div_i;
    for ( div_i = 0u; div_i < divisor.size(); ++div_i )
    {
      if ( ( c & divisor[div_i] ) == divisor[div_i] )
        break;
    }

    /* cube is not divisible -> reminder */
    if ( div_i >= divisor.size() )
      continue;

    /* extract quotient */
    uint64_t c_quotient = c & ~divisor[div_i];

    /* find if c_quotient can be obtained for all the divisors */
    bool found = true;
    for ( auto const& div : divisor )
    {
      if ( div == divisor[div_i] )
        continue;

      found = false;
      for ( auto const& c2 : divident )
      {
        /* cube has been already covered */
        if ( detail::cube_has_lit( c2, 63 ) )
          continue;

        if ( ( ( c2 & div ) == div ) && ( c_quotient == ( c2 & ~div ) ) )
        {
          found = true;
          break;
        }
      }

      if ( !found )
        break;
    }

    if ( !found )
      continue;
    
    /* valid divisor, select covered cubes */
    quotient.push_back( c_quotient );

    divident[i] |= static_cast<uint64_t>( 1 ) << 63;
    for ( auto const& div : divisor )
    {
      if ( div == divisor[div_i] )
        continue;

      for ( auto& c2 : divident )
      {
        /* cube has been already covered */
        if ( detail::cube_has_lit( c2, 63 ) )
          continue;

        if ( ( ( c2 & div ) == div ) && ( c_quotient == ( c2 & ~div ) ) )
        {
          c2 |= static_cast<uint64_t>( 1 ) << 63;
          break;
        }
      }
    }
  }

  /* add remainder */
  for ( auto& c : divident )
  {
    if ( !detail::cube_has_lit( c, 63 ) )
    {
      reminder.push_back( c );
    }
    else
    {
      /* unmark */
      c &= ~( static_cast<uint64_t>( 1 ) << 63 );
    }
  }
}

/*! \brief Boolean division
*
* This method divides a SOP (divident) by the divisor
* and stores the resulting quotient and reminder.
* 
* \param Divident
* \param Divisor
* \param Quotient
* \param Reminder
*/
inline void sop_bool_divide( std::vector<uint64_t>& divident, std::vector<uint64_t> const& divisor, std::vector<uint64_t>& quotient, std::vector<uint64_t>& reminder, uint32_t const num_lit )
{
  /* algebraic division */
  sop_divide( divident, divisor, quotient, reminder );

  /* compute don't care set g*x' + g'*x */
  define_cube_size( static_cast<int>( num_lit / 2 + 1 ) );
  pset_family divisor_set = detail::sop_to_espresso( divisor, num_lit / 2 + 1 );
  pset_family n_divisor_set = complement( cube1list( divisor_set ) );
  pset_family dcset = sf_new( divisor_set->count + n_divisor_set->count, cube.size );
  pset set, set_dc;

  /* add g*x' */
  for ( int i = 0; i < divisor_set->count; ++i )
  {
    set = GETSET( divisor_set, i );
    set_dc = GETSET( dcset, dcset->count++ );
    set_copy( set_dc, set );
    set_remove( set_dc, num_lit + 1 );
  }

  /* add g'*x */
  for ( int i = 0; i < n_divisor_set->count; ++i )
  {
    set = GETSET( n_divisor_set, i );
    set_dc = GETSET( dcset, dcset->count++ );
    set_copy( set_dc, set );
    set_remove( set_dc, num_lit );
  }

  /* compute onset: h*x + e */
  pset_family onset = sf_new( quotient.size() + reminder.size(), cube.size );

  /* add h*x */
  for ( auto const& c : quotient )
  {
    set = GETSET( onset, onset->count++ );
    set_copy( set, cube.fullset );
    for ( int i = 0; i < num_lit / 2; ++i )
    {
      if ( ( ( c >> ( 2 * i ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i + 1 );
      }
      else if ( ( ( c >> ( 2 * i + 1 ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i );
      }
    }
    set_remove( set, num_lit );
  }

  /* add e */
  for ( auto const& c : reminder )
  {
    set = GETSET( onset, onset->count++ );
    set_copy( set, cube.fullset );
    for ( int i = 0; i < num_lit / 2; ++i )
    {
      if ( ( ( c >> ( 2 * i ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i + 1 );
      }
      else if ( ( ( c >> ( 2 * i + 1 ) ) & 1 ) == 1 )
      {
        set_remove( set, 2 * i );
      }
    }
  }

  pset_family offset;
  offset = complement( cube2list( onset, dcset ) );
  onset = complement( cube2list( offset, dcset ) );

  /* smoothing x' from the onset */
  onset = minimize_exact( onset, dcset, offset, 0 );

  /* TODO: bug fix reminder could contains x' (return algebraic division if this is the case) */
  for ( int i = 0; i < onset->count; ++i )
  {
    set = GETSET( onset, i );
    int lit = GETINPUT( set, num_lit / 2 );
    if ( lit == ESPRESSO_ZERO )
    {
      sf_free( offset );
      sf_free( onset );
      sf_free( dcset );
      sf_free( n_divisor_set );
      sf_free( divisor_set );
      return;
    }
  }

  auto f = detail::espresso_to_sop( onset );

  /* extract quotient and reminder */
  sop_divide_by_cube( f, {static_cast<uint64_t>( 1 << ( num_lit + 1) )}, quotient, reminder );

  sf_free( offset );
  sf_free( onset );
  sf_free( dcset );
  sf_free( n_divisor_set );
  sf_free( divisor_set );
}

/*! \brief Extracts all the kernels
*
* This method is used to identify and collect all
* the kernels.
*
* \param sop
* \param kernels
* \param j
* \param num_lit
*/
inline void sop_kernels_rec( std::vector<uint64_t> const& sop, std::vector<std::vector<uint64_t>>& kernels, uint32_t const j, uint32_t const num_lit )
{
  std::vector<uint64_t> kernel;

  for ( uint32_t i = j; i < num_lit; ++i )
  {
    uint64_t c;
    if ( detail::sop_maximal_cube_literal( sop, c, i ) )
    {
      /* cube has been visited already */
      if ( ( c & ( ( static_cast<uint64_t>( 1 ) << i ) - 1 ) ) > 0u )
        continue;

      sop_divide_by_cube_no_reminder( sop, c, kernel );
      sop_kernels_rec( kernel, kernels, i + 1, num_lit );
    }
  }

  kernels.push_back( sop );
}

/*! \brief Extracts the best factorizing kernel
*
* This method is used to identify the best kernel
* according to the algebraic factorization value.
*
* \param sop
* \param kernel
* \param best_kernel
* \param j
* \param best_cost
* \param num_lit
*/
inline uint32_t sop_best_kernel_rec( std::vector<uint64_t>& sop, std::vector<uint64_t>& kernel, std::vector<uint64_t>& best_kernel, uint32_t const j, uint32_t& best_cost, uint32_t const num_lit )
{
  std::vector<uint64_t> new_kernel;
  std::vector<uint64_t> quotient;
  std::vector<uint64_t> reminder;

  /* evaluate kernel */
  sop_divide( sop, kernel, quotient, reminder );
  uint32_t division_cost = sop_count_literals( quotient ) + sop_count_literals( reminder );
  uint32_t best_fact_cost = sop_count_literals( kernel );

  for ( uint32_t i = j; i < num_lit; ++i )
  {
    uint64_t c;
    if ( detail::sop_maximal_cube_literal( kernel, c, i ) )
    {
      /* cube has been visited already */
      if ( ( c & ( ( static_cast<uint64_t>( 1 ) << i ) - 1 ) ) > 0u )
        continue;

      /* extract the new kernel */
      sop_divide_by_cube( kernel, {c}, new_kernel, reminder );
      uint32_t fact_cost_rec = detail::cube_count_literals( c ) + sop_count_literals( reminder );
      uint32_t fact_cost = sop_best_kernel_rec( sop, new_kernel, best_kernel, i + 1, best_cost, num_lit );

      /* compute the factorization value for kernel */
      if ( ( fact_cost + fact_cost_rec ) < best_fact_cost )
        best_fact_cost = fact_cost + fact_cost_rec;
    }
  }

  if ( best_kernel.empty() || ( division_cost + best_fact_cost ) < best_cost )
  {
    best_kernel = kernel;
    best_cost = division_cost + best_fact_cost;
  }

  return best_fact_cost;
}

/*! \brief Extracts the best factorizing kernel
*
* This method is used to identify the best kernel
* according to the Boolean factorization value.
*
* \param sop
* \param kernel
* \param best_kernel
* \param j
* \param best_cost
* \param num_lit
*/
inline uint32_t sop_best_kernel_bool_rec( std::vector<uint64_t>& sop, std::vector<uint64_t>& kernel, std::vector<uint64_t>& best_kernel, uint32_t const j, uint32_t& best_cost, uint32_t const num_lit )
{
  std::vector<uint64_t> new_kernel;
  std::vector<uint64_t> quotient;
  std::vector<uint64_t> reminder;

  /* evaluate kernel */
  sop_bool_divide( sop, kernel, quotient, reminder, num_lit );
  uint32_t division_cost = sop_count_literals( quotient ) + sop_count_literals( reminder );
  uint32_t best_fact_cost = sop_count_literals( kernel );

  for ( uint32_t i = j; i < num_lit; ++i )
  {
    uint64_t c;
    if ( detail::sop_maximal_cube_literal( kernel, c, i ) )
    {
      /* cube has been visited already */
      if ( ( c & ( ( static_cast<uint64_t>( 1 ) << i ) - 1 ) ) > 0u )
        continue;

      /* extract the new kernel */
      sop_divide_by_cube( kernel, {c}, new_kernel, reminder );
      uint32_t fact_cost_rec = detail::cube_count_literals( c ) + sop_count_literals( reminder );
      uint32_t fact_cost = sop_best_kernel_bool_rec( sop, new_kernel, best_kernel, i + 1, best_cost, num_lit );

      /* compute the factorization value for kernel */
      if ( ( fact_cost + fact_cost_rec ) < best_fact_cost )
        best_fact_cost = fact_cost + fact_cost_rec;
    }
  }

  if ( best_kernel.empty() || ( division_cost + best_fact_cost ) < best_cost )
  {
    best_kernel = kernel;
    best_cost = division_cost + best_fact_cost;
  }

  return best_fact_cost;
}

/*! \brief Extracts a one level-0 kernel
*
* This method is used to identify and return a one
* level-0 kernel.
* 
* \param sop
* \param num_lit
*/
inline void sop_one_level_zero_kernel_rec( std::vector<uint64_t>& sop, uint32_t const num_lit )
{
  /* find least occurring leteral which occurs more than once. TODO: test other metrics */
  int64_t min_lit = detail::sop_least_occurrent_literal( sop, num_lit );

  if ( min_lit == -1 )
    return;
  
  sop_divide_by_literal_inplace( sop, static_cast<uint64_t>( min_lit ) );
  sop_make_cube_free( sop );

  sop_one_level_zero_kernel_rec( sop, num_lit );
}

/*! \brief Finds a quick divisor for a SOP
*
* This method is used to identify and return a quick
* divisor for the SOP.
* 
* \param sop
* \param num_lit
*/
inline bool sop_quick_divisor( std::vector<uint64_t> const& sop, std::vector<uint64_t>& res, uint32_t const num_lit )
{
  if ( sop.size() <= 1 )
    return false;

  /* each literal appears no more than once */
  if ( detail::sop_literals_occurrences( sop, num_lit ) < 0 )
    return false;

  /* one level 0-kernel */
  res = sop;
  sop_one_level_zero_kernel_rec( res, num_lit );

  assert( res.size() );
  return true;
}

/*! \brief Finds a good divisor for a SOP
*
* This method is used to identify and return a good
* divisor for the SOP.
* 
* \param sop
* \param num_lit
*/
inline bool sop_good_divisor( std::vector<uint64_t>& sop, std::vector<uint64_t>& res, uint32_t const num_lit )
{
  if ( sop.size() <= 1 )
    return false;

  /* each literal appears no more than once */
  if ( detail::sop_literals_occurrences( sop, num_lit ) < 0 )
    return false;

  std::vector<uint64_t> kernel = sop;

  /* compute all the kernels and return the one with the best factorization value */
  uint32_t best_cost = 0;
  sop_best_kernel_rec( sop, kernel, res, 0, best_cost, num_lit );

  return true;
}

/*! \brief Finds a good Boolean divisor for a SOP
*
* This method is used to identify and return a good
* Boolean divisor for the SOP.
* 
* \param sop
* \param num_lit
*/
inline bool sop_good_divisor_bool( std::vector<uint64_t>& sop, std::vector<uint64_t>& res, uint32_t const num_lit )
{
  if ( sop.size() <= 1 )
    return false;

  /* each literal appears no more than once */
  if ( detail::sop_literals_occurrences( sop, num_lit ) < 0 )
    return false;

  std::vector<uint64_t> kernel = sop;

  /* compute all the kernels and return the one with the best factorization value */
  uint32_t best_cost = 0;
  sop_best_kernel_bool_rec( sop, kernel, res, 0, best_cost, num_lit );

  return true;
}

/*! \brief Translates cubes into products
*
* This method translate SOP of kitty::cubes (bits + mask)
* into SOP of products represented by literals.
* Example for: ab'c*
*   - cube: _bits = 1010; _mask = 1110
*   - product: 10011000
* 
* \param cubes Sum-of-products described using cubes
* \param num_vars Number of variables
*/
inline std::vector<uint64_t> cubes_to_sop( std::vector<kitty::cube> const& cubes, uint32_t const num_vars )
{
  using sop_t = std::vector<uint64_t>;

  sop_t sop( cubes.size() );

  /* Represent literals instead of variables as a a' b b' */
  /* Bit 63 is reserved, up to 31 varibles support */
  auto it = sop.begin();
  for ( auto const& c : cubes )
  {
    uint64_t& product = *it++;
    for ( auto i = 0; i < num_vars; ++i )
    {
      if ( c.get_mask( i ) )
        product |= static_cast<uint64_t>( 1 ) << ( 2*i + static_cast<unsigned>( c.get_bit( i ) ) );
    }
  }

  return sop;
}

inline void minimize_sop( std::vector<uint64_t>& sop, uint32_t const num_vars, sop_simplify_type type = sop_simplify_type::simp_espresso )
{
  pset_family cover;

  define_cube_size( static_cast<int>( num_vars ) );
  cover = detail::sop_to_espresso( sop, num_vars );

  if ( type == sop_simplify_type::simp_espresso )
    cover = detail::espresso_minimize( cover, NULL );
  else
    cover = detail::exact_minimize( cover, NULL );

  sop = detail::espresso_to_sop( cover );

  sf_free( cover );
}

} //namespace mockturtle