/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
  \file experiments.hpp
  \brief Framework for simple experimental evaluation

  \author Mathias Soeken
*/

#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <fmt/color.h>
#include <fmt/format.h>
#include <mockturtle/io/write_bench.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <nlohmann/json.hpp>

namespace crypto_experiments
{

struct json_table
{
  explicit json_table( nlohmann::json const& data, std::vector<std::string> const& columns )
      : columns_( columns )
  {
    for ( auto const& column : columns )
    {
        max_widths_.push_back( static_cast<uint32_t>( column.size() ) );
      //max_widths_.push_back( column.size() );
    }
    entries_.push_back( columns );
    for ( auto const& row : data )
    {
      add_row( row );
    }
  }

  void print( std::ostream& os )
  {
    for ( const auto& entry : entries_ )
    {
      os << "|";
      for ( auto i = 0u; i < entry.size(); ++i )
      {
        os << fmt::format( " {1:>{0}} |", max_widths_[i], entry[i] );
      }
      os << "\n";
    }
  }

private:
  void add_row( nlohmann::json const& row )
  {
    std::vector<std::string> entry;
    uint32_t ctr{0u};
    for ( auto const& key : columns_ )
    {
      auto const& data = row[key];
      std::string cell;

      if ( data.is_string() )
      {
        cell = static_cast<std::string>( data );
      }
      else if ( data.is_number_integer() )
      {
        cell = std::to_string( static_cast<int>( data ) );
      }
      else if ( data.is_number() )
      {
        cell = fmt::format( "{:.2f}", static_cast<float>( data ) );
      }
      else if ( data.is_boolean() )
      {
        cell = fmt::format( "{}", static_cast<bool>( data ) );
      }

      max_widths_[ctr] = std::max( max_widths_[ctr], static_cast<uint32_t>( cell.size() ) );
      //max_widths_[ctr] = std::max<uint32_t>( max_widths_[ctr], cell.size() );
      ++ctr;
      entry.push_back( cell );
    }
    entries_.push_back( entry );
  }

private:
  std::vector<uint32_t> max_widths_;
  std::vector<std::string> columns_;
  std::vector<std::vector<std::string>> entries_;
};

static constexpr const char* use_github_revision = "##GITHUB##";

template<typename T, typename... Ts>
struct first_type
{
  using type = T;
};

template<typename... Ts>
using first_type_t = typename first_type<Ts...>::type;

template<typename... ColumnTypes>
class experiment
{
public:
  template<typename... T>
  explicit experiment( std::string_view name, T... column_names )
      : name_( name )
  {
    static_assert( ( sizeof...( ColumnTypes ) > 0 ), "at least one column must be specified" );
    static_assert( ( sizeof...( ColumnTypes ) == sizeof...( T ) ), "number of column names must match column types" );
    static_assert( ( std::is_constructible_v<std::string, T> && ... ), "all column names must be strings" );
    ( column_names_.push_back( column_names ), ... );

#ifndef EXPERIMENTS_PATH
    filename_ = fmt::format( "{}.json", name );
#else
    filename_ = fmt::format( "{}{}.json", EXPERIMENTS_PATH, name );
#endif

    std::ifstream in( filename_, std::ifstream::in );
    if ( in.good() )
    {
      data_ = nlohmann::json::parse( in );
    }
  }

  void save( std::string_view version = use_github_revision )
  {
    nlohmann::json entries;
    for ( auto const& row : rows_ )
    {
      auto it = column_names_.begin();
      nlohmann::json entry;
      std::apply(
          [&]( auto&&... args ) {
            ( ( entry[*it++] = args ), ... );
          },
          row );
      entries.push_back( entry );
    }

    std::string version_;
    version_ = version;
//#ifdef GIT_SHORT_REVISION
//    if ( version == crypto_experiments::use_github_revision )
//    {
//      version_ = GIT_SHORT_REVISION;
//    }
//#endif

    if ( !data_.empty() && data_.back()["version"] == version_ )
    {
      data_.erase( data_.size() - 1u );
    }

    data_.push_back( {{"version", version_},
                      {"entries", entries}} );

    std::ofstream os( filename_, std::ofstream::out );
    os << data_.dump( 2 ) << "\n";
  }

  void operator()( ColumnTypes... args )
  {
    rows_.emplace_back( args... );
  }

  nlohmann::json const& dataset( std::string const& version, nlohmann::json const& def ) const
  {
    if ( version.empty() )
    {
      return def;
    }
    else
    {
      if ( const auto it = std::find_if( data_.begin(), data_.end(), [&]( auto const& entry ) { return entry["version"] == version; } ); it != data_.end() )
      {
        return *it;
      }
      else
      {
        throw std::exception();
      }
    }
  }

  bool table( std::string const& version = {}, std::ostream& os = std::cout ) const
  {
    if ( data_.empty() )
    {
      fmt::print( "[w] no data available\n" );
      return false;
    }

    try
    {
      auto const& data = dataset( version, data_.back() );
      fmt::print( "[i] dataset " );
      fmt::print( fg( fmt::terminal_color::blue ), "{}\n", data["version"] );

      json_table( data["entries"], column_names_ ).print( os );
    }
    catch ( ... )
    {
      fmt::print( "[w] version {} not found\n", version );
      return false;
    }

    return true;
  }

  bool compare( std::string const& old_version = {},
                std::string const& current_version = {},
                std::vector<std::string> const& track_columns = {},
                std::ostream& os = std::cout )
  {
    if ( data_.size() < 2u )
    {
      fmt::print( "[w] dataset contains less than two entry sets\n" );
      return false;
    }

    try
    {
      auto const& data_old = dataset( old_version, data_[data_.size() - 2u] );
      auto const& data_cur = dataset( current_version, data_.back() );

      auto const& entries_old = data_old["entries"];
      auto const& entries_cur = data_cur["entries"];

      fmt::print( "[i] compare " );
      fmt::print( fg( fmt::terminal_color::blue ), "{}", data_old["version"] );
      fmt::print( " to " );
      fmt::print( fg( fmt::terminal_color::blue ), "{}\n", data_cur["version"] );

      /* collect keys */
      using first_t = first_type_t<ColumnTypes...>;
      std::vector<first_t> keys;
      for ( auto const& entry : entries_cur )
      {
        nlohmann::json const& j = entry[column_names_.front()];
        keys.push_back( j.get<first_t>() );
      }
      for ( auto const& entry : entries_old )
      {
        nlohmann::json const& j = entry[column_names_.front()];
        auto value = j.get<first_t>();
        if ( std::find( keys.begin(), keys.end(), value ) == keys.end() )
        {
          keys.push_back( value );
        }
      }

      /* track differences */
      std::unordered_map<std::string, uint32_t> differences;
      for ( auto const& column : track_columns )
      {
        differences[column] = 0u;
      }

      /* prepare entries */
      auto find_key = [&]( nlohmann::json const& entries, first_t const& key ) {
        return std::find_if( entries.begin(), entries.end(), [&]( auto const& entry ) {
          nlohmann::json const& j = entry[column_names_.front()];
          return j.get<first_t>() == key;
        } );
      };

      auto compare_columns = column_names_;
      std::transform( column_names_.begin() + 1, column_names_.end(), std::back_inserter( compare_columns ),
                      []( auto const& name ) { return name + "'"; } );

      nlohmann::json compare_entries;
      for ( auto const& key : keys )
      {
        nlohmann::json row;
        const auto it_old = find_key( entries_old, key );
        if ( it_old != entries_old.end() )
        {
          row = *it_old;
        }
        if ( auto const it = find_key( entries_cur, key ); it != entries_cur.end() )
        {
          if ( it_old == entries_old.end() )
          {
            row[column_names_[0]] = (*it)[column_names_[0]];
          }
          for ( auto i = 1u; i < column_names_.size(); ++i )
          {
            row[column_names_[i] + "'"] = (*it)[column_names_[i]];

            if ( it_old != entries_old.end() )
            {
              if ( const auto it_diff = differences.find( column_names_[i] ); it_diff != differences.end() && row[column_names_[i]] != row[column_names_[i] + "'"] )
              {
                it_diff->second++;
              }
            }
          }
        }
        compare_entries.push_back( row );
      }

      json_table( compare_entries, compare_columns ).print( os );

      for ( const auto& [k, v] : differences ) {
        if ( v == 0u ) {
          os << fmt::format( "[i] no differences in column '{}'\n", k );
        } else {
          os << fmt::format( "[i] {} differences in column '{}'\n", v, k );
        }
      }
    }
    catch ( ... )
    {
      fmt::print( "[w] dataset not found\n" );
      return false;
    }

    return true;
  }

private:
  std::string name_;
  std::string filename_;
  std::vector<std::string> column_names_;
  std::vector<std::tuple<ColumnTypes...>> rows_;

  nlohmann::json data_;
};

// clang-format off
static constexpr uint32_t adder      = 0b00000000000000000001;
static constexpr uint32_t bar        = 0b00000000000000000010;
static constexpr uint32_t div        = 0b00000000000000000100;
static constexpr uint32_t hyp        = 0b00000000000000001000;
static constexpr uint32_t log2       = 0b00000000000000010000;
static constexpr uint32_t max        = 0b00000000000000100000;
static constexpr uint32_t multiplier = 0b00000000000001000000;
static constexpr uint32_t sin        = 0b00000000000010000000;
static constexpr uint32_t sqrt       = 0b00000000000100000000;
static constexpr uint32_t square     = 0b00000000001000000000;
static constexpr uint32_t arbiter    = 0b00000000010000000000;
static constexpr uint32_t cavlc      = 0b00000000100000000000;
static constexpr uint32_t ctrl       = 0b00000001000000000000;
static constexpr uint32_t dec        = 0b00000010000000000000;
static constexpr uint32_t i2c        = 0b00000100000000000000;
static constexpr uint32_t int2float  = 0b00001000000000000000;
static constexpr uint32_t mem_ctrl   = 0b00010000000000000000;
static constexpr uint32_t priority   = 0b00100000000000000000;
static constexpr uint32_t router     = 0b01000000000000000000;
static constexpr uint32_t voter      = 0b10000000000000000000;
static constexpr uint32_t arithmetic = 0b00000000001111111111;
static constexpr uint32_t random     = 0b11111111110000000000;
static constexpr uint32_t all        = 0b11111111111111111111;
// clang-format on

static const char* epfl_benchmark_names[] = {
    "adder", "bar", "div", "hyp", "log2", "max", "multiplier", "sin", "sqrt", "square",
    "arbiter", "cavlc", "ctrl", "dec", "i2c", "int2float", "mem_ctrl", "priority", "router", "voter"};

static const char* crypto_benchmark_names[] = {
  "AES-expanded_untilsat",
  "AES-non-expanded_unstilsat",
  "DES-expanded_untilsat",
  "DES-non-expanded_untilsat",
  "adder_32bit_untilsat",
  "adder_64bit_untilsat",
  "adder_untilsat",
  // "arbiter_untilsat",
  // "bar_untilsat",
  // "cavlc_untilsat",
  "comparator_32bit_signed_lt_untilsat",
  "comparator_32bit_signed_lteq_untilsat",
  "comparator_32bit_unsigned_lt_untilsat",
  "comparator_32bit_unsigned_lteq_untilsat",
  // "ctrl_untilsat",
  // "dec_untilsat",
  // "div_untilsat",
  // "i2c_untilsat",
  // "int2float_untilsat",
  // "log2_untilsat",
  // "max_untilsat",
  "md5_untilsat",
  // "mem_ctrl_untilsat",
  "mult_32x32_untilsat",
  // "multiplier_untilsat",
  // "priority_untilsat",
  // "router_untilsat",
  "sha-1_untilsat",
  "sha-256_untilsat",
  // "sin_untilsat",
  // "sqrt_untilsat",
  // "square_untilsat",
  // "voter_untilsat"
};

std::vector<std::string> epfl_benchmarks( uint32_t selection = all )
{
  std::vector<std::string> result;
  for ( uint32_t i = 0u; i < 20u; ++i )
  {
    if ( ( selection >> i ) & 1 )
    {
      result.push_back( epfl_benchmark_names[i] );
    }
  }
  return result;
}

std::vector<std::string> crypto_benchmarks()
{
  std::vector<std::string> result;
  for ( uint32_t i = 0u; i < 15u; ++i )
  {
    result.push_back( crypto_benchmark_names[i] );
  }
  return result;
}

std::string benchmark_path( std::string const& benchmark_name, std::string const& path_type = "", std::string const& filetype = "v" ) 
{
#ifndef EXPERIMENTS_PATH
    return fmt::format( "{}.v", benchmark_name );
#else
  return fmt::format( "{}crypto_benchmarks{}/{}.{}", EXPERIMENTS_PATH, path_type, benchmark_name, filetype );
#endif
}

template<class Ntk>
bool abc_cec( Ntk const& ntk, std::string const& benchmark, std::string const& path_type = "", std::string const& filetype = "v" )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/test.bench\"", benchmark_path( benchmark, path_type, filetype ) );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  return result.size() >= 23 && result.substr( 0u, 23u ) == "Networks are equivalent";
}
template <class Ntk>
float abc_techmap( Ntk const& ntk, std::string const& genlib_path )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string const command = fmt::format( "abc -q \"read /tmp/test.bench; read_genlib {} ;map; print_gates\"", genlib_path );

  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std::cout << "results for mapping ============" <<  std::endl << result << std::endl;

  std::string const total_str = result.substr( result.find( "TOTAL " ) + 1 );
  uint32_t const sp = total_str.find( "Area" );
  uint32_t const lp = total_str.find( "100 \%" );
  std::string const str1 = total_str.substr( ( sp + 6 ), ( lp - sp - 6 ) ); // 6 as to ignore "="

  return std::stof( str1 );
}

template <class Ntk>
float abc_map_dc2 (Ntk const& ntk, std::string const& genlib_path )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"read /tmp/test.bench; read_genlib {} ;strash; dc2; map; print_gates\"", genlib_path);
  
  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std :: cout << "results for mapping ============" <<  std::endl << result << std::endl;

  std::string total_str = result.substr ( result.find( "TOTAL " ) + 1);
  uint32_t sp = total_str.find( "Area" );
  uint32_t lp = total_str.find( "100 \%" );
  //std::cout << " sp "<< sp << " lp " << lp << std::endl;
  std::string str1 = total_str.substr ( ( sp + 6 ), ( lp - sp - 6 ) ); // 6 as to ignore "=" 
  //std::cout << " value of str1 = " << str1 << std::endl;

  return std::stof( str1 );
}

template <class Ntk>
float abc_map_compress2rs (Ntk const& ntk, std::string const& genlib_path )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"read /tmp/test.bench; read_genlib {} ;strash; compress2rs; map; print_gates\"", genlib_path);

  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std :: cout << "results for mapping ============" <<  std::endl << result << std::endl;

  std::string total_str = result.substr ( result.find( "TOTAL " ) + 1);
  uint32_t sp = total_str.find( "Area" );
  uint32_t lp = total_str.find( "100 \%" );
  //std::cout << " sp "<< sp << " lp " << lp << std::endl;
  std::string str1 = total_str.substr ( ( sp + 6 ), ( lp - sp - 6 ) ); // 6 as to ignore "=" 
  //std::cout << " value of str1 = " << str1 << std::endl;

  return std::stof( str1 );
}

template <class Ntk>
float abc_map_dch (Ntk const& ntk, std::string const& genlib_path )
{ 
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"read /tmp/test.bench; read_genlib {} ;strash; dch; map; print_gates\"", genlib_path);
  
  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  { 
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  { 
    result += buffer.data();
  }
  
  std :: cout << "results for mapping ============" <<  std::endl << result << std::endl;
  
  std::string total_str = result.substr ( result.find( "TOTAL " ) + 1);
  uint32_t sp = total_str.find( "Area" );
  uint32_t lp = total_str.find( "100 \%" );
  //std::cout << " sp "<< sp << " lp " << lp << std::endl; 
  std::string str1 = total_str.substr ( ( sp + 6 ), ( lp - sp - 6 ) ); // 6 as to ignore "=" 
  //std::cout << " value of str1 = " << str1 << std::endl;
  
  return std::stof( str1 );
}

template <class Ntk>
float abc_map_rsrw (Ntk const& ntk, std::string const& genlib_path )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"read /tmp/test.bench; read_genlib {} ;strash; rw; rs; map; print_gates\"", genlib_path);

  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std :: cout << "results for mapping ============" <<  std::endl << result << std::endl;

  std::string total_str = result.substr ( result.find( "TOTAL " ) + 1);
  uint32_t sp = total_str.find( "Area" );
  uint32_t lp = total_str.find( "100 \%" );
  //std::cout << " sp "<< sp << " lp " << lp << std::endl;
  std::string str1 = total_str.substr ( ( sp + 6 ), ( lp - sp - 6 ) ); // 6 as to ignore "=" 
  //std::cout << " value of str1 = " << str1 << std::endl;

  return std::stof( str1 );
}


void abc_lut_reader_if( std::string const& benchmark )
{
  std::string command = fmt::format( "abc -q \"read {}; if -K 4; print_stats; write_bench {}\"", benchmark_path( benchmark ), benchmark_path( benchmark, "_if_bench", "bench") );

  std::array<char, 1024> buffer;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  std::string result;

  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }
  std::cout << "result LUT if-mapped ===============" << std::endl <<  std::endl;
  std::cout << result << std::endl;
}

void abc_lut_reader_mf( std::string const& benchmark )
{
  std::string command = fmt::format( "abc -q \"read {};&get; &mf -K 4;&put print_stats; write_bench {}\"", benchmark_path( benchmark ), benchmark_path( benchmark, "_mf_bench", "bench") );

  std::array<char, 1024> buffer;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  std::string result;

  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }
  std::cout << "result LUT mf-mapped ===============" << std::endl <<  std::endl;
  std::cout << result << std::endl;
}

struct lut_info
{
    uint32_t depth;
    uint32_t size;
};

template <class Ntk>
lut_info abc_lut_mapper_if( Ntk const& ntk )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"read /tmp/test.bench; if -K 6 ; print_stats\"" );

  lut_info ldata;
  std::array<char, 1024> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  std :: cout << "results for mapping ============" <<  std::endl << result << std::endl;

  std::string node_str = result.substr ( result.find( "nd =" ) + 5, 5);
  std::string lev_str = result.substr ( result.find( "lev" ) + 6 ); 
  std::cout << "node_str ---------> " << node_str << std::endl;
  std::cout << "lev_str ------>  " << lev_str << std::endl;
  ldata.depth = std::stoi( lev_str );
  ldata.size = std::stoi( node_str );
  //uint32_t sp = total_str.find( "Area" );
  //uint32_t lp = total_str.find( "100 \%" );
  //std::string str1 = total_str.substr ( ( sp + 6 ), ( lp - sp - 6 ) );
  return ldata;

}

} // namespace crypto_experiments
