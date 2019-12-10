// //////////////////////////////////////////////////////////
// 
// Copyright (c) 2014,2015 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//
#pragma once
//#include "hash.h"
#include <string>
// define fixed size integer types
#ifdef _MSC_VER
// Windows
typedef unsigned __int8  uint8_t;
typedef unsigned __int64 uint64_t;
#else
// GCC
#include <stdint.h>
#endif
/// compute SHA3 hash
/** Usage:
    SHA3 sha3;
    std::string myHash  = sha3("Hello World");     // std::string
    std::string myHash2 = sha3("How are you", 11); // arbitrary data, 11 bytes
    // or in a streaming fashion:
    SHA3 sha3;
    while (more data available)
      sha3.add(pointer to fresh data, number of new bytes);
    std::string myHash3 = sha3.getHash();
  */
class SHA3 //: public Hash
{
public:
  /// algorithm variants
  enum Bits { Bits224 = 224, Bits256 = 256, Bits384 = 384, Bits512 = 512 };
  /// same as reset()
  explicit SHA3(Bits bits = Bits256);
  /// compute hash of a memory block
  std::string operator()(const void* data, size_t numBytes);
  /// compute hash of a string, excluding final zero
  std::string operator()(const std::string& text);
  /// add arbitrary number of bytes
  void add(const void* data, size_t numBytes);
  /// return latest hash as hex characters
  std::string getHash();
  /// restart
  void reset();
private:
  /// process a full block
  void processBlock(const void* data);
  /// process everything left in the internal buffer
  void processBuffer();
  /// 1600 bits, stored as 25x64 bit, BlockSize is no more than 1152 bits (Keccak224)
  enum { StateSize    = 1600 / (8 * 8),
         MaxBlockSize =  200 - 2 * (224 / 8) };
  /// hash
  uint64_t m_hash[StateSize];
  /// size of processed data in bytes
  uint64_t m_numBytes;
  /// block size (less or equal to MaxBlockSize)
  size_t   m_blockSize;
  /// valid bytes in m_buffer
  size_t   m_bufferSize;
  /// bytes not processed yet
  uint8_t  m_buffer[MaxBlockSize];
  /// variant
  Bits     m_bits;
};


// //////////////////////////////////////////////////////////
// sha3.cpp
// Copyright (c) 2014,2015 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//
//#include "
"
// big endian architectures need #define __BYTE_ORDER __BIG_ENDIAN
#ifndef _MSC_VER
#include <endian.h>
#endif
/// same as reset()
SHA3::SHA3(Bits bits)
: m_blockSize(200 - 2 * (bits / 8)),
  m_bits(bits)
{
  reset();
}
/// restart
void SHA3::reset()
{
  for (size_t i = 0; i < StateSize; i++)
    m_hash[i] = 0;
  m_numBytes   = 0;
  m_bufferSize = 0;
}
/// constants and local helper functions
namespace
{
  const unsigned int Rounds = 24;
  const uint64_t XorMasks[Rounds] =
  {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
  };
  /// rotate left and wrap around to the right
  inline uint64_t rotateLeft(uint64_t x, uint8_t numBits)
  {
    return (x << numBits) | (x >> (64 - numBits));
  }
  /// convert litte vs big endian
  inline uint64_t swap(uint64_t x)
  {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(x);
#endif
#ifdef _MSC_VER
    return _byteswap_uint64(x);
#endif
    return  (x >> 56) |
           ((x >> 40) & 0x000000000000FF00ULL) |
           ((x >> 24) & 0x0000000000FF0000ULL) |
           ((x >>  8) & 0x00000000FF000000ULL) |
           ((x <<  8) & 0x000000FF00000000ULL) |
           ((x << 24) & 0x0000FF0000000000ULL) |
           ((x << 40) & 0x00FF000000000000ULL) |
            (x << 56);
  }
  /// return x % 5 for 0 <= x <= 9
  unsigned int mod5(unsigned int x)
  {
    if (x < 5)
      return x;
    return x - 5;
  }
}
/// process a full block
void SHA3::processBlock(const void* data)
{
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && (__BYTE_ORDER == __BIG_ENDIAN)
#define LITTLEENDIAN(x) swap(x)
#else
#define LITTLEENDIAN(x) (x)
#endif
  const uint64_t* data64 = (const uint64_t*) data;
  // mix data into state
  for (unsigned int i = 0; i < m_blockSize / 8; i++)
    m_hash[i] ^= LITTLEENDIAN(data64[i]);
  // re-compute state
  for (unsigned int round = 0; round < Rounds; round++)
  {
    // Theta
    uint64_t coefficients[5];
    for (unsigned int i = 0; i < 5; i++)
      coefficients[i] = m_hash[i] ^ m_hash[i + 5] ^ m_hash[i + 10] ^ m_hash[i + 15] ^ m_hash[i + 20];
    for (unsigned int i = 0; i < 5; i++)
    {
      uint64_t one = coefficients[mod5(i + 4)] ^ rotateLeft(coefficients[mod5(i + 1)], 1);
      m_hash[i     ] ^= one;
      m_hash[i +  5] ^= one;
      m_hash[i + 10] ^= one;
      m_hash[i + 15] ^= one;
      m_hash[i + 20] ^= one;
    }
    // temporary
    uint64_t one;
    // Rho Pi
    uint64_t last = m_hash[1];
    one = m_hash[10]; m_hash[10] = rotateLeft(last,  1); last = one;
    one = m_hash[ 7]; m_hash[ 7] = rotateLeft(last,  3); last = one;
    one = m_hash[11]; m_hash[11] = rotateLeft(last,  6); last = one;
    one = m_hash[17]; m_hash[17] = rotateLeft(last, 10); last = one;
    one = m_hash[18]; m_hash[18] = rotateLeft(last, 15); last = one;
    one = m_hash[ 3]; m_hash[ 3] = rotateLeft(last, 21); last = one;
    one = m_hash[ 5]; m_hash[ 5] = rotateLeft(last, 28); last = one;
    one = m_hash[16]; m_hash[16] = rotateLeft(last, 36); last = one;
    one = m_hash[ 8]; m_hash[ 8] = rotateLeft(last, 45); last = one;
    one = m_hash[21]; m_hash[21] = rotateLeft(last, 55); last = one;
    one = m_hash[24]; m_hash[24] = rotateLeft(last,  2); last = one;
    one = m_hash[ 4]; m_hash[ 4] = rotateLeft(last, 14); last = one;
    one = m_hash[15]; m_hash[15] = rotateLeft(last, 27); last = one;
    one = m_hash[23]; m_hash[23] = rotateLeft(last, 41); last = one;
    one = m_hash[19]; m_hash[19] = rotateLeft(last, 56); last = one;
    one = m_hash[13]; m_hash[13] = rotateLeft(last,  8); last = one;
    one = m_hash[12]; m_hash[12] = rotateLeft(last, 25); last = one;
    one = m_hash[ 2]; m_hash[ 2] = rotateLeft(last, 43); last = one;
    one = m_hash[20]; m_hash[20] = rotateLeft(last, 62); last = one;
    one = m_hash[14]; m_hash[14] = rotateLeft(last, 18); last = one;
    one = m_hash[22]; m_hash[22] = rotateLeft(last, 39); last = one;
    one = m_hash[ 9]; m_hash[ 9] = rotateLeft(last, 61); last = one;
    one = m_hash[ 6]; m_hash[ 6] = rotateLeft(last, 20); last = one;
                      m_hash[ 1] = rotateLeft(last, 44);
    // Chi
    for (unsigned int j = 0; j < 25; j += 5)
    {
      // temporaries
      uint64_t one = m_hash[j];
      uint64_t two = m_hash[j + 1];
      m_hash[j]     ^= m_hash[j + 2] & ~two;
      m_hash[j + 1] ^= m_hash[j + 3] & ~m_hash[j + 2];
      m_hash[j + 2] ^= m_hash[j + 4] & ~m_hash[j + 3];
      m_hash[j + 3] ^=      one      & ~m_hash[j + 4];
      m_hash[j + 4] ^=      two      & ~one;
    }
    // Iota
    m_hash[0] ^= XorMasks[round];
  }
}
/// add arbitrary number of bytes
void SHA3::add(const void* data, size_t numBytes)
{
  const uint8_t* current = (const uint8_t*) data;
  // copy data to buffer
  if (m_bufferSize > 0)
  {
    while (numBytes > 0 && m_bufferSize < m_blockSize)
    {
      m_buffer[m_bufferSize++] = *current++;
      numBytes--;
    }
  }
  // full buffer
  if (m_bufferSize == m_blockSize)
  {
    processBlock((void*)m_buffer);
    m_numBytes  += m_blockSize;
    m_bufferSize = 0;
  }
  // no more data ?
  if (numBytes == 0)
    return;
  // process full blocks
  while (numBytes >= m_blockSize)
  {
    processBlock(current);
    current    += m_blockSize;
    m_numBytes += m_blockSize;
    numBytes   -= m_blockSize;
  }
  // keep remaining bytes in buffer
  while (numBytes > 0)
  {
    m_buffer[m_bufferSize++] = *current++;
    numBytes--;
  }
}
/// process everything left in the internal buffer
void SHA3::processBuffer()
{
  // add padding
  size_t offset = m_bufferSize;
  // add a "1" byte
  m_buffer[offset++] = 0x06;
  // fill with zeros
  while (offset < m_blockSize)
    m_buffer[offset++] = 0;
  // and add a single set bit
  m_buffer[offset - 1] |= 0x80;
  processBlock(m_buffer);
}
/// return latest hash as 16 hex characters
std::string SHA3::getHash()
{
  // process remaining bytes
  processBuffer();
  // convert hash to string
  static const char dec2hex[16 + 1] = "0123456789abcdef";
  // number of significant elements in hash (uint64_t)
  unsigned int hashLength = m_bits / 64;
  std::string result;
  result.reserve(m_bits / 4);
  for (unsigned int i = 0; i < hashLength; i++)
    for (unsigned int j = 0; j < 8; j++) // 64 bits => 8 bytes
    {
      // convert a byte to hex
      unsigned char oneByte = (unsigned char) (m_hash[i] >> (8 * j));
      result += dec2hex[oneByte >> 4];
      result += dec2hex[oneByte & 15];
    }
  // SHA3-224's last entry in m_hash provides only 32 bits instead of 64 bits
  unsigned int remainder = m_bits - hashLength * 64;
  unsigned int processed = 0;
  while (processed < remainder)
  {
    // convert a byte to hex
    unsigned char oneByte = (unsigned char) (m_hash[hashLength] >> processed);
    result += dec2hex[oneByte >> 4];
    result += dec2hex[oneByte & 15];
    processed += 8;
  }
  return result;
}
/// compute SHA3 of a memory block
std::string SHA3::operator()(const void* data, size_t numBytes)
{
  reset();
  add(data, numBytes);
  return getHash();
}
/// compute SHA3 of a string, excluding final zero
std::string SHA3::operator()(const std::string& text)
{
  reset();
  add(text.c_str(), text.size());
  return getHash();
}


// //////////////////////////////////////////////////////////
// digest.cpp
// Copyright (c) 2014,2015 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
//
// g++ -O3 sha3.cpp -o sha3
#include <iostream>
#include <fstream>
int main(int argc, char** argv)
{
  // syntax check
  if (argc < 2 || argc > 3)
  {
    std::cout << "./digest filename [--crc|--md5|--sha1|--sha256|--keccak|--sha3]" << std::endl;
    return 1;
  }
  // parameters
  std::string filename  = argv[1];
  std::string algorithm = argc == 3 ? argv[2] : "";
  bool computeSha3      = algorithm.empty() || algorithm == "--sha3";
  SHA3   digestSha3  (SHA3  ::Bits256);
  // each cycle processes about 1 MByte (divisible by 144 => improves Keccak/SHA3 performance)
  const size_t BufferSize = 144*7*1024;
  char* buffer = new char[BufferSize];
  // select input source: either file or standard-in
  std::ifstream file;
  std::istream* input = NULL;
  // accept std::cin, syntax will be: "./digest - --sha3 < data"
  if (filename == "-")
  {
    input = &std::cin;
  }
  else
  {
    // open file
    file.open(filename.c_str(), std::ios::in | std::ios::binary);
    if (!file)
    {
      std::cerr << "Can't open '" << filename << "'" << std::endl;
      return 2;
    }
    input = &file;
  }
  // process file
  while (*input)
  {
    input->read(buffer, BufferSize);
    std::size_t numBytesRead = size_t(input->gcount());
    if (computeSha3)
      digestSha3  .add(buffer, numBytesRead);
  }
  // clean up
  file.close();
  delete[] buffer;
  // show results
  if (computeSha3)
    std::cout << "SHA3/256:   " << digestSha3  .getHash() << std::endl;
  return 0;
}
