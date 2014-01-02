//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2013 RIKEN
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Spatiocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Spatiocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Spatiocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
//

#ifndef __Common_hpp
#define __Common_hpp

#include <iostream>
#include <vector>
#include <bitset>
#include <climits>
#include <emmintrin.h> //SSE2 intrinsics
#include <immintrin.h> //AVX2 intrinsics

class Compartment;
class Diffuser;
class Model;
class Species;
class Stepper;
class VisualLogger;

#define TYPE unsigned
#define WORD (sizeof(TYPE)*8)
#define ADJS 12

struct Vector
{
  Vector(const double a=0, const double b=0, const double c=0):
    x(a),
    y(b),
    z(c) {}
  double x;
  double y;
  double z;
};

typedef union
{
  __m128i m128i[2];
  __m256i m256i;
  uint8_t uint8[32];
  uint16_t uint16[16];
  uint32_t uint32[8];
  int8_t int8[32];
  int16_t int16[16];
  int32_t int32[8];
} union256;

typedef union
{
  __m128i m128i;
  uint8_t uint8[16];
  uint16_t uint16[8];
  uint32_t uint32[4];
  int8_t int8[16];
  int16_t int16[8];
  int32_t int32[4];
} union128;

template<typename T>
void cout_binary(const T& a)
{
  const char* beg(reinterpret_cast<const char*>(&a));
  const char* end(beg + sizeof(a));
  while(1) {
    std::cout << std::bitset<CHAR_BIT>(*--end) << ' ';
    if(end == beg) {
      break;
    }
  }
  std::cout << std::endl;
}

#endif /* __Common_hpp */
