//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2014 RIKEN
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
// based on the SFMT code by Agner Fog
//

#include <Random.hpp>

void Random::RandomInit(int seed) {
   // Re-seed
   uint32_t i;                         // Loop counter
   uint32_t y = seed;                  // Temporary
   uint32_t statesize = SFMT_N*4;      // Size of state vector

   // Fill state vector with random numbers from seed
   ((uint32_t*)state)[0] = y;
   const uint32_t factor = 1812433253U;// Multiplication factor

   for (i = 1; i < statesize; i++) {
      y = factor * (y ^ (y >> 30)) + i;
      ((uint32_t*)state)[i] = y;
   }

   // Further initialization and period certification
   Init2();
}

// Functions used by Random::RandomInitByArray
static uint32_t func1(uint32_t x) {
    return (x ^ (x >> 27)) * 1664525U;
}

static uint32_t func2(uint32_t x) {
    return (x ^ (x >> 27)) * 1566083941U;
}

void Random::RandomInitByArray(int const seeds[], int NumSeeds) {
   // Seed by more than 32 bits
   uint32_t i, j, count, r, lag;

   if (NumSeeds < 0) NumSeeds = 0;

   const uint32_t size = SFMT_N*4; // number of 32-bit integers in state

   // Typecast state to uint32_t *
   uint32_t * sta = (uint32_t*)state;

   if (size >= 623) {
      lag = 11;} 
   else if (size >= 68) {
      lag = 7;}
   else if (size >= 39) {
      lag = 5;}
   else {
      lag = 3;
   }
   const uint32_t mid = (size - lag) / 2;

   if ((uint32_t)NumSeeds + 1 > size) {
      count = (uint32_t)NumSeeds;
   }
   else {
      count = size - 1;
   }
#if 0
   // Original code. Argument to func1 is constant!
   for (i = 0; i < size; i++) sta[i] = 0x8B8B8B8B;
   r = func1(sta[0] ^ sta[mid] ^ sta[size - 1]);
   sta[mid] += r;
   r += NumSeeds;
   sta[mid + lag] += r;
   sta[0] = r;
#else
   // 1. loop: Fill state vector with random numbers from NumSeeds
   const uint32_t factor = 1812433253U;// Multiplication factor
   r = (uint32_t)NumSeeds;
   for (i = 0; i < SFMT_N*4; i++) {
      r = factor * (r ^ (r >> 30)) + i;
      sta[i] = r;
   }

#endif

   // 2. loop: Fill state vector with random numbers from seeds[]
   for (i = 1, j = 0; j < count; j++) {
      r = func1(sta[i] ^ sta[(i + mid) % size] ^ sta[(i + size - 1) % size]);
      sta[(i + mid) % size] += r;
      if (j < (uint32_t)NumSeeds) r += (uint32_t)seeds[j];
      r += i;
      sta[(i + mid + lag) % size] += r;
      sta[i] = r;
      i = (i + 1) % size;
   }

   // 3. loop: Randomize some more
   for (j = 0; j < size; j++) {
      r = func2(sta[i] + sta[(i + mid) % size] + sta[(i + size - 1) % size]);
      sta[(i + mid) % size] ^= r;
      r -= i;
      sta[(i + mid + lag) % size] ^= r;
      sta[i] = r;
      i = (i + 1) % size;
   }
   
   // Further initialization and period certification
   Init2();
}


void Random::Init2() {
   // Various initializations and period certification
   uint32_t i, j, temp;

   // Initialize mask
   static const uint32_t maskinit[4] = {SFMT_MASK};
   mask = _mm_loadu_si128((__m128i*)maskinit);

   // Period certification
   // Define period certification vector
   static const uint32_t parityvec[4] = {SFMT_PARITY};

   // Check if parityvec & state[0] has odd parity
   temp = 0;
   for (i = 0; i < 4; i++) {
      temp ^= parityvec[i] & ((uint32_t*)state)[i];
   }
   for (i = 16; i > 0; i >>= 1) temp ^= temp >> i;
   if (!(temp & 1)) {
      // parity is even. Certification failed
      // Find a nonzero bit in period certification vector
      for (i = 0; i < 4; i++) {
         if (parityvec[i]) {
            for (j = 1; j; j <<= 1) {
               if (parityvec[i] & j) {
                  // Flip the corresponding bit in state[0] to change parity
                  ((uint32_t*)state)[i] ^= j;
                  // Done. Exit i and j loops
                  i = 5;  break;
               }
            }
         }
      }
   }
   // Generate first random numbers and set ix = 0
   Generate();
}


// Subfunction for the sfmt algorithm
static inline __m128i sfmt_recursion(__m128i const &a, __m128i const &b, 
__m128i const &c, __m128i const &d, __m128i const &mask) {
    __m128i a1, b1, c1, d1, z1, z2;
    b1 = _mm_srli_epi32(b, SFMT_SR1);
    a1 = _mm_slli_si128(a, SFMT_SL2);
    c1 = _mm_srli_si128(c, SFMT_SR2);
    d1 = _mm_slli_epi32(d, SFMT_SL1);
    b1 = _mm_and_si128(b1, mask);
    z1 = _mm_xor_si128(a, a1);
    z2 = _mm_xor_si128(b1, d1);
    z1 = _mm_xor_si128(z1, c1);
    z2 = _mm_xor_si128(z1, z2);
    return z2;
}

void Random::Generate() {
   // Fill state array with new random numbers
   int i;
   __m128i r, r1, r2;

   r1 = state[SFMT_N - 2];
   r2 = state[SFMT_N - 1];
   for (i = 0; i < SFMT_N - SFMT_M; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   for (; i < SFMT_N; i++) {
      r = sfmt_recursion(state[i], state[i + SFMT_M - SFMT_N], r1, r2, mask);
      state[i] = r;
      r1 = r2;
      r2 = r;
   }
   ix = 0;
}

int  Random::IRan(int min, int max) {
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32

   interval = (uint32_t)(max - min + 1);
   longran  = (uint64_t)BRan() * interval;
   iran = (uint32_t)(longran >> 32);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}

uint32_t Random::RanUint32_12() {
   return (uint32_t)(((uint64_t)BRan()*12) >> 32);
}

uint32_t Random::BRan() {
   // Output 32 random bits
   uint32_t y;

   if (ix >= SFMT_N*4) {
      Generate();
   }
   y = ((uint32_t*)state)[ix++];
   return y;
}

/*
union256i_d Random::BRan8() {
   // Output 256 random bits
   union256i_d y;
   if (ix >= SFMT_N/2) {
      Generate();
   }
   y.x = ((__m256i*)state)[ix++];
   return y;
}

union256i_d Random::IRan8_12() {
   uint32_t interval;
   uint64_t longran;
   uint32_t iran;
   interval = 13;
   longran  = (uint64_t)BRan() * interval;
   iran = (uint32_t)(longran >> 32);
   return (int32_t)iran;
}
*/

int  Random::IRanX (int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Each output value has exactly the same probability.
   // This is obtained by rejecting certain bit values so that the number
   // of possible bit values is divisible by the interval length
   if (max <= min) {
      if (max == min) {
         return min;                   // max == min. Only one possible value
      }
      else {
         return 0x80000000;            // max < min. Error output
      }
   }
   // Assume 64 bit integers supported. Use multiply and shift method
   uint32_t interval;                  // Length of interval
   uint64_t longran;                   // Random bits * interval
   uint32_t iran;                      // Longran / 2^32
   uint32_t remainder;                 // Longran % 2^32

   interval = (uint32_t)(max - min + 1);
   if (interval != LastInterval) {
      // Interval length has changed. Must calculate rejection limit
      // Reject when remainder = 2^32 / interval * interval
      // RLimit will be 0 if interval is a power of 2. No rejection then.
      RLimit = (uint32_t)(((uint64_t)1 << 32) / interval) * interval - 1;
      LastInterval = interval;
   }
   do { // Rejection loop
      longran  = (uint64_t)BRan() * interval;
      iran = (uint32_t)(longran >> 32);
      remainder = (uint32_t)longran;
   } while (remainder > RLimit);
   // Convert back to signed and return result
   return (int32_t)iran + min;
}

double Random::Ran() {
   // Output random floating point number
   if (ix >= SFMT_N*4-1) {
      // Make sure we have at least two 32-bit numbers
      Generate();
   }
   uint64_t r = *(uint64_t*)((uint32_t*)state+ix);
   ix += 2;
   // 53 bits resolution:
   // return (int64_t)(r >> 11) * (1./(67108864.0*134217728.0)); 
   // (r >> 11)*2^(-53)
   // 52 bits resolution for compatibility with assembly version:
   // (r >> 12)*2^(-52)
   return (int64_t)(r >> 12) * (1./(67108864.0*67108864.0));
}
