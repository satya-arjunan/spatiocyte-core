/*****************************   randomc.h   **********************************
* Author:        Agner Fog
* Date created:  1997
* Last modified: 2008-11-16
* Project:       randomc.h
* Source URL:    www.agner.org/random
*
* Description:
* This header file contains class declarations and other definitions for the 
* randomc class library of uniform random number generators in C++ language.
*
* Overview of classes:
* ====================
*
* class CRandomMersenne:
* Random number generator of type Mersenne twister.
* Source file mersenne.cpp
*
* class CRandomMother:
* Random number generator of type Mother-of-All (Multiply with carry).
* Source file mother.cpp
*
* class CRandomSFMT:
* Random number generator of type SIMD-oriented Fast Mersenne Twister.
* The class definition is not included here because it is not
* portable to all platforms. See sfmt.h and sfmt.cpp for details.
*
* Member functions (methods):
* ===========================
*
* All these classes have identical member functions:
*
* Constructor(int seed):
* The seed can be any integer. The time may be used as seed.
* Executing a program twice with the same seed will give the same sequence 
* of random numbers. A different seed will give a different sequence.
*
* void RandomInit(int seed);
* Re-initializes the random number generator with a new seed.
*
* void RandomInitByArray(int const seeds[], int NumSeeds);
* In CRandomMersenne and CRandomSFMT only: Use this function if you want 
* to initialize with a seed with more than 32 bits. All bits in the seeds[]
* array will influence the sequence of random numbers generated. NumSeeds 
* is the number of entries in the seeds[] array.
*
* double Random();
* Gives a floating point random number in the interval 0 <= x < 1.
* The resolution is 32 bits in CRandomMother and CRandomMersenne, and
* 52 bits in CRandomSFMT.
*
* int IRandom(int min, int max);
* Gives an integer random number in the interval min <= x <= max.
* (max-min < MAXINT).
* The precision is 2^-32 (defined as the difference in frequency between 
* possible output values). The frequencies are exact if max-min+1 is a
* power of 2.
*
* int IRandomX(int min, int max);
* Same as IRandom, but exact. In CRandomMersenne and CRandomSFMT only.
* The frequencies of all output values are exactly the same for an 
* infinitely long sequence. (Only relevant for extremely long sequences).
*
* uint32_t BRandom();
* Gives 32 random bits. 
*
*
* Example:
* ========
* The file EX-RAN.CPP contains an example of how to generate random numbers.
*
*
* Library version:
* ================
* Optimized versions of these random number generators are provided as function
* libraries in randoma.zip. These function libraries are coded in assembly
* language and support only x86 platforms, including 32-bit and 64-bit
* Windows, Linux, BSD, Mac OS-X (Intel based). Use randoma.h from randoma.zip
*
*
* Non-uniform random number generators:
* =====================================
* Random number generators with various non-uniform distributions are 
* available in stocc.zip (www.agner.org/random).
*
*
* Further documentation:
* ======================
* The file ran-instructions.pdf contains further documentation and 
* instructions for these random number generators.
*
* Copyright 1997-2008 by Agner Fog. 
* GNU General Public License http://www.gnu.org/licenses/gpl.html
*******************************************************************************/

#ifndef __Random_hpp
#define __Random_hpp

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
  // Compilers supporting C99 or C++0x have stdint.h defining these integer types
  #include <stdint.h>
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS) 
  // 16 bit systems use long int for 32 bit integer.
  typedef   signed long int int32_t;
  typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
  // Older Microsoft compilers have their own definition
  typedef   signed __int32  int32_t;
  typedef unsigned __int32 uint32_t;
  typedef   signed __int64  int64_t;
  typedef unsigned __int64 uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#else
  // This works with most compilers
  typedef signed int          int32_t;
  typedef unsigned int       uint32_t;
  typedef long long           int64_t;
  typedef unsigned long long uint64_t;
  #define INT64_SUPPORTED // Remove this if the compiler doesn't support 64-bit integers
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(const char *ErrorText);// System-specific error reporting (userintf.cpp)

#if defined(__cplusplus)               // class definitions only in C++
/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
// Choose which version of Mersenne Twister you want:
#if 0 
// Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else    
// or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif

public:
   CRandomMersenne(int seed) {         // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(int seed);          // Re-seed
   void RandomInitByArray(int const seeds[], int NumSeeds); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32_t BRandom();                 // Output random bits
private:
   void Init0(int seed);               // Basic initialization procedure
   uint32_t mt[MERS_N];                // State vector
   int mti;                            // Index into mt
   uint32_t LastInterval;              // Last interval length for IRandomX
   uint32_t RLimit;                    // Rejection limit used by IRandomX
};    


class CRandomMother {                  // Encapsulate random number generator
public:
   void RandomInit(int seed);          // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32_t BRandom();                 // Output random bits
   CRandomMother(int seed) {           // Constructor
      RandomInit(seed);}
protected:
   uint32_t x[5];                      // History buffer
};

#endif // __cplusplus


#include <emmintrin.h>                 // Define SSE2 intrinsics

// Choose one of the possible Mersenne exponents.
// Higher values give longer cycle length and use more memory:
//#define MEXP   607
//#define MEXP  1279
//#define MEXP  2281
//#define MEXP  4253
  #define MEXP 11213
//#define MEXP 19937
//#define MEXP 44497

// Define constants for the selected Mersenne exponent:
#if MEXP == 44497
#define SFMT_N    348                  // Size of state vector
#define SFMT_M    330                  // Position of intermediate feedback
#define SFMT_SL1    5                  // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	  3                  // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1    9                  // Right shift of W[M], 32-bit words
#define SFMT_SR2	  3                  // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xeffffffb,0xdfbebfff,0xbfbf7bef,0x9ffd7bff // AND mask
#define SFMT_PARITY 1,0,0xa3ac4000,0xecc1327a   // Period certification vector

#elif MEXP == 19937
#define SFMT_N    156                  // Size of state vector
#define SFMT_M    122                  // Position of intermediate feedback
#define SFMT_SL1   18                  // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	  1                  // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1   11                  // Right shift of W[M], 32-bit words
#define SFMT_SR2	  1                  // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xdfffffef,0xddfecb7f,0xbffaffff,0xbffffff6 // AND mask
#define SFMT_PARITY 1,0,0,0x13c9e684   // Period certification vector

#elif MEXP == 11213
#define SFMT_N    88                   // Size of state vector
#define SFMT_M    68                   // Position of intermediate feedback
#define SFMT_SL1	14                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 7                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 3                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0xeffff7fb,0xffffffef,0xdfdfbfff,0x7fffdbfd // AND mask
#define SFMT_PARITY 1,0,0xe8148000,0xd0c7afa3 // Period certification vector

#elif MEXP == 4253
#define SFMT_N    34                   // Size of state vector
#define SFMT_M    17                   // Position of intermediate feedback
#define SFMT_SL1	20                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 1                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 7                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0x9f7bffff, 0x9fffff5f, 0x3efffffb, 0xfffff7bb // AND mask
#define SFMT_PARITY 0xa8000001, 0xaf5390a3, 0xb740b3f8, 0x6c11486d // Period certification vector

#elif MEXP == 2281
#define SFMT_N    18                   // Size of state vector
#define SFMT_M    12                   // Position of intermediate feedback
#define SFMT_SL1	19                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 1                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 5                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	 0xbff7ffbf, 0xfdfffffe, 0xf7ffef7f, 0xf2f7cbbf // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x41dfa600  // Period certification vector

#elif MEXP == 1279
#define SFMT_N    10                   // Size of state vector
#define SFMT_M     7                   // Position of intermediate feedback
#define SFMT_SL1	14                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	 5                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 1                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xf7fefffd, 0x7fefcfff, 0xaff3ef3f, 0xb5ffff7f  // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x20000000  // Period certification vector

#elif MEXP == 607
#define SFMT_N     5                   // Size of state vector
#define SFMT_M     2                   // Position of intermediate feedback
#define SFMT_SL1	15                   // Left shift of W[N-1], 32-bit words
#define SFMT_SL2	 3                   // Left shift of W[0], *8, 128-bit words
#define SFMT_SR1	13                   // Right shift of W[M], 32-bit words
#define SFMT_SR2	 3                   // Right shift of W[N-2], *8, 128-bit words
#define SFMT_MASK	  0xfdff37ff, 0xef7f3f7d, 0xff777b7d, 0x7ff7fb2f  // AND mask
#define SFMT_PARITY 0x00000001, 0x00000000, 0x00000000, 0x5986f054  // Period certification vector
#endif

// Class for SFMT generator with or without Mother-Of-All generator
class CRandomSFMT {                              // Encapsulate random number generator
public:
   CRandomSFMT(int seed, int IncludeMother = 0) {// Constructor
      UseMother = IncludeMother; 
      LastInterval = 0;
      RandomInit(seed);}
   void RandomInit(int seed);                    // Re-seed
   void RandomInitByArray(int const seeds[], int NumSeeds); // Seed by more than 32 bits
   int  IRandom  (int min, int max);             // Output random integer
   int  IRandomX (int min, int max);             // Output random integer, exact
   double Random();                              // Output random floating point number
   uint32_t BRandom();                           // Output random bits
private:
   void Init2();                                 // Various initializations and period certification
   void Generate();                              // Fill state array with new random numbers
   uint32_t MotherBits();                        // Get random bits from Mother-Of-All generator
   uint32_t ix;                                  // Index into state array
   uint32_t LastInterval;                        // Last interval length for IRandom
   uint32_t RLimit;                              // Rejection limit used by IRandom
   uint32_t UseMother;                           // Combine with Mother-Of-All generator
   __m128i  mask;                                // AND mask
   __m128i  state[SFMT_N];                       // State vector for SFMT generator
   uint32_t MotherState[5];                      // State vector for Mother-Of-All generator
};

// Class for SFMT generator without Mother-Of-All generator
// Derived from CRandomSFMT
class CRandomSFMT0 : public CRandomSFMT {
public:
   CRandomSFMT0(int seed) : CRandomSFMT(seed,0) {}
};

// Class for SFMT generator combined with Mother-Of-All generator
// Derived from CRandomSFMT
class CRandomSFMT1 : public CRandomSFMT {
public:
   CRandomSFMT1(int seed) : CRandomSFMT(seed,1) {}
};

#endif /* __Random_hpp */ 
