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
//

#include <math.h>
#include <Compartment.hpp>
#include <Species.hpp>
#include <Model.hpp>
#include <Optimization.hpp>


Compartment::Compartment(std::string name, const double len_x,
    const double len_y, const double len_z, Model& model,
    const unsigned num_box)
  : name_(name),
    model_(model),
    lattice_(NUM_VOXEL, Vector<unsigned>(NUM_COL, NUM_ROW, NUM_LAY),
        num_box),
    dimensions_(get_lattice().get_dimensions().x*2*VOXEL_RADIUS, 
        get_lattice().get_dimensions().y*2*VOXEL_RADIUS,
        get_lattice().get_dimensions().z*2*VOXEL_RADIUS),
    volume_species_("volume", 0, 0, model, *this, volume_species_, true),
    surface_species_("surface", 0, 0, model, *this, volume_species_, true) {}

void Compartment::initialize() {
  lattice_.initialize();
  set_offsets();
  set_volume_structure();
  set_surface_structure();
  umol_t multiplier_colrow, multiplier_row, nshift_colrow, nshift_row;
  set_const_division_param(NUM_COLROW, &multiplier_colrow, &nshift_colrow);
  multiplier_colrow_ = _mm256_set1_epi16(multiplier_colrow);
  nshift_colrow_ = _mm_set_epi64x((uint64_t)0, (uint64_t)nshift_colrow);
  set_const_division_param(NUM_ROW, &multiplier_row, &nshift_row);
  multiplier_row_ = _mm256_set1_epi16(multiplier_row);
  nshift_row_ = _mm_set_epi64x((uint64_t)0, (uint64_t)nshift_row);
}

umol_t Compartment::get_tar(const umol_t vdx, const unsigned nrand) const {
  const bool odd_lay((vdx/NUM_COLROW)&1);
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  //return vdx+offsets_[nrand+odd_lay*24+odd_col*12];
  //return vdx+offsets_[nrand+odd_lay+odd_col];
  int val(int(vdx)+offsets_[nrand+(24&(-odd_lay))+(12&(-odd_col))]);
  if(val < 0 || val > NUM_VOXEL) {
    return vdx;
  }
  return val;
  /*
  The registers:
    rax = 64 bit
      eax = lower 32 bits
    rbx = 64 bit
      ebx = lower 32 bits
    rcx = 64 bit
      ecx = lower 32 bits
    rdx = 64 bit
      edx = lower 32 bits
    rsi = 64 bit
      esi = lower 32 bits
    rsp = 64 bit
      esp = lower 32 bits
    r8 = 64 bit
      r8d = lower 32 bits
    r9 = 64 bit
      r9d = lower 32 bits
  
  Function parameters:
  System V AMD64 ABI[11] is followed on Solaris, GNU/Linux, FreeBSD, Mac OS X,
    and other UNIX-like or POSIX-compliant operating systems.
    The first six integer arguments (from the left) are passed in:
      RDI, RSI, RDX, RCX, R8, and R9, in that order.
      Additional integer arguments are passed on the stack.
      These registers, plus RAX, R10 and R11 are destroyed by function calls.
      For system calls, R10 is used instead of RCX.
    Integer return values are passed in RAX and RDX, in that order.
    Floating point is done using SSE registers, except for long double.
      Floating-point arguments are passed in XMM0 to XMM7. 
    Floating point return is XMM0 and XMM1. 
    Long double are passed on the stack, and returned in ST0 and ST1.
    On 64-bit Unix, long is 64 bits.
    Integer and SSE register arguments are counted separately, 
      so for the case of
        void foo(long a, double b, int c)
        a is passed in RDI, b in XMM0, and c in ESI.
    In C++ classes, the "this" pointer is passed as the first integer parameter
      the next three parameters are passed in the registers, while the rest
      are passed on the stack

  GNU assembler coding:
  mnemonic source, destination
  There are up to 4 parameters of an address operand that are presented in:
    displacement(base register, offset register, scalar multiplier)
    which is equivalent to:
      base register + displacement + offset register*scalar multiplier
  suffixes:
    b = byte (8 bit)
    s = short (16 bit integer) or single (32-bit floating point)
    w = word (16 bit)
    l = long (32 bit integer or 64-bit floating point)
    q = quad (64 bit)
    t = ten bytes (80-bit floating point)
  prefixes:
    when referencing registers, prefix it with "%"
    when using constant numbers, prefix it with "$"
  examples:
  movq [rbx], rax: moves 8 bytes beginning at rbx into rax
  movl -4(%ebp, %edx, 4), %eax: load *(ebp - 4 + (edx * 4)) into eax
  movl -4(%ebp), %eax: load a stack variable into eax
  movl (%ecx), %edx: copy the target of a pointer into a register
  leal 8(,%eax,4), %eax: multiply eax by 4 and add 8
  leal (%eax,%eax,2), %eax: multiply eax by 2 and add eax (i.e. multiply by 3)

  EVEN row numbers: 
  unsigned get_tar(const unsigned vdx, const unsigned nrand)
  const bool odd_lay((vdx/47066)&1);
  const bool odd_col((vdx%47066/202)&1);
  return vdx+offsets_[nrand+(24&(-odd_lay))+(12&(-odd_col))];

  Corresponding assembly:
    rdi: contains the "this" pointer
    esi: contains vdx
    edx: contains nrand
    eax: will contain the return value

	movl	%esi, %eax         : eax contains vdx
	movl	$747553905, %r8d   : r8d contains the magic number of 47066: 747553905
	movl	%edx, %r9d         : r9d contains nrand
	mull	%r8d               : mul is unsigned multiply (page 3-586 Vol2A)
                             edx:eax = eax*r8d, so edx:eax = vdx*747553905
	movl	%esi, %eax         : eax contains vdx
	shrl	$13, %edx          : unsigned divide source by 2^13 (page 4-333 Vol2B)
                             or shift logical right by 13 times
                             we shift by 13 times to get division result
                             edx = edx >> 13
                             edx = MostSignificant32bits(vdx*747553905) >> 13
                             edx = vdx/47066
	imull	$47066, %edx, %r8d : signed multiply (page 3-387 Vol2A)
                             result is truncated to the size of destination
                             register, r8d which is 32 bits
                             r8d = 47066*edx
                             r8d = 47066*(vdx/47066)
	movl	%edx, %ecx         : ecx contains vdx/47066
	movl	$680390859, %edx   : edx contains the magic number of 202: 680390859
	andl	$1, %ecx           : ecx = ecx&1
                             ecx = (vdx/47066)&1
                             ecx = odd_lay
	negl	%ecx               : ecx = -odd_lay
	subl	%r8d, %eax         : eax = eax-r8d
                             eax = vdx-47066*(vdx/47066)
                             eax = vdx%47066
	andl	$24, %ecx          : ecx = 24&ecx
                             ecx = 24&(-odd_lay)
	shrl	%eax               : eax = eax/2
                             eax = (vdx%47066)/2
	addl	%r9d, %ecx         : ecx = ecx+r9d
                             ecx = 24&(-odd_lay)+nrand
	mull	%edx               : edx:eax = eax*edx
                             edx:eax = [(vdx%47066)/2]*680390859
	movq	22056(%rdi), %rax  : rax = this+22056
                             rax = offset_
	sall	$27, %edx          : shift arithmetic left
                             multiply edx by 2 for 27 times
                             edx = 27 << edx 
                                 = edx*2^27
                                 = 27 << MS32{[(vdx%47066)/2]*680390859} 
	sarl	$31, %edx          : signed divide edx by 2 for 31 times
                             edx = edx >> 31
                                 = (27 << MS32{[(vdx%47066)/2]*680390859}) >> 31
                                 = -odd_col
	andl	$12, %edx          : edx = 12&edx
                                 = 12&(-odd_col)
	addl	%ecx, %edx         : edx = edx + ecx
                                 = 12&(-odd_col)+24&(-odd_lay)+nrand
	addl	(%rax,%rdx,4), %esi: esi = (rax+rdx*4)+esi
                               = offset_[12&(-odd_col)+24&(-odd_lay)+nrand]+vdx
	movl	%esi, %eax         : eax = esi
                             return esi
	ret
  */
  /*
  num_colrow = 1440

 	movzwl	%si, %r8d
	movl	%esi, %eax
	movq	22056(%rdi), %rdi  : rdi = this+22056
	imull	$11651, %r8d, %r8d
	shrl	$24, %r8d
	movl	%r8d, %ecx
	imulw	$1440, %r8w, %r8w :
	andl	$1, %ecx
	negl	%ecx
	subl	%r8d, %eax
	andl	$24, %ecx
	movzwl	%ax, %eax
	addl	%edx, %ecx
	movl	%esi, %edx
	imull	$58255, %eax, %eax
	sall	$10, %eax
	sarl	$31, %eax
	andl	$12, %eax
	addl	%ecx, %eax
	addw	(%rdi,%rax,4), %dx
	movl	%edx, %eax
	ret
  */
  //Move integer values from an aligned memory location (32 bytes total)
  //The integers could be uint8_t, umol_t or uint32_t
  //We first use uint32_t, so there are 8 vdx values:
  //vdx.m256i = _mm256_load_si256(mvdx);
}

/*AVX2 SIMD implementation of:
 *  const bool odd_lay((vdx/NUM_COLROW)&1);
 *  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
 *  return vdx+offsets_[nrand+(24&(-odd_lay))+(12&(-odd_col))];
 */
void Compartment::set_tars(const __m256i vdx, __m256i nrand,
    uint32_t* tars) const { 
  //[mull] multiply unsigned and store high 16 bit result
  //vdx*multiplier
  __m256i quotient_colrow(_mm256_mulhi_epu16(vdx, multiplier_colrow_));
  //[shrl] shift right logical (set the new bits on the left as 0)
  //vdx/num_colrow = vdx*multiplier/2^nshift_colrow_
  quotient_colrow = _mm256_srl_epi16(quotient_colrow, nshift_colrow_);
  __m256i odd_lay(_mm256_and_si256(quotient_colrow, _mm256_set1_epi16(1)));
  //[imull] signed multiply and store low 16 bit result
  //(vdx/num_colrow)*num_colrow
  quotient_colrow = _mm256_mullo_epi16(quotient_colrow,
                                       _mm256_set1_epi16(NUM_COLROW));
  //[subl] a-b 
  //remainder = vdx-(vdx/num_colrow)*num_colrow
  __m256i remainder_colrow(_mm256_sub_epi16(vdx, quotient_colrow));
  //remainder*multiplier
  __m256i quotient_row(_mm256_mulhi_epu16(remainder_colrow, multiplier_row_));
  //remainder*multiplier/2^nshift_row_
  quotient_row = _mm256_srl_epi16(quotient_row, nshift_row_);
  __m256i odd_col(_mm256_and_si256(quotient_row, _mm256_set1_epi16(1)));
  //Option 2 faster
  odd_lay = _mm256_sign_epi16(odd_lay, _mm256_set1_epi64x(-1));
  odd_lay = _mm256_and_si256(odd_lay, _mm256_set1_epi16(24));
  odd_col = _mm256_sign_epi16(odd_col, _mm256_set1_epi64x(-1));
  odd_col = _mm256_and_si256(odd_col, _mm256_set1_epi16(12));
  nrand = _mm256_add_epi16(nrand, odd_lay);
  nrand = _mm256_add_epi16(nrand, odd_col);
  //cast first 8 indices from uint16_t to uint32_t
  __m256i index(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(nrand)));
  //get the first 8 offsets
  index = _mm256_i32gather_epi32(offsets_, index, 4);
  //cast first 8 vdx from uint16_t to uint32_t and add with offset
  *(__m256i*)(tars) = 
    _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(vdx)), index);
  //cast second 8 indices from uint16_t to uint32_t
  index = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(nrand, 1));
  //get the second 8 offsets
  index =  _mm256_i32gather_epi32(offsets_, index, 4);
  //cast second 8 vdx from uint16_t to uint32_t and add with offset
  *(__m256i*)(&tars[8]) =
    _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_extractf128_si256(vdx, 1)),
                     index);
}

/*AVX2 SIMD implementation of:
 *  const bool odd_lay((vdx/NUM_COLROW)&1);
 *  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
 *  return vdx+offsets_[nrand+(24&(-odd_lay))+(12&(-odd_col))];
 */
__m256i Compartment::get_tars(const __m256i vdx, __m256i nrand) const { 
  /*
  union256 arand;
  arand.m256i = nrand;
  */

  //[mull] multiply unsigned and store high 16 bit result
  //vdx*multiplier
  __m256i quotient_colrow(_mm256_mulhi_epu16(vdx, multiplier_colrow_));
  //[shrl] shift right logical (set the new bits on the left as 0)
  //vdx/num_colrow = vdx*multiplier/2^nshift_colrow_
  quotient_colrow = _mm256_srl_epi16(quotient_colrow, nshift_colrow_);
  __m256i odd_lay(_mm256_and_si256(quotient_colrow, _mm256_set1_epi16(1)));
  //[imull] signed multiply and store low 16 bit result
  //(vdx/num_colrow)*num_colrow
  quotient_colrow = _mm256_mullo_epi16(quotient_colrow,
                                       _mm256_set1_epi16(NUM_COLROW));
  //[subl] a-b 
  //remainder = vdx-(vdx/num_colrow)*num_colrow
  __m256i remainder_colrow(_mm256_sub_epi16(vdx, quotient_colrow));
  //remainder*multiplier
  __m256i quotient_row(_mm256_mulhi_epu16(remainder_colrow, multiplier_row_));
  //remainder*multiplier/2^nshift_row_
  quotient_row = _mm256_srl_epi16(quotient_row, nshift_row_);
  __m256i odd_col(_mm256_and_si256(quotient_row, _mm256_set1_epi16(1)));
  /*
  odd_lay = _mm256_mullo_epi16(odd_lay, _mm256_set1_epi16(24));
  odd_col = _mm256_mullo_epi16(odd_col, _mm256_set1_epi16(12));
  */
  //combine 4 instructions below into:
  //left shift 8 bits (1 byte) of odd_lay
  //odd_lay = _mm256_slli_epi16(odd_lay, 8);
  //Option 1
  /*
  odd_lay = _mm256_slli_si256(odd_lay, 1);
  //OR odd_lay with odd_col
  odd_lay = _mm256_or_si256(odd_lay, odd_col);
  odd_lay = _mm256_maddubs_epi16(odd_lay, _mm256_set1_epi16(24));
  nrand = _mm256_add_epi16(nrand, odd_lay);
  */
  //Option 2 faster
  odd_lay = _mm256_sign_epi16(odd_lay, _mm256_set1_epi64x(-1));
  odd_lay = _mm256_and_si256(odd_lay, _mm256_set1_epi16(24));
  odd_col = _mm256_sign_epi16(odd_col, _mm256_set1_epi64x(-1));
  odd_col = _mm256_and_si256(odd_col, _mm256_set1_epi16(12));
  nrand = _mm256_add_epi16(nrand, odd_lay);
  nrand = _mm256_add_epi16(nrand, odd_col);
  /*
  odd_lay = _mm256_maddubs_epi16(odd_lay, _mm256_set1_epi16(24));
  odd_col = _mm256_maddubs_epi16(odd_col, _mm256_set1_epi16(12));
  nrand.m256i = _mm256_add_epi16(nrand.m256i, odd_col);
  nrand.m256i = _mm256_add_epi16(nrand.m256i, odd_lay);
  */
  /*
  union256 mvdx;
  mvdx.m256i = vdx;
  for(unsigned i(0); i != 16; ++i)
    {
      std::cout << "index:" << nrand.uint16[i];
      const bool bodd_lay((mvdx.uint16[i]/NUM_COLROW)&1);
      const bool bodd_col((mvdx.uint16[i]%NUM_COLROW/NUM_ROW)&1);
      std::cout << " actual:" << 
        arand.uint16[i]+(24&(-bodd_lay))+(12&(-bodd_col)) << std::endl;
    }
  exit(0);
  */
  //cast first 8 indices from uint16_t to uint32_t
  __m256i index(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(nrand)));
  //get the first 8 offsets
  index = _mm256_i32gather_epi32(offsets_, index, 4);
  //cast first 8 vdx from uint16_t to uint32_t and add with offset
  __m256i tar1 = _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(vdx)), index);
                                                     
  //cast second 8 indices from uint16_t to uint32_t
  index = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(nrand, 1));
  //get the second 8 offsets
  index =  _mm256_i32gather_epi32(offsets_, index, 4);
  //cast second 8 vdx from uint16_t to uint32_t and add with offset
  __m256i tar2 =
    _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_extractf128_si256(vdx, 1)),
                     index);
  tar2 = _mm256_packus_epi32(tar1, tar2);
  return _mm256_permute4x64_epi64(tar2, 216);
  /*
   Option 2: slower
  __m256i index(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(nrand)));
  //get the first 8 offsets
  index = _mm256_i32gather_epi32(offsets_, index, 4);
  //cast first 8 vdx from uint16_t to uint32_t and add with offset
  __m256i tar1 = _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(vdx)), index);
  __m128i tar1high = _mm256_extractf128_si256(tar1, 1);
                                                     
  //cast second 8 indices from uint16_t to uint32_t
  index = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(nrand, 1));
  //get the second 8 offsets
  index =  _mm256_i32gather_epi32(offsets_, index, 4);
  //cast second 8 vdx from uint16_t to uint32_t and add with offset
  __m256i tar2 =
    _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_extractf128_si256(vdx, 1)),
                     index);
  tar1 = _mm256_insertf128_si256(tar1, _mm256_castsi256_si128(tar2), 1);
  tar2 = _mm256_insertf128_si256(tar2, tar1high, 0);
  return _mm256_packus_epi32(tar1, tar2);
  */
  /*

  union256 mvdx;
  mvdx.m256i = vdx;
  for(unsigned i(0); i != 16; ++i)
    {
      std::cout << "tar:" << ((uint16_t*)&tar2)[i];
      const bool bodd_lay((mvdx.uint16[i]/NUM_COLROW)&1);
      const bool bodd_col((mvdx.uint16[i]%NUM_COLROW/NUM_ROW)&1);
      std::cout << " actual:" << mvdx.uint16[i]+offsets_[
        arand.uint16[i]+(24&(-bodd_lay))+(12&(-bodd_col))] << std::endl;
    }
  exit(0);
  */
}

__m256i Compartment::get_tars_exp(const __m256i vdx, __m256i nrand) const { 
  /*
  __m256i quotient_colrow(_mm256_mulhi_epu16(vdx, multiplier_colrow_));
  quotient_colrow = _mm256_srl_epi16(quotient_colrow, nshift_colrow_);
  __m256i odd_lay(_mm256_and_si256(quotient_colrow, _mm256_set1_epi16(1)));
  quotient_colrow = _mm256_mullo_epi16(quotient_colrow,
                                       _mm256_set1_epi16(NUM_COLROW));
  __m256i remainder_colrow(_mm256_sub_epi16(vdx, quotient_colrow));
  __m256i quotient_row(_mm256_mulhi_epu16(remainder_colrow, multiplier_row_));
  quotient_row = _mm256_srl_epi16(quotient_row, nshift_row_);
  __m256i odd_col(_mm256_and_si256(quotient_row, _mm256_set1_epi16(1)));
  odd_lay = _mm256_sign_epi16(odd_lay, _mm256_set1_epi64x(-1));
  odd_lay = _mm256_and_si256(odd_lay, _mm256_set1_epi16(24));
  odd_col = _mm256_sign_epi16(odd_col, _mm256_set1_epi64x(-1));
  odd_col = _mm256_and_si256(odd_col, _mm256_set1_epi16(12));
  nrand = _mm256_add_epi16(nrand, odd_lay);
  nrand = _mm256_add_epi16(nrand, odd_col);
  __m256i index(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(nrand)));
  index = _mm256_i32gather_epi32(offsets_, index, 4);
  __m256i tar1 = _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(vdx)), index);
  index = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(nrand, 1));
  index =  _mm256_i32gather_epi32(offsets_, index, 4);
  __m256i tar2 =
    _mm256_add_epi32(_mm256_cvtepu16_epi32(_mm256_extractf128_si256(vdx, 1)),
                     index);
  tar2 = _mm256_packus_epi32(tar1, tar2);
  return _mm256_permute4x64_epi64(tar2, 216);
  */
  //__m256i tar1 = _mm256_i32gather_epi32(offsets_, nrand, 2);
 // __m256i tar2 = _mm256_i32gather_epi32(offsets_, _mm256_srli_epi16(nrand, 2), 4);
 //index contains odd_col:odd_lay:odd_nrand x 16 
 //first get odd_nrand
 __m256i odd_nrand(_mm256_and_si256(nrand, _mm256_set1_epi16(1)));
 //set last bit of nrand as 0
 //nrand = [10, 8, 6, 4, 2, 0]
 nrand = _mm256_and_si256(nrand, _mm256_set1_epi16(0xfffe));
 //nrand*3 = [30, 24, 18, 12, 6, 0]
 __m256i nrand2(_mm256_add_epi16(nrand, nrand));
 nrand = _mm256_add_epi16(nrand2, nrand);
 //1.66 s

 __m256i odd_lay(_mm256_srli_epi16(
     _mm256_and_si256(vdx, _mm256_set1_epi16(32)), 4));
 __m256i odd_col(_mm256_srli_epi16(
     _mm256_and_si256(vdx, _mm256_set1_epi16(32)), 8));
 __m256i index(_mm256_or_si256(odd_nrand, odd_lay));
 index = _mm256_or_si256(index, odd_col);
 //1.72 s
 __m256i clr(_mm256_setr_epi32(1292979281,
      3429915664, 1339051024, 4231036115, 1276988432, 3480243411, 1288723537,
      4231285776));
  __m256i tar1 = _mm256_permutevar8x32_epi32(clr,
      _mm256_and_si256(index, _mm256_set1_epi32(7)));
  //1.80
  tar1 = _mm256_srlv_epi32(tar1, _mm256_and_si256(nrand, _mm256_set1_epi32(0xffff)));
  //1.91
  __m256i tar2 = _mm256_permutevar8x32_epi32(clr, _mm256_srli_epi32(index, 16));
  //2.00
  tar2 = _mm256_slli_epi32(
      _mm256_srlv_epi32(tar2, _mm256_srli_epi32(nrand, 16)), 16);
  //2.14
  __m256i clr6(_mm256_or_si256(tar1, tar2));
  //2.14
  clr6 = _mm256_xor_si256(clr6, _mm256_set1_epi16(21));
  //2.21
  __m256i lay(_mm256_slli_epi16(_mm256_and_si256(clr6, _mm256_set1_epi16(12)), 3));
  //2.28
  __m256i col(_mm256_slli_epi16(_mm256_and_si256(clr6, _mm256_set1_epi16(48)), 6));
  //2.39
  clr6 = _mm256_and_si256(clr6, _mm256_set1_epi16(3));
  clr6 = _mm256_or_si256(clr6, lay);
  //2.48
  clr6 = _mm256_or_si256(clr6, col);
  __m256i vdx2 = _mm256_sub_epi16(vdx, _mm256_set1_epi16(1));
  //2.50
  return _mm256_add_epi16(vdx2, clr6);

  /*
  cout_binary(test, "test");
  for(unsigned i(0); i != 8; ++i)
    {
      std::cout << ((uint32_t*)&test)[i] << " " << std::endl;
    }

  exit(0);
  */


  //return _mm256_add_epi16(clr6, _mm256_maddubs_epi16(tar1, tar2));
  //return clr6;
}
//_mm256_i32gather_epi32 = 2.84 s
//_mm256_shuffle_epi8 = 1.63 s (so gather takes about 1.2 s)
//2x_mm256_i32gather_epi32 = 4.1 s (1.63+2*1.2 = 4.1 s)


/*
  1 xor 1 = 0
  0 xor 1 = 1
  3 xor 1 = 2
  so 
    ori [-1, 0, 1] ->
    add [0, 1, 2] (to avoid -1 which requires signed addition) -> 
    xor with 1 [1, 0, 3] (to convert col:row from 01:01 to 00:00 to support the
                          missing 4 bits that become 0000 after right logical
                          shift)
     -1 -> 0 -> 1
      0 -> 1 -> 0
      1 -> 2 -> 3
  odd_col:odd_lay:odd_nrand = 0:0:0 [col=even:layer=even:nrand=even]
  offsets_[2].clr = (1, 0, 1) = 010001; 6 bits 
  offsets_[4].clr = (3, 0, 1) = 110001; 6 bits
  offsets_[6].clr = (1, 1, 0) = 010100; 6 bits
  offsets_[8].clr = (0, 1, 0) = 000100; 6 bits
  offsets_[10].clr = (0, 3, 1) = 001101; 6 bits
  offsets_[0].clr = (0, 0, 1) = 01; 2 bits (first 4 bits are always 0000) 
  total = 2+5*6 = 32 bits = 01001101000100010100110001010001 = 1292979281

  col:lay:nrand = 0:0:1 [col=even:layer=even:nrand=odd]
  offsets_[3].clr = (1, 0, 0) = 010000;
  offsets_[5].clr = (3, 0, 0) = 110000;
  offsets_[7].clr = (0, 1, 1) = 000101;
  offsets_[9].clr = (1, 3, 0) = 011100;
  offsets_[11].clr = (0, 3, 0) = 001100;
  offsets_[1].clr = (0, 0, 3) = 11;
  total = 11001100011100000101110000010000 = 3429915664

  col:lay:nrand = 0:1:0 [col=even:layer=odd:nrand=even]
  offsets_[26].clr = (1, 0, 0) = 010000;
  offsets_[28].clr = (3, 0, 0) = 110000;
  offsets_[30].clr = (0, 1, 0) = 000100;
  offsets_[32].clr = (3, 1, 0) = 110100;
  offsets_[34].clr = (0, 3, 3) = 001111;
  offsets_[24].clr = (0, 0, 1) = 01;
  total = 01001111110100000100110000010000 = 1339051024 

  col:lay:nrand = 0:1:1 [col=even:layer=odd:nrand=odd]
  offsets_[27].clr = (1, 0, 3) = 010011;
  offsets_[29].clr = (3, 0, 3) = 110011;
  offsets_[31].clr = (0, 1, 3) = 000111;
  offsets_[33].clr = (0, 3, 0) = 001100;
  offsets_[35].clr = (3, 3, 0) = 111100;
  offsets_[25].clr = (0, 0, 3) = 11;
  total = 11111100001100000111110011010011 = 4231036115

  col:lay:nrand = 1:0:0 [col=odd:layer=even:nrand=even]
  offsets_[14].clr = (1, 0, 0) = 010000;
  offsets_[16].clr = (3, 0, 0) = 110000;
  offsets_[18].clr = (1, 1, 0) = 010100;
  offsets_[20].clr = (0, 1, 3) = 000111;
  offsets_[22].clr = (0, 3, 0) = 001100;
  offsets_[12].clr = (0, 0, 1) = 01;
  total = 01001100000111010100110000010000 = 1276988432

  col:lay:nrand = 1:0:1 [col=odd:layer=even:nrand=odd]
  offsets_[15].clr = (1, 0, 3) = 010011;
  offsets_[17].clr = (3, 0, 3) = 110011;
  offsets_[19].clr = (0, 1, 0) = 000100;
  offsets_[21].clr = (1, 3, 0) = 011100;
  offsets_[23].clr = (0, 3, 3) = 001111;
  offsets_[13].clr = (0, 0, 3) = 11;
  total = 11001111011100000100110011010011 = 3480243411

  col:lay:nrand = 1:1:0 [col=odd:layer=odd:nrand=even]
  offsets_[38].clr = (1, 0, 1) = 010001;
  offsets_[40].clr = (3, 0, 1) = 110001;
  offsets_[42].clr = (0, 1, 1) = 000101;
  offsets_[44].clr = (3, 1, 0) = 110100;
  offsets_[46].clr = (0, 3, 0) = 001100;
  offsets_[36].clr = (0, 0, 1) = 01;
  total = 01001100110100000101110001010001 = 1288723537

  col:lay:nrand = 1:1:1 [col=odd:layer=odd:nrand=odd]
  offsets_[39].clr = (1, 0, 0) = 010000;
  offsets_[41].clr = (3, 0, 0) = 110000;
  offsets_[43].clr = (0, 1, 0) = 000100;
  offsets_[45].clr = (0, 3, 1) = 001101;
  offsets_[47].clr = (3, 3, 0) = 111100;
  offsets_[37].clr = (0, 0, 3) = 11;
  total = 11111100001101000100110000010000 = 4231285776
  _mm256_setr_epi32(1292979281, 3429915664, 1339051024, 4231036115, 1276988432, 3480243411, 1288723537, 4231285776);
*/


void Compartment::set_offsets() {
  //col=even, layer=even
  offsets_ = new int[ADJS*4];
  offsets_[0] = -1;
  offsets_[1] = 1;
  offsets_[2] = -NUM_ROW-1;
  offsets_[3] = -NUM_ROW;
  offsets_[4] = NUM_ROW-1;
  offsets_[5] = NUM_ROW;
  offsets_[6] = -NUM_COLROW-NUM_ROW;
  offsets_[7] = -NUM_COLROW-1;
  offsets_[8] = -NUM_COLROW;
  offsets_[9] = NUM_COLROW-NUM_ROW;
  offsets_[10] = NUM_COLROW-1;
  offsets_[11] = NUM_COLROW;

  //col=even, layer=odd +24 = %layer*24
  offsets_[24] = -1;
  offsets_[25] = 1;
  offsets_[26] = -NUM_ROW;
  offsets_[27] = -NUM_ROW+1;
  offsets_[28] = NUM_ROW;
  offsets_[29] = NUM_ROW+1;
  offsets_[30] = -NUM_COLROW;
  offsets_[31] = -NUM_COLROW+1;
  offsets_[32] = -NUM_COLROW+NUM_ROW;
  offsets_[33] = NUM_COLROW;
  offsets_[34] = NUM_COLROW+1;
  offsets_[35] = NUM_COLROW+NUM_ROW;

  //col=odd, layer=even +12 = %col*12
  offsets_[12] = -1;
  offsets_[13] = 1;
  offsets_[14] = -NUM_ROW;
  offsets_[15] = -NUM_ROW+1;
  offsets_[16] = NUM_ROW;
  offsets_[17] = NUM_ROW+1;
  offsets_[18] = -NUM_COLROW-NUM_ROW;
  offsets_[19] = -NUM_COLROW;
  offsets_[20] = -NUM_COLROW+1;
  offsets_[21] = NUM_COLROW-NUM_ROW;
  offsets_[22] = NUM_COLROW;
  offsets_[23] = NUM_COLROW+1;

  //col=odd, layer=odd +36 = %col*12 + %layer*24
  offsets_[36] = -1;
  offsets_[37] = 1;
  offsets_[38] = -NUM_ROW-1;
  offsets_[39] = -NUM_ROW;
  offsets_[40] = NUM_ROW-1;
  offsets_[41] = NUM_ROW;
  offsets_[42] = -NUM_COLROW-1;
  offsets_[43] = -NUM_COLROW; //a
  offsets_[44] = -NUM_COLROW+NUM_ROW;
  offsets_[45] = NUM_COLROW-1;
  offsets_[46] = NUM_COLROW;
  offsets_[47] = NUM_COLROW+NUM_ROW;
}

const Vector<double>& Compartment::get_dimensions() const {
  return dimensions_;
}

Vector<double> Compartment::get_center() const {
  return dimensions_/2;
}

Species& Compartment::get_surface_species() {
  return surface_species_;
}

Species& Compartment::get_volume_species() {
  return volume_species_;
}

Model& Compartment::get_model() {
  return model_;
}

const std::string& Compartment::get_name() const {
  return name_;
}

Lattice& Compartment::get_lattice() {
  return lattice_;
}

void Compartment::set_volume_structure() {
  for(unsigned box(0); box != get_lattice().get_num_box(); ++box) {
    for(umol_t i(0); i != get_lattice().get_num_box_voxel(); ++i) {
      get_volume_species().populate_mol(box, i);
    }
  }
}

void Compartment::set_surface_structure() {
  for(unsigned box(0); box != get_lattice().get_num_box(); ++box) {
    //row_col xy-plane
    for (umol_t i(0); i != NUM_COLROW; ++i) {
      get_surface_species().populate_mol(box, i);
      get_surface_species().populate_mol(box, NUM_VOXEL-1-i);
    }
    for (umol_t i(1); i != NUM_LAY-1; ++i) {
      //layer_row yz-plane
      for (umol_t j(0); j != NUM_ROW; ++j) {
        get_surface_species().populate_mol(box, i*NUM_COLROW+j);
        get_surface_species().populate_mol(box, i*NUM_COLROW+j+NUM_ROW*
                                           (NUM_COL-1));
      }
      //layer_col xz-plane
      for (umol_t j(1); j != NUM_COL-1; ++j) {
        get_surface_species().populate_mol(box, i*NUM_COLROW+j*NUM_ROW);
        get_surface_species().populate_mol(box, i*NUM_COLROW+j*NUM_ROW+
                                           NUM_ROW-1);
      }
    }
  }
}


