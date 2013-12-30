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

#include <math.h>
#include <Compartment.hpp>
#include <Species.hpp>
#include <Model.hpp>


Compartment::Compartment(std::string name, const double len_x,
                         const double len_y, const double len_z, Model& model)
  : name_(name),
    length_(len_x, len_y, len_z),
    center_(len_x/2, len_y/2, len_z/2),
    model_(model),
    volume_species_("volume", 0, 0, model, *this, volume_species_, true),
    surface_species_("surface", 0, 0, model, *this, surface_species_, true) {}

void Compartment::initialize() {
  setOffsets();
  nbit_ = model_.get_nbit();
  sur_xor_ = surface_species_.get_id()^volume_species_.get_id();
  lattice_.resize(ceil(double(NUM_VOXEL)*nbit_/WORD), 0);
  set_surface();
  std::cout << "nrow:" << NUM_ROW << " ncol:" << NUM_COL << " nlay:" <<
    NUM_LAY << " latticeSize:" << lattice_.size() << " memory:" << 
    lattice_.size()*sizeof(unsigned)/(1024*1024.0) << " MB" << std::endl;
}

unsigned Compartment::get_num_col() const {
  return NUM_COL;
}

unsigned Compartment::get_num_lay() const {
  return NUM_LAY;
}

unsigned Compartment::get_num_row() const {
  return NUM_ROW;
}

unsigned Compartment::get_num_voxel() const {
  return NUM_VOXEL;
}

/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  const bool odd_lay((vdx/NUM_COLROW)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-NUM_ROW-1 ;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx+NUM_ROW*(odd_lay-NUM_COL-1)-(odd_col&odd_lay);
    case 7:
      return vdx+!odd_col*(odd_lay-!odd_lay)-NUM_COLROW;
    case 8:
      return vdx+NUM_ROW*(odd_lay-NUM_COL)+(odd_col&!odd_lay);
    case 9:
      return vdx+NUM_ROW*(NUM_COL-!odd_lay)-(odd_col&odd_lay);
    case 10:
      return vdx+NUM_COLROW+!odd_col*(odd_lay-!odd_lay);
    case 11:
      return vdx+NUM_ROW*(NUM_COL+odd_lay)+(odd_col&!odd_lay);
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  unsigned value(0);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      value = -1;
    case 3:
      return vdx+(odd_col^odd_lay)-202+value;
    case 4:
      value = -1;
    case 5:
      return vdx+(odd_col^odd_lay)+202+value;
    case 6:
      value = 94132;
    case 7:
      return vdx+value-47066-(odd_col&odd_lay)-202*(!odd_lay);
    case 8:
      value = 94132;
    case 9:
      return vdx+value-47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 10:
      value = 94132;
    case 11:
      return vdx+value-47066+(odd_col&!odd_lay)+202*odd_lay;
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-nrowm_;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+nrowm_;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx-NUM_COLROW-(odd_col&odd_lay)-NUM_ROW*(!odd_lay);
    case 7:
      return vdx-NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-NUM_COLROW+(odd_col&!odd_lay)+NUM_ROW*odd_lay;
    case 9:
      return vdx+NUM_COLROW-(odd_col&odd_lay)+NUM_ROW*(!odd_lay);
    case 10:
      return vdx+NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+NUM_COLROW+(odd_col&!odd_lay)+NUM_ROW*odd_lay;
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-201;
    case 3:
      return vdx+(odd_col^odd_lay)-202;
    case 4:
      return vdx+(odd_col^odd_lay)+201;
    case 5:
      return vdx+(odd_col^odd_lay)+202;
    case 6:
      return vdx-47066-(odd_col&odd_lay)-(202&(-!odd_lay));
    case 7:
      return vdx-47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-47066+(odd_col&!odd_lay)+(202&(-odd_lay));
    case 9:
      return vdx+47066-(odd_col&odd_lay)-(202&(-!odd_lay));
    case 10:
      return vdx+47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+47066+(odd_col&!odd_lay)+(202&(-odd_lay));
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const {
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  const bool odd_lay((vdx/NUM_COLROW)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-NUM_ROW-1;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx-NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 7:
      return vdx-NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    case 9:
      return vdx+NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 10:
      return vdx+NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    }
  return vdx-1;
}
*/

unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const {
  const bool odd_lay((vdx/NUM_COLROW)&1);
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  //return vdx+offsets_[nrand+odd_lay*24+odd_col*12];
  return vdx+offsets_[nrand+(24&(-odd_lay))+(12&(-odd_col))];
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
}

void Compartment::set_tars(const unsigned vdx, union256& nrand) const {
  nrand.m256i = _mm256_i32gather_epi32(offsets_, nrand.m256i, 4);
  std::cout << "num_colrow:" << NUM_COLROW << " num_row:" << NUM_ROW << std::endl;
  //num_colrow: 47066
  //num_row: 202
}

void Compartment::setOffsets() {
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

const Vector& Compartment::get_center() const {
  return center_;
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

std::vector<unsigned>& Compartment::get_lattice() {
  return lattice_;
}

void Compartment::set_surface() {
  //row_col xy-plane
  for (unsigned i(0); i != NUM_COLROW; ++i) {
      populate_mol(i);
      populate_mol(NUM_VOXEL-1-i);
    }
  for (unsigned i(1); i != NUM_LAY-1; ++i) {
      //layer_row yz-plane
      for (unsigned j(0); j != NUM_ROW; ++j) {
          populate_mol(i*NUM_COLROW+j);
          populate_mol(i*NUM_COLROW+j+NUM_ROW*(NUM_COL-1));
        }
      //layer_col xz-plane
      for (unsigned j(1); j != NUM_COL-1; ++j) {
          populate_mol(i*NUM_COLROW+j*NUM_ROW);
          populate_mol(i*NUM_COLROW+j*NUM_ROW+NUM_ROW-1);
        }
    }
  //std::cout << "surface size:" << mols.size() << " actual size:" <<
  //NUM_COLROW*2 + NUM_ROW*(NUM_LAY-2)*2 + (NUM_COL-2)*(NUM_LAY-2)*2 <<
  //std::endl;
}

void Compartment::populate_mol(const unsigned vdx) {
  lattice_[vdx*nbit_/WORD] ^= sur_xor_ << vdx*nbit_%WORD;
  surface_species_.get_mols().push_back(vdx);
}

