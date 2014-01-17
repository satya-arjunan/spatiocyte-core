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

#include <time.h>
#include <Diffuser.hpp>
#include <Compartment.hpp>
#include <Model.hpp>

Diffuser::Diffuser(const double D, Species& species):
  D_(D),
  species_(species),
  compartment_(species_.get_compartment()),
  box_mols_(species_.get_box_mols()),
  species_id_(species_.get_id()),
  vac_id_(species_.get_vac_id()),
  vac_xor_(species_.get_vac_xor()),
  rng_(time(0)) {}

void Diffuser::initialize() {
  box_voxels_ = compartment_.get_lattice().get_box_voxels();
  nbit_ = species_.get_compartment().get_model().get_nbit();
  one_nbit_ = pow(2, nbit_)-1;
  std::cout << species_.get_name_id() << std::endl;
  std::cout << "nbit:" << nbit_ << std::endl;
  std::cout << "one_nbit:" << one_nbit_ << std::endl;
  std::cout << "vac_id:" << vac_id_ << std::endl;
  std::cout << "vac_xor:" << vac_xor_ << std::endl;
}

//t(10000, 16) = 10.97
void Diffuser::walk() {
  for(unsigned box(0), n(box_mols_.size()); box != n; ++box) {
    walk((__m256i*)(&box_mols_[box][0]), box_mols_[box].size());
  }
}

void Diffuser::walk(__m256i* base, const unsigned size) {
  __m256i mols_m256i(_mm256_load_si256(base));
  __m256i tars(compartment_.get_tars(mols_m256i, rng_.Ran16()));
  __m256i cmps(_mm256_cmpeq_epi16(tars, mols_m256i));
  for (unsigned i(0); i != 15; ++i) {
    tars = _mm256_or_si256(_mm256_slli_si256(tars, 2),
                           _mm256_srli_si256(tars, 30)); 
    cmps = _mm256_or_si256(_mm256_slli_si256(cmps, 2),
                           _mm256_srli_si256(cmps, 30)); 
    cmps = _mm256_or_si256(cmps, _mm256_cmpeq_epi16(tars, mols_m256i));
  }
  cmps = _mm256_or_si256(_mm256_slli_si256(cmps, 2),
                         _mm256_srli_si256(cmps, 30)); 
  mols_m256i = _mm256_and_si256(mols_m256i, cmps);
  cmps = _mm256_andnot_si256(cmps, cmps);
  tars = _mm256_and_si256(tars, cmps);
  mols_m256i = _mm256_or_si256(mols_m256i, tars);
  _mm256_store_si256(base, cmps);
}

/*
// t = 43 s
void Diffuser::walk(__m256i* base, const unsigned size) {
  __m256i mols_m256i(_mm256_load_si256(base));
  __m256i tars(compartment_.get_tars(mols_m256i, rng_.Ran16()));
  __m256i cmps(_mm256_cmpeq_epi16(tars, mols_m256i));
  for (unsigned i(0); i != 15; ++i) {
    //rotate tars
    const umol_t tar(((umol_t*)&tars)[0]);
    tars = _mm256_slli_si256(tars, 2); 
    //((umol_t*)&tars)[0] = tar;
    ((umol_t*)&tars)[0] |= tar;
    //rotate cmps
    const umol_t cmp(((umol_t*)&cmps)[0]);
    cmps = _mm256_slli_si256(cmps, 2); 
    //((umol_t*)&cmps)[0] = cmp;
    ((umol_t*)&cmps)[0] |= cmp;
    cmps = _mm256_or_si256(cmps, _mm256_cmpeq_epi16(tars, mols_m256i));
  }
  const umol_t cmp(((umol_t*)&cmps)[0]);
  cmps = _mm256_slli_si256(cmps, 2); 
  //((umol_t*)&cmps)[0] = cmp;
  ((umol_t*)&cmps)[0] |= cmp;
  mols_m256i = _mm256_and_si256(mols_m256i, cmps);
  cmps = _mm256_andnot_si256(cmps, cmps);
  tars = _mm256_and_si256(tars, cmps);
  mols_m256i = _mm256_or_si256(mols_m256i, tars);
  _mm256_store_si256(base, mols_m256i);
}
*/

/*
t = 27 s
void Diffuser::walk() {
  for(unsigned box(0), n(box_mols_.size()); box != n; ++box) {
    walk((__m256i*)(&box_mols_[box][0]), box_mols_[box].size());
  }
}

void Diffuser::walk(__m256i* base, const unsigned size) {
  __m256i mols_m256i(_mm256_load_si256(base));
  for (unsigned j(0); j != 16; ++j) {
    umol_t tar(compartment_.get_tar(((umol_t*)&mols_m256i)[j], rng_.Ran16_12()));
    __m256i tars(_mm256_set1_epi16(tar));
    __m256i res(_mm256_cmpeq_epi16(tars, mols_m256i));
    if (!_mm256_movemask_epi8(res)) {
      ((umol_t*)&mols_m256i)[j] = tar;
    }
  }
  _mm256_store_si256(base, mols_m256i);
}
*/

/*
void Diffuser::walk(__m256i* base, const unsigned size) {
  const unsigned n(size/16);
  __m256i* mols(base);
  for (unsigned k(0); k != n; ++k, ++mols) {
    __m256i mols_m256i(_mm256_load_si256(mols));
    for (unsigned j(0); j != 16; ++j) {
      __m256i* mols2(base);
      umol_t tar(compartment_.get_tar(((umol_t*)&mols)[j], rng_.Ran16_12()));
      __m256i tars(_mm256_set1_epi16(tar));
      for (unsigned i(0); i != n; ++i, ++mols2) {
        __m256i res(_mm256_cmpeq_epi16(tars, _mm256_load_si256(mols2)));
        if (_mm256_movemask_epi8(res)) {
          goto next;
        }
      }
      ((umol_t*)&mols_m256i)[j] = tar;
next:
    }
    _mm256_store_si256(mols, mols_m256i);
  }
}
*/
//VPMASKMOVD void _mm256_maskstore_epi32(int *a, __m256i mask, __m256i b)
//VMOVDQA void _mm256_store_si256(__m256i *a, __m256i b)
//VPCMPEQB __m256i _mm256_cmpeq_epi8 ( __m256i a, __m256i b)
//VPCMPEQW __m256i _mm256_cmpeq_epi16 ( __m256i a, __m256i b)


/*
//t(10000, 16) = 56 s
//t(10000, 1) = 1.87 s
void Diffuser::walk() {
  for(unsigned box(0), n(box_mols_.size()); box != n; ++box) {
    walk(&box_mols_[box][0], box_mols_[box].size());
  }
}

void Diffuser::walk(umol_t* mols, const unsigned size) {
  for (unsigned i(0); i != size; ++i) {
    umol_t tar(compartment_.get_tar(mols[i], rng_.Ran16_12()));
    for(unsigned j(0); j != size; ++j) {
      if(mols[j] == tar) {
        goto next;
      }
    }
    mols[i] = tar;
next:
    continue;
  }
}
*/


/*
//t(10000, 16) = 94 s
//t(1, 10000) = 1.11 s
void Diffuser::walk() {
  for(unsigned box(0), n(box_mols_.size()); box != n; ++box) {
    walk(box_voxels_[box], (__m256i*)(&box_mols_[box][0]),
        box_mols_[box].size());
  }
}

void Diffuser::walk(voxel_t* voxels, __m256i* mols, const unsigned size) {
  const unsigned n(size/16);
  const unsigned m(size%16);
  uint32_t tars[16];
  for (unsigned k(0); k != n; ++k, ++mols) {
    __m256i mols_m256i(_mm256_load_si256(mols));
    compartment_.set_tars(mols_m256i, rng_.Ran16(), tars);
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (voxels[vdx] == vac_id_) {
        voxels[((umol_t*)&mols_m256i)[j]] = vac_id_;
        voxels[vdx] = species_id_;
        ((umol_t*)&mols_m256i)[j] = vdx;
      }
    }
    _mm256_store_si256(mols, mols_m256i);
  }
  if(m) {
    compartment_.set_tars(_mm256_load_si256(mols), rng_.Ran16(), tars);
    for (unsigned j(0); j != m; ++j) {
      const uint32_t vdx(tars[j]);
      if (voxels[vdx] == vac_id_) {
        voxels[((umol_t*)mols)[j]] = vac_id_;
        voxels[vdx] = species_id_;
        ((umol_t*)mols)[j] = vdx;
      }
    }
  }
}
*/

/*
static uint8_t data[H][W]; for(int i=1; i<H-1; ++i) for(int j=1; j<W-1; ++j) {
border = 1; for(int k=-1;k<=1;++k) for(int m=-1;m<=1;++m) border |= data[i+k][j+m] == 0;
if(border)
}
ProcessElement(i, j);

// Do for every line
__m256i l2 = _mm256_loadu_si256((__m256i*)src_ptr);
int mdl = _mm256_movemask_epi8( _mm256_cmpeq_epi8(l2, vzero) );
// Check if mdl[1..30] == 0 - all pixels are masked out
if( (~mdl & 0x7FFFFFFE) != 0) {
  __m256i l1 = _mm256_loadu_si256((__m256i *) (src_ptr-W));
  __m256i l3 = _mm256_loadu_si256((__m256i *) (src_ptr+W));
  int up = _mm256_movemask_epi8(_mm256_cmpeq_epi8(l1, vzero));
  int dwn = _mm256_movemask_epi8(_mm256_cmpeq_epi8(l3, vzero));
  int border = up | mdl | dwn;
  border = (border << 1) | border | (border >> 1);
  // Process non border elements
  for(int t=1; t <= 30; ++t)
    if(border & (1 << t))
      ProcessElement(i, j+t);
}
src_ptr += 30;
*/

/*
//t = 1.115
void Diffuser::walk() {
  walk(lattice_, (__m256i*)(&mols_[0]), mols_.size());
}

void Diffuser::walk(voxel_t* voxels, __m256i* base, const unsigned size) {
  const unsigned n(size/16);
  const unsigned m(size%16);
  uint32_t tars[16];
  for (unsigned k(0); k != n; ++k, ++base) {
    __m256i mols_m256i(_mm256_load_si256(base));
    compartment_.set_tars(mols_m256i, rng_.Ran16(), tars);
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (voxels[vdx] == vac_id_) {
        voxels[((umol_t*)&mols_m256i)[j]] = vac_id_;
        voxels[vdx] = species_id_;
        ((umol_t*)&mols_m256i)[j] = vdx;
      }
    }
    _mm256_store_si256(base, mols_m256i);
  }
  compartment_.set_tars(_mm256_load_si256(base), rng_.Ran16(), tars);
  for (unsigned j(0); j != m; ++j) {
    const uint32_t vdx(tars[j]);
    if (voxels[vdx] == vac_id_) {
      voxels[((umol_t*)base)[j]] = vac_id_;
      voxels[vdx] = species_id_;
      ((umol_t*)base)[j] = vdx;
    }
  }
}
*/

/*
//t = 1.124
void Diffuser::walk() {
  walk(lattice_, mols_);
}

void Diffuser::walk(voxel_t* voxels, std::vector<umol_t>& mols) {
  const unsigned n(mols.size()/16);
  const unsigned m(mols.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols[0]));
  for (unsigned k(0); k != n; ++k, ++base) {
    __m256i mols_m256i(_mm256_load_si256(base));
    compartment_.set_tars(mols_m256i, rng_.Ran16(), tars);
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (voxels[vdx] == vac_id_) {
        voxels[((umol_t*)&mols_m256i)[j]] = vac_id_;
        voxels[vdx] = species_id_;
        ((umol_t*)&mols_m256i)[j] = vdx;
      }
    }
    _mm256_store_si256(base, mols_m256i);
  }
  compartment_.set_tars(_mm256_load_si256(base), rng_.Ran16(), tars);
  for (unsigned j(0); j != m; ++j) {
    const uint32_t vdx(tars[j]);
    if (voxels[vdx] == vac_id_) {
      voxels[((umol_t*)base)[j]] = vac_id_;
      voxels[vdx] = species_id_;
      ((umol_t*)base)[j] = vdx;
    }
  }
}
*/

/*
//t = 1.165 s
void Diffuser::walk() {
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  for (unsigned k(0); k != n; ++k, ++base) {
    __m256i mols(_mm256_load_si256(base));
    compartment_.set_tars(mols, rng_.Ran16(), tars);
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (lattice_[vdx] == vac_id_) {
        lattice_[((umol_t*)&mols)[j]] = vac_id_;
        lattice_[vdx] = species_id_;
        ((umol_t*)&mols)[j] = vdx;
      }
    }
    _mm256_store_si256(base, mols);
  }
  __m256i mols(_mm256_load_si256(base));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for (unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i) {
    const uint32_t vdx(tars[j]);
    if (lattice_[vdx] == vac_id_) {
      lattice_[mols_[i]] = vac_id_;
      lattice_[vdx] = species_id_;
      mols_[i] = vdx;
    }
  }
}
*/

/*
t = 1.20
void Diffuser::walk() {
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  const int* base2((int*)(&lattice_[0]));
  __m256i vox[2];
  for (unsigned k(0); k != n; ++k, ++base) {
    __m256i mols(_mm256_load_si256(base));
    compartment_.set_tars(mols, rng_.Ran16(), tars);
    vox[0] = (_mm256_i32gather_epi32(base2, *(__m256i*)(&tars[0]), 1));
    vox[1] = (_mm256_i32gather_epi32(base2, *(__m256i*)(&tars[8]), 1));
    for (unsigned j(0); j != 16; ++j) {
      if(((uint8_t*)&vox[0])[j*4] == vac_id_) {
        lattice_[((umol_t*)&mols)[j]] = vac_id_;
        ((umol_t*)&mols)[j] = tars[j];
      }
    }
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (lattice_[vdx] == vac_id_) {
        lattice_[vdx] = species_id_;
      }
    }
    _mm256_store_si256(base, mols);
  }
  __m256i mols(_mm256_load_si256(base));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for (unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i) {
    const uint32_t vdx(tars[j]);
    if (lattice_[vdx] == vac_id_) {
      lattice_[vdx] = species_id_;
      lattice_[mols_[i]] = vac_id_;
      mols_[i] = vdx;
    }
  }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  uint32_t tars[16];
  const int* base((int*)(&lattice_[0]));
  for(unsigned k(0); k != n; ++k)
    {
      compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
      __m256i vox(_mm256_i32gather_epi32(base, *(__m256i*)(&tars[0]), 1));
      for(unsigned j(0); j != 8; ++j, ++i)
        {
          const uint8_t voxel(((uint8_t*)&vox)[j*4]);
          if(voxel == vac_id_)
            {
              const uint32_t vdx(tars[j]);
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
      vox = (_mm256_i32gather_epi32(base, *(__m256i*)(&tars[8]), 1));
      for(unsigned j(0); j != 8; ++j, ++i)
        {
          const uint8_t voxel(((uint8_t*)&vox)[j*4]);
          if(voxel == vac_id_)
            {
              const uint32_t vdx(tars[j+8]);
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
//t = 1.18 s
void Diffuser::walk() {
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  for (unsigned k(0); k != n; ++k, ++base) {
    __m256i mols(_mm256_load_si256(base));
    compartment_.set_tars(mols, rng_.Ran16(), tars);
    for (unsigned j(0); j != 16; ++j) {
      const uint32_t vdx(tars[j]);
      if (lattice_[vdx] == vac_id_) {
        lattice_[vdx] = species_id_;
        lattice_[((umol_t*)&mols)[j]] = vac_id_;
        ((umol_t*)&mols)[j] = vdx;
      }
    }
    _mm256_store_si256(base, mols);
  }
  __m256i mols(_mm256_load_si256(base));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for (unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i) {
    const uint32_t vdx(tars[j]);
    if (lattice_[vdx] == vac_id_) {
      lattice_[vdx] = species_id_;
      lattice_[mols_[i]] = vac_id_;
      mols_[i] = vdx;
    }
  }
}
*/

/*
//t = 1.23 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  //uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  for(unsigned k(0); k != n; ++k, ++base)
    {
      __m256i mols(_mm256_load_si256(base));
      const __m256i tars(compartment_.get_tars(mols, rng_.Ran16()));
      for(unsigned j(0); j != 16; ++j)
        {
          if(lattice_[((uint16_t*)&tars)[j]] == vac_id_)
            {
              lattice_[((uint16_t*)&tars)[j]] = species_id_;
              lattice_[((umol_t*)&mols)[j]] = vac_id_;
              ((umol_t*)&mols)[j] = ((uint16_t*)&tars)[j];
            }
        }
      _mm256_store_si256(base, mols);
    }
  const __m256i mols(_mm256_load_si256(base));
  const __m256i tars(compartment_.get_tars(mols, rng_.Ran16()));
  for(unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i)
    {
      if(lattice_[((uint16_t*)&tars)[j]] == vac_id_)
        {
          lattice_[((uint16_t*)&tars)[j]] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = ((uint16_t*)&tars)[j];
        }
    }
}
*/

/*
t = 1.20 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  for(unsigned k(0); k != n; ++k, ++base)
    {
      __m256i mols(_mm256_load_si256(base));
      compartment_.set_tars(mols, rng_.Ran16(), tars);
      for(unsigned j(0); j != 16; ++j)
        {
          const uint32_t vdx(tars[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[((umol_t*)&mols)[j]] = vac_id_;
              ((umol_t*)&mols)[j] = vdx;
            }
        }
      _mm256_store_si256(base, mols);
    }
  __m256i mols(_mm256_load_si256(base));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for(unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
t = 1.24 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  uint32_t tars[16];
  __m256i* base((__m256i*)(&mols_[0]));
  for(unsigned k(0); k != n; ++k, ++base)
    {
      __m256i mols(_mm256_load_si256(base));
      compartment_.set_tars(mols, rng_.Ran16(), tars);
      for(unsigned j(0); j != 16; ++j)
        {
          const uint32_t vdx(tars[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[((umol_t*)&mols)[j]] = vac_id_;
              ((umol_t*)&mols)[j] = vdx;
            }
        }
      _mm256_store_si256(base, mols);
    }
  __m256i mols(_mm256_load_si256(base));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for(unsigned j(0), i(mols_.size()-m); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
t = 1.41 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  uint32_t tars[16];
  for(unsigned k(0); k != n; ++k)
    {
      __m256i mols(_mm256_load_si256((__m256i*)(&mols_[i])));
      compartment_.set_tars(mols, rng_.Ran16(), tars);
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint32_t vdx(tars[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  __m256i mols(_mm256_load_si256((__m256i*)(&mols_[i])));
  compartment_.set_tars(mols, rng_.Ran16(), tars);
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
t = 1.48 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  uint32_t tars[16];
  const int* base((int*)(&lattice_[0]));
  for(unsigned k(0); k != n; ++k)
    {
      compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
      __m256i vox(_mm256_i32gather_epi32(base, *(__m256i*)(&tars[0]), 1));
      for(unsigned j(0); j != 8; ++j, ++i)
        {
          const uint8_t voxel(((uint8_t*)&vox)[j*4]);
          if(voxel == vac_id_)
            {
              const uint32_t vdx(tars[j]);
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
      vox = (_mm256_i32gather_epi32(base, *(__m256i*)(&tars[8]), 1));
      for(unsigned j(0); j != 8; ++j, ++i)
        {
          const uint8_t voxel(((uint8_t*)&vox)[j*4]);
          if(voxel == vac_id_)
            {
              const uint32_t vdx(tars[j+8]);
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
t = 1.47 s
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  uint32_t tars[16];
  for(unsigned k(0); k != n; ++k)
    {
      compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint32_t vdx(tars[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  compartment_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
      __m256i quos(_mm256_srli_epi16(tars, 4));
      __m256i rems(_mm256_slli_epi16(_mm256_sub_epi16(tars, _mm256_slli_epi16(quos, 4)), 1));
      //cast first 8 quos from uint16_t to uint32_t 
      __m256i rem(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(rems)));
      __m256i quo(_mm256_cvtepu16_epi32(_mm256_castsi256_si128(quos)));
      __m256i vox(_mm256_i32gather_epi32(lattice_, quo, 4));
      __m256i val(_mm256_srlv_epi32(vox, rem));
      val = _mm256_and_si256(val,_mm256_set1_epi32(3));
      for(unsigned j(0); j != 8; ++j, ++i)
        {
          const uint32_t voxel(((uint32_t*)&val)[j]);
          if(!voxel)
            {
              const uint16_t quo1(((uint16_t*)&quos)[j]);
              const uint16_t rem1(((uint16_t*)&rems)[j]);
              const uint16_t vdx1(((uint16_t*)&tars)[j]);
              lattice_[quo1] ^= 2 << rem1;
              lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
              //mols_[i] = ((uint16_t*)&tars)[j];
              mols_[i] = vdx1;
            }
        }
      //cast second 8 quos from uint16_t to uint32_t 
      rem = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(rems, 1));
      quo = _mm256_cvtepu16_epi32(_mm256_extractf128_si256(quos, 1));
      vox = _mm256_i32gather_epi32(lattice_, quo, 4);
      val = _mm256_srlv_epi32(vox, rem);
      val = _mm256_and_si256(val,_mm256_set1_epi32(3));
      for(unsigned j(8); j != 16; ++j, ++i)
        {
          const uint32_t voxel(((uint32_t*)&val)[j-8]);
          if(!voxel)
            {
              const uint16_t quo1(((uint16_t*)&quos)[j]);
              const uint16_t rem1(((uint16_t*)&rems)[j]);
              const uint16_t vdx1(((uint16_t*)&tars)[j]);
              lattice_[quo1] ^= 2 << rem1;
              lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
              //mols_[i] = ((uint16_t*)&tars)[j];
              mols_[i] = vdx1;
            }
        }
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
        {
          lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
          lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
      __m256i quos(_mm256_srli_epi16(tars, 4));
      __m256i rems(_mm256_slli_epi16(_mm256_sub_epi16(tars, _mm256_slli_epi16(quos, 4)), 1));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint16_t quo(((uint16_t*)&quos)[j]);
          const uint16_t rem(((uint16_t*)&rems)[j]);
          //const uint16_t vdx(((uint16_t*)&tars)[j]);
          if(0 == ((lattice_[quo] >> rem) & 3))
            {
              lattice_[quo] ^= 2 << rem;
              lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
              mols_[i] = ((uint16_t*)&tars)[j];
              //mols_[i] = vdx;
            }
        }
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
        {
          lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
          lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint16_t vdx(((uint16_t*)&tars)[j]);
          if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
            {
              lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
              lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
              mols_[i] = vdx;
            }
        }
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
        {
          lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
          lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint16_t vdx(((uint16_t*)&tars)[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint16_t vdx(((uint16_t*)&tars)[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      union256 rans;
      for(unsigned j(0); j != 16; ++j)
        {
          rans.uint16[j] = rng_.IRan(0, ADJS-1);
        }
      __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rans.m256i));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint16_t vdx(((uint16_t*)&tars)[j]);
          if(lattice_[vdx] == vac_id_)
            {
              lattice_[vdx] = species_id_;
              lattice_[mols_[i]] = vac_id_;
              mols_[i] = vdx;
            }
        }
    }
  union256 rans;
  for(unsigned j(0); j != 16; ++j)
    {
      rans.uint16[j] = rng_.IRan(0, ADJS-1);
    }
  __m256i tars(compartment_.get_tars(*(__m256i*)(&mols_[i]), rans.m256i));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint16_t vdx(((uint16_t*)&tars)[j]);
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

/*
void Diffuser::walk()
{
  for(unsigned i(0), n(mols_.size()); i != n; ++i)
    { 
      const unsigned vdx(compartment_.get_tar(mols_[i], rng_.IRan(0, ADJS-1)));
      if(lattice_[vdx] == vac_id_)
        {
          lattice_[vdx] = species_id_;
          lattice_[mols_[i]] = vac_id_;
          mols_[i] = vdx;
        }
    }
}
*/

