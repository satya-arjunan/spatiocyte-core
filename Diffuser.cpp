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

#include <Diffuser.hpp>
#include <Compartment.hpp>
#include <Model.hpp>
#include <time.h>

Diffuser::Diffuser(const double D, Species& species):
  D_(D),
  species_(species),
  comp_(species_.get_comp()),
  mols_(species_.get_mols()),
  species_id_(species_.get_id()),
  vac_id_(species_.get_vac_id()),
  vac_xor_(species_.get_vac_xor()),
  rng_(time(0)) {}

void Diffuser::initialize()
{
  lattice_ = comp_.get_lattice();
  nbit_ = species_.get_comp().get_model().get_nbit();
  one_nbit_ = pow(2, nbit_)-1;
  std::cout << species_.get_name_id() << std::endl;
  std::cout << "nbit:" << nbit_ << std::endl;
  std::cout << "one_nbit:" << one_nbit_ << std::endl;
  std::cout << "vac_id:" << vac_id_ << std::endl;
  std::cout << "vac_xor:" << vac_xor_ << std::endl;
}

void Diffuser::walk()
{
  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  uint32_t tars[16];
  for(unsigned k(0); k != n; ++k)
    {
      __m256i mols(_mm256_load_si256((__m256i*)(&mols_[i])));
      comp_.set_tars(mols, rng_.Ran16(), tars);
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
  comp_.set_tars(mols, rng_.Ran16(), tars);
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
      comp_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
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
  comp_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
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
      comp_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
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
  comp_.set_tars(*(__m256i*)(&mols_[i]), rng_.Ran16(), tars);
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
      __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
  __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
      __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
  __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
      __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
  __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
      __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
  __m256i tars(comp_.get_tars(*(__m256i*)(&mols_[i]), rng_.Ran16()));
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
