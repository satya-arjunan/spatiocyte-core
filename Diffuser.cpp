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
  lattice_(comp_.get_lattice()),
  vac_id_(species_.get_vac_id()),
  vac_xor_(species_.get_vac_xor()),
  rng_(time(0)) {}

void Diffuser::initialize()
{
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
  /*
  union256 ran(rng_.Ran16());
  for(unsigned i(0); i != 16; ++i)
    {
      std::cout << "ran:" << ran.uint16[i] << std::endl;
    }
  comp_.set_tars((__m256i*)(&mols_[0]), ran);
  exit(0);
  */
  //cast uint16_t to uint32_t
  //idx.m256i = _mm256_cvtepu16_epi32(ran.m128i[0]);
  //res = idx;
  /*
  for(unsigned i(0); i != 8; ++i)
    {
      std::cout << "i:" << i << " " << ran.uint16[i] << " " << 
        idx.uint32[i] << " " << res.int32[i] << std::endl;
    } 
    */
  

  const unsigned n(mols_.size()/16);
  const unsigned m(mols_.size()%16);
  unsigned i(0);
  for(unsigned k(0); k != n; ++k)
    {
      union512 tars(comp_.get_tars((__m256i*)(&mols_[i]), rng_.Ran16()));
      for(unsigned j(0); j != 16; ++j, ++i)
        {
          const uint32_t vdx(tars.uint32[j]);
          if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
            {
              lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
              lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
              mols_[i] = vdx;
            }
        }
    }
  union512 tars(comp_.get_tars((__m256i*)(&mols_[i]), rng_.Ran16()));
  for(unsigned j(0); j != m; ++j, ++i)
    {
      const uint32_t vdx(tars.uint32[j]);
      if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
        {
          lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
          lattice_[mols_[i]*2/WORD] ^= 2 << mols_[i]*2%WORD;
          mols_[i] = vdx;
        }
    }
}

