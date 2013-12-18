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

Diffuser::Diffuser(const double D, Species& species):
  D_(D),
  species_(species),
  comp_(species_.get_comp()),
  mols_(species_.get_mols()),
  lattice_(comp_.get_lattice()),
  vac_id_(species_.get_vac_id()),
  vac_xor_(species_.get_vac_xor()) {}

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
  nboxdiv_ = species_.get_nboxdiv();
  for(unsigned i(0), m(mols_.size()); i < m; ++i)
    { 
      for(unsigned j(0), n(mols_[i].size()); j < n; ++j)
        { 
          const unsigned vdx(comp_.get_tar(mols_[i][j], rng_.IntegerC(ADJS-1)));
          if(0 == ((lattice_[vdx*2/WORD] >> vdx*2%WORD) & 3))
            {
              lattice_[vdx*2/WORD] ^= 2 << vdx*2%WORD;
              lattice_[mols_[i][j]*2/WORD] ^= 2 << mols_[i][j]*2%WORD;
              if(i == vdx/nboxdiv_)
                {
                  mols_[i][j] = vdx;
                }
              else
                {
                  mols_[i][j] = mols_[i].back();
                  mols_[i].pop_back();
                  mols_[vdx/nboxdiv_].push_back(vdx);
                  n--;
                }
            }
        }
    }
}

