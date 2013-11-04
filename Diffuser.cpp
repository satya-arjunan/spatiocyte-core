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

Diffuser::Diffuser(const double D, Species& species, Compartment& comp,
                   std::vector<unsigned>& mols):
  D_(D),
  species_(species),
  comp_(comp),
  mols_(mols),
  lattice_(comp.get_lattice()) {}

void Diffuser::walk()
{
  for(unsigned i(0), n(mols_.size()); i != n; ++i)
    { 
      const unsigned vdx(comp_.get_tar(mols_[i], rng_.IntegerC(ADJS-1)));
      if(!(lattice_[vdx/WORD] & (1 << vdx%WORD)))
        {
          lattice_[vdx/WORD] |= 1 << vdx%WORD;
          lattice_[mols_[i]/WORD] &= ~(1 << mols_[i]%WORD);
          mols_[i] = vdx;
        }
    }
}

