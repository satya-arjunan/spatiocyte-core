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

void Diffuser::walk()
{
  for(unsigned i(0); i != 100000; ++i)
  { 
    //const unsigned aTar(_comp.getTar2(_mols[i], _rng.IntegerC(ADJS-1)));
    const unsigned index(7);
    /*
    const unsigned aTar(_comp.getTar(_mols[i], index));
    std::cout << "aTar1:" << aTar << std::endl;
    */
    const unsigned aTar(_comp.getTar2(_mols[i], index));
    //std::cout << "aTar2:" << aTar2 << std::endl;
    if(!(_lattice[aTar/WORD] & (1 << aTar%WORD)))
      {
        _lattice[aTar/WORD] |= 1 << aTar%WORD;
        _lattice[_mols[i]/WORD] &= ~(1 << _mols[i]%WORD);
        _mols[i] = aTar;
      }
    /*
    const unsigned aTar(_comp.getTar(_mols[i], _rng.IntegerC(ADJS-1)));
    if(!_lattice[aTar])
      {
        _lattice[aTar] = i+1;
        _lattice[_mols[i]] = 0;
        _mols[i] = aTar;
      }
      */
  }
}

