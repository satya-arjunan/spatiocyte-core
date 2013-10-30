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


#ifndef __Species_hpp
#define __Species_hpp

#include <RandomLib/Random.hpp>
#include <Common.hpp>
#include <Diffuser.hpp>
#include <Compartment.hpp>

class Species
{ 
public: 
  Species(const unsigned nmols, const double D, Compartment& comp):
    _isCompVacant(false),
    _comp(comp),
    _diffuser(D, *this, comp, _mols),
    _lattice(comp.getLattice())
  {
    _mols.resize(nmols);
  }
  ~Species() {}
  std::vector<unsigned>& getMols()
    {
      return _mols;
    }
  void populate()
    {
      for(unsigned short i(0),  j(_mols.size()); i != j; ++i)
        {
          unsigned coord(_rng.IntegerC(_lattice.size()*WORD-1));
          while(_lattice[coord/WORD] & (1 << coord%WORD))
            {
              coord = _rng.IntegerC(_lattice.size()-1);
            }
          _mols[i] = coord;
          _lattice[coord/WORD] |= 1 << coord%WORD;
        }
    }
  Diffuser& getDiffuser()
    {
      return _diffuser;
    }
  bool getIsCompVacant() const
    {
      return _isCompVacant;
    }
private:
  bool _isCompVacant;
  RandomLib::Random _rng;
  Compartment& _comp;
  Diffuser _diffuser;
  std::vector<unsigned> _mols;
  std::vector<unsigned>& _lattice;
};

#endif /* __Species_hpp */

