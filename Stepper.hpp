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


#ifndef __Stepper_hpp
#define __Stepper_hpp

#include <Compartment.hpp>

class Stepper
{ 
public: 
  Stepper(Compartment& comp):
    _lattice(comp.getLattice()),
    _mols(comp.getMols()),
    _comp(comp),
    _cols(comp.getCols()),
    _lays(comp.getLays()),
    _rows(comp.getRows()),
    _voxs(comp.getVoxs()) {}
  virtual ~Stepper() {}
  virtual void step();
  unsigned getTar(const unsigned, const unsigned) const;
private:
  std::vector<unsigned>& _lattice;
  std::vector<unsigned>& _mols;
  RandomLib::Random _rng;
  Compartment& _comp;
  const int _cols;
  const int _lays;
  const int _rows;
  const unsigned _voxs;
};

#endif /* __Stepper_hpp */

