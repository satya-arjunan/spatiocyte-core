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


#ifndef __SpatiocyteStepper_hpp
#define __SpatiocyteStepper_hpp

#include <RandomLib/Random.hpp>

class SpatiocyteStepper
{ 
public: 
  SpatiocyteStepper():
    nMols(10000),
    nDim(pow(11722464,1.0/3)),
    nCols(nDim),
    nLayers(nDim),
    nRows(nDim),
    nVoxels(nCols*nLayers*nRows) {}
  virtual ~SpatiocyteStepper() {}
  virtual void initialize();
  virtual void step();
private:
  void setOffsets(std::vector<int>&);
  unsigned getTar(const unsigned, const unsigned);
private:
  const short nMols;
  const unsigned nDim;
  const int nCols;
  const int nLayers;
  const int nRows;
  const unsigned nVoxels;
  std::vector<int> theOffsets;
  RandomLib::Random theRng;
  std::vector<unsigned> theLattice;
  std::vector<unsigned> theMols;
};

#endif /* __SpatiocyteStepper_hpp */

