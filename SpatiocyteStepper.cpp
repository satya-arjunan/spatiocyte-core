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

#include <SpatiocyteStepper.hpp>
#include <SpatiocyteCommon.hpp>

void SpatiocyteStepper::initialize()
{
  theLattice.resize(nVoxels/WORD+1, 0);
  theMols.resize(nMols);
  setOffsets(theOffsets);

  //populate lattice:
  for(unsigned short i(0); i != nMols; ++i)
    {
      unsigned aCoord(theRng.IntegerC(nVoxels-1));
      while(theLattice[aCoord/WORD] & (1 << aCoord%WORD))
        {
          aCoord = theRng.IntegerC(nVoxels-1);
        }
      theMols[i] = aCoord;
      theLattice[aCoord/WORD] |= 1 << aCoord%WORD;
    }
}

void SpatiocyteStepper::setOffsets(std::vector<int>& anOffsets)
{
  anOffsets.resize(ADJS);
  anOffsets[0] = -nRows*nCols;
  anOffsets[1] = -1;
  anOffsets[2] = -nRows;
  anOffsets[3] = +nRows;
  anOffsets[4] = 1;
  anOffsets[5] = nRows*nCols;
}

unsigned SpatiocyteStepper::getTar(const unsigned curr, const int anOffset)
{
  const long ret(curr+anOffset);
  if(ret < 0 || ret >= nVoxels)
    {
      return curr;
    }
  return ret;
}


void SpatiocyteStepper::step()
{
  for(unsigned short i(0); i != 10000; ++i)
    { 
      const unsigned aTar(getTar(theMols[i], 
                                 theOffsets[theRng.IntegerC(ADJS-1)]));
      if(!(theLattice[aTar/WORD] & (1 << aTar%WORD)))
        {
          theLattice[aTar/WORD] |= 1 << aTar%WORD;
          theLattice[theMols[i]/WORD] &= ~(1 << theMols[i]%WORD);
          theMols[i] = aTar;
        }
    }
}

