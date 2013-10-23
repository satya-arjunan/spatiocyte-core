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
#include <climits>

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
  /*
  anOffsets.resize(ADJS);
  anOffsets[0] = -nRows*nCols;
  anOffsets[1] = -1;
  anOffsets[2] = -nRows;
  anOffsets[3] = +nRows;
  anOffsets[4] = 1;
  anOffsets[5] = nRows*nCols;
  */
  //col=even, layer=even
  //col=0, layer=0
  /*
  anOffsets.resize(ADJS*4);
  anOffsets[0] = -1;
  anOffsets[1] = 1;
  anOffsets[2] = -nRows-1+col*layer;
  anOffsets[3] = -nRows+col*layer;
  anOffsets[4] = nRows-1+(col-layer)*(col-layer);
  anOffsets[5] = nRows+(col-layer)*(col-layer);
  anOffsets[6] = -nRows*nCols-nRows+nRows*layer-col*layer;
  anOffsets[7] = -nRows*nCols+layer+(1-layer)*(col-1); //-1-layer*(col-2)+col
  anOffsets[8] = -nRows*nCols+nRows*layer+(1-layer)*col;
  anOffsets[9] = nRows*nCols-nRows+layer*nRows-col*layer;
  anOffsets[10] = nRows*nCols+(1-col)*(2*layer-1);
  anOffsets[11] = nRows*nCols+nRows*layer+(1-layer)*col;


  /*
  anOffsets.resize(ADJS*4);
  anOffsets[0] = -1;
  anOffsets[1] = 1;
  anOffsets[2] = -nRows-1;
  anOffsets[3] = -nRows;
  anOffsets[4] = nRows-1;
  anOffsets[5] = nRows;
  anOffsets[6] = -nRows*nCols-nRows;
  anOffsets[7] = -nRows*nCols-1;
  anOffsets[8] = -nRows*nCols;
  anOffsets[9] = nRows*nCols-nRows;
  anOffsets[10] = nRows*nCols-1;
  anOffsets[11] = nRows*nCols;

  //col=even, layer=odd +24 = %layer*24
  //col=0, layer=1
  anOffsets[24] = -1;
  anOffsets[25] = 1;
  anOffsets[26] = -nRows;
  anOffsets[27] = -nRows+1;
  anOffsets[28] = nRows;
  anOffsets[29] = nRows+1;
  anOffsets[30] = -nRows*nCols;
  anOffsets[31] = -nRows*nCols+1;
  anOffsets[32] = -nRows*nCols+nRows;
  anOffsets[33] = nRows*nCols;
  anOffsets[34] = nRows*nCols+1;
  anOffsets[35] = nRows*nCols+nRows;

  //col=odd, layer=even +12 = %col*12
  //col=1, layer=0
  anOffsets[12] = -1;
  anOffsets[13] = 1;
  anOffsets[14] = -nRows;
  anOffsets[15] = -nRows+1;
  anOffsets[16] = nRows;
  anOffsets[17] = nRows+1;
  anOffsets[18] = -nRows*nCols-nRows;
  anOffsets[19] = -nRows*nCols;
  anOffsets[20] = -nRows*nCols+1;
  anOffsets[21] = nRows*nCols-nRows;
  anOffsets[22] = nRows*nCols;
  anOffsets[23] = nRows*nCols+1;


  //col=odd, layer=odd +36 = %col*12 + %layer*24
  //col=1, layer=1
  anOffsets[36] = -1;
  anOffsets[37] = 1;
  anOffsets[38] = -nRows;
  anOffsets[39] = -nRows+1;
  anOffsets[40] = nRows-1;
  anOffsets[41] = nRows;
  anOffsets[42] = -nRows*nCols-1;
  anOffsets[43] = -nRows*nCols+1;
  anOffsets[44] = -nRows*nCols+nRows;
  anOffsets[45] = nRows*nCols-1;
  anOffsets[46] = nRows*nCols;
  anOffsets[47] = nRows*nCols+nRows;
  */
}

unsigned SpatiocyteStepper::getTar(const unsigned curr, const unsigned aRand)
{
  const unsigned col((curr%(nRows*nCols)/nRows)%2);
  const unsigned layer((curr/(nRows*nCols))%2);
  int ret(0);
  switch(aRand)
    {
    case 0:
      ret = -1;
      break;
    case 1:
      ret = 1;
      break;
    case 2:
      ret = -nRows-1+col*layer;
      break;
    case 3:
      ret = -nRows+col*layer;
      break;
    case 4:
      ret = nRows-1+(col-layer)*(col-layer);
      break;
    case 5:
      ret = nRows+(col-layer)*(col-layer);
      break;
    case 6:
      ret = -nRows*nCols-nRows+nRows*layer-col*layer;
      break;
    case 7:
      ret = -nRows*nCols+layer+(1-layer)*(col-1); //-1-layer*(col-2)+col
      break;
    case 8:
      ret = -nRows*nCols+nRows*layer+(1-layer)*col;
      break;
    case 9:
      ret = nRows*nCols-nRows+layer*nRows-col*layer;
      break;
    case 10:
      ret = nRows*nCols+(1-col)*(2*layer-1);
      break;
    case 11:
      ret = nRows*nCols+nRows*layer+(1-layer)*col;
      break;
    }
  if(ret < curr || ret >= nVoxels-curr)
    {
      return curr;
    }
  return ret+curr;
}


void SpatiocyteStepper::step()
{
  for(unsigned i(0); i != 10000; ++i)
    { 
      const unsigned aTar(getTar(theMols[i], theRng.IntegerC(ADJS-1)));
      if(!(theLattice[aTar/WORD] & (1 << aTar%WORD)))
        {
          theLattice[aTar/WORD] |= 1 << aTar%WORD;
          theLattice[theMols[i]/WORD] &= ~(1 << theMols[i]%WORD);
          theMols[i] = aTar;
        }
    }
}

