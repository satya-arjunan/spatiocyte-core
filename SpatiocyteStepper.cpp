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

void SpatiocyteStepper::setOffsets(
                     std::vector<std::vector<std::vector<int> > >& anOffsets)
{
  anOffsets.resize(2);
  anOffsets[0].resize(2);
  anOffsets[1].resize(2);
  //col=even, layer=even
  anOffsets.resize(ADJS*4);
  anOffsets[0][0].resize(12);
  anOffsets[0][0][0] = -1;
  anOffsets[0][0][1] = 1;
  anOffsets[0][0][2] = -nRows-1;
  anOffsets[0][0][3] = -nRows;
  anOffsets[0][0][4] = nRows-1;
  anOffsets[0][0][5] = nRows;
  anOffsets[0][0][6] = -nRows*nCols-nRows;
  anOffsets[0][0][7] = -nRows*nCols-1;
  anOffsets[0][0][8] = -nRows*nCols;
  anOffsets[0][0][9] = nRows*nCols-nRows;
  anOffsets[0][0][10] = nRows*nCols-1;
  anOffsets[0][0][11] = nRows*nCols;

  //col=odd, layer=even +12 = %col*12
  anOffsets[1][0].resize(12);
  anOffsets[1][0][0] = -1;
  anOffsets[1][0][1] = 1;
  anOffsets[1][0][2] = -nRows;
  anOffsets[1][0][3] = -nRows+1;
  anOffsets[1][0][4] = nRows;
  anOffsets[1][0][5] = nRows+1;
  anOffsets[1][0][6] = -nRows*nCols-nRows;
  anOffsets[1][0][7] = -nRows*nCols;
  anOffsets[1][0][8] = -nRows*nCols+1;
  anOffsets[1][0][9] = nRows*nCols-nRows;
  anOffsets[1][0][10] = nRows*nCols;
  anOffsets[1][0][11] = nRows*nCols+1;

  //col=even, layer=odd +24 = %layer*24
  anOffsets[0][1].resize(12);
  anOffsets[0][1][0] = -1;
  anOffsets[0][1][1] = 1;
  anOffsets[0][1][2] = -nRows;
  anOffsets[0][1][3] = -nRows+1;
  anOffsets[0][1][4] = nRows;
  anOffsets[0][1][5] = nRows+1;
  anOffsets[0][1][6] = -nRows*nCols;
  anOffsets[0][1][7] = -nRows*nCols+1;
  anOffsets[0][1][8] = -nRows*nCols+nRows;
  anOffsets[0][1][9] = nRows*nCols;
  anOffsets[0][1][10] = nRows*nCols+1;
  anOffsets[0][1][11] = nRows*nCols+nRows;

  //col=odd, layer=odd +36 = %col*12 + %layer*24
  anOffsets[1][1].resize(12);
  anOffsets[1][1][0] = -1;
  anOffsets[1][1][1] = 1;
  anOffsets[1][1][2] = -nRows;
  anOffsets[1][1][3] = -nRows+1;
  anOffsets[1][1][4] = nRows-1;
  anOffsets[1][1][5] = nRows;
  anOffsets[1][1][6] = -nRows*nCols-1;
  anOffsets[1][1][7] = -nRows*nCols+1;
  anOffsets[1][1][8] = -nRows*nCols+nRows;
  anOffsets[1][1][9] = nRows*nCols-1;
  anOffsets[1][1][10] = nRows*nCols;
  anOffsets[1][1][11] = nRows*nCols+nRows;

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

  //col=odd, layer=even +12 = %col*12
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

  //col=even, layer=odd +24 = %layer*24
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

  //col=odd, layer=odd +36 = %col*12 + %layer*24
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
      const unsigned aTar(getTar(theMols[i], theOffsets[
                                 (theMols[i]%(nRows*nCols)/nRows)%2][
                                 (theMols[i]/(nRows*nCols))%2][
                                 theRng.IntegerC(ADJS-1)]));
      /*
      const unsigned aTar(getTar(theMols[i], theOffsets[
                                 theRng.IntegerC(ADJS-1) + 
                                 (theMols[i]%(nRows*nCols)/nRows)%2*12 +
                                 (theMols[i]/(nRows*nCols))%2*24]));
                                 */
      if(!(theLattice[aTar/WORD] & (1 << aTar%WORD)))
        {
          theLattice[aTar/WORD] |= 1 << aTar%WORD;
          theLattice[theMols[i]/WORD] &= ~(1 << theMols[i]%WORD);
          theMols[i] = aTar;
        }
    }
}

