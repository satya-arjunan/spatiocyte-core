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

#include <Compartment.hpp>
#include <math.h>

Compartment::Compartment(const double voxRadius, const double lenX,
                         const double lenY, const double lenZ):
  _hcpX(voxRadius*sqrt(3)),
  _hcpO(voxRadius/sqrt(3)), //protruding lenX at an odd numbered layer
  _hcpZ(voxRadius*sqrt(8.0/3)),
  /*
  _cols((int)rint(lenX/_hcpX)),
  _lays((int)rint(lenZ/_hcpZ)),
  _rows((int)rint(lenY/voxRadius/2)),
  */
  _cols(8),
  _lays(10),
  _rows(5),
  _voxs(_cols*_lays*_rows),
  _lenX(lenX),
  _lenY(lenY),
  _lenZ(lenZ),
  _voxRadius(voxRadius),
  _center(lenX/2, lenY/2, lenZ/2)
{
  _lattice.resize(ceil(double(_voxs)/WORD), 0);
  setOffsets();
  //_lattice.resize(_voxs, 0);
}

unsigned Compartment::getTar2(const unsigned curr, const unsigned aRand) const
{
  const bool col((curr%(_rows*_cols)/_rows)&1);
  const bool layer((curr/(_rows*_cols))&1);
  int ret(0);
  int x,y;
  switch(aRand)
    {
    case 0:
      ret = -1;
      break;
    case 1:
      ret = 1;
      break;
    case 2:
      //ret = -_rows - 1 + (col-layer)*(col-layer);
      ret = (col^layer) -_rows - 1 ;
      break;
    case 3:
      //ret = -_rows + (col-layer)*(col-layer);
      ret = (col^layer) -_rows;
      break;
    case 4:
      //ret = _rows - 1 + (col-layer)*(col-layer);
      ret = (col^layer) + _rows - 1;
      break;
    case 5:
      //ret = _rows + (col-layer)*(col-layer);
      ret = (col^layer) + _rows;
      break;
    case 6:
      //ret = _rows*(layer-_cols-1) - col*layer;
      ret = _rows*(layer-_cols-1) - (col&layer);
      //ret = -_rows*(_cols+!layer) - (col&layer);
      break;
    case 7:
      //ret = -_rows*_cols - 1 + 2*layer*(1-col) + col;
      //ret = -_rows*_cols-!(col|layer)+(!col&layer);
      ret = -_rows*_cols+!col*(layer-!layer);
      break;
    case 8:
      //ret = -_rows*_cols + layer*(_rows-col) + col;
      ret = _rows*(layer-_cols)+(col&!layer);
      break;
    case 9:
      //ret = _rows*(_cols-1+layer) - col*layer;
      //ret = _rows*(_cols-1+layer) - (col&layer);
      ret = _rows*(_cols-!layer)-(col&layer);
      break;
    case 10:
      //ret = _rows*_cols + (1-col)*(2*layer-1);
      ret = _rows*_cols+!col*(layer-!layer);
      break;
    case 11:
      //ret = _rows*(_cols+layer) + (1-layer)*col;
      ret = _rows*(_cols+layer)+(col&!layer);
      break;
    }
  /*
  int ret(0);
  switch(aRand)
    {
    case 0:
      ret = -_rows*_cols;
      break;
    case 1:
      ret = -1;
      break;
    case 2:
      ret = -_rows;
      break;
    case 3:
      ret = _rows;
      break;
    case 4:
      ret = 1;
      break;
    case 5:
      ret = _rows*_cols;
      break;
    }
    */

  if(long(curr)+ret < 0 || ret+long(curr) >= _voxs)
    {
      return curr;
    }
  return ret+curr;
}

unsigned Compartment::getTar(const unsigned curr, const unsigned aRand) const
{
  const unsigned col((curr%(_rows*_cols)/_rows)%2);
  const unsigned layer((curr/(_rows*_cols))%2);
  const unsigned index(aRand+layer*24+col*12);
  const long ret(curr+_offsets[index]);
  if(ret < 0 || ret >= _voxs)
    {
      return curr;
    }
  return ret;
}

void Compartment::setOffsets()
{
  /*
  _offsets.resize(ADJS);
  _offsets[0] = -_rows*_cols;
  _offsets[1] = -1;
  _offsets[2] = -_rows;
  _offsets[3] = +_rows;
  _offsets[4] = 1;
  _offsets[5] = _rows*_cols;
  */


  /*
      ret = _rows*_cols+_rows*layer+(col&!layer)
      ret = _rows*(_cols+layer)+(col&!layer)


  00 = 0
  _offsets[11] = _rows*_cols;
  01 = 1
  _offsets[35] = _rows*_cols+_rows;
  10 = 1
  _offsets[23] = _rows*_cols+1;
  11 = 2
  _offsets[47] = _rows*_cols+_rows;
  */


  //col=even, layer=even
  _offsets.resize(ADJS*4);
  _offsets[0] = -1;
  _offsets[1] = 1;
  _offsets[2] = -_rows-1;
  _offsets[3] = -_rows;
  _offsets[4] = _rows-1;
  _offsets[5] = _rows;
  _offsets[6] = -_rows*_cols-_rows;
  _offsets[7] = -_rows*_cols-1;
  _offsets[8] = -_rows*_cols;
  _offsets[9] = _rows*_cols-_rows;
  _offsets[10] = _rows*_cols-1;
  _offsets[11] = _rows*_cols;


  //col=even, layer=odd +24 = %layer*24
  _offsets[24] = -1;
  _offsets[25] = 1;
  _offsets[26] = -_rows;
  _offsets[27] = -_rows+1;
  _offsets[28] = _rows;
  _offsets[29] = _rows+1;
  _offsets[30] = -_rows*_cols;
  _offsets[31] = -_rows*_cols+1;
  _offsets[32] = -_rows*_cols+_rows;
  _offsets[33] = _rows*_cols;
  _offsets[34] = _rows*_cols+1;
  _offsets[35] = _rows*_cols+_rows;

  //col=odd, layer=even +12 = %col*12
  _offsets[12] = -1;
  _offsets[13] = 1;
  _offsets[14] = -_rows;
  _offsets[15] = -_rows+1;
  _offsets[16] = _rows;
  _offsets[17] = _rows+1;
  _offsets[18] = -_rows*_cols-_rows;
  _offsets[19] = -_rows*_cols;
  _offsets[20] = -_rows*_cols+1;
  _offsets[21] = _rows*_cols-_rows;
  _offsets[22] = _rows*_cols;
  _offsets[23] = _rows*_cols+1;


  //col=odd, layer=odd +36 = %col*12 + %layer*24
  _offsets[36] = -1;
  _offsets[37] = 1;
  _offsets[38] = -_rows-1;
  _offsets[39] = -_rows;
  _offsets[40] = _rows-1;
  _offsets[41] = _rows;
  _offsets[42] = -_rows*_cols-1;
  _offsets[43] = -_rows*_cols; //a
  _offsets[44] = -_rows*_cols+_rows;
  _offsets[45] = _rows*_cols-1;
  _offsets[46] = _rows*_cols;
  _offsets[47] = _rows*_cols+_rows;
}
