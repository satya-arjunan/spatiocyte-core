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
  _cols(227),
  _lays(227),
  _rows(227),
  _voxs(_cols*_lays*_rows),
  _lenX(lenX),
  _lenY(lenY),
  _lenZ(lenZ),
  _voxRadius(voxRadius),
  _center(lenX/2, lenY/2, lenZ/2)
{
  _lattice.resize(_voxs/WORD, 0);
}

unsigned Compartment::getTar(const unsigned curr, const unsigned aRand) const
{
  const unsigned col((curr%(_rows*_cols)/_rows)%2);
  const unsigned layer((curr/(_rows*_cols))%2);
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
      ret = -_rows-1+col*layer;
      break;
    case 3:
      ret = -_rows+col*layer;
      break;
    case 4:
      ret = _rows-1+(col-layer)*(col-layer);
      break;
    case 5:
      ret = _rows+(col-layer)*(col-layer);
      break;
    case 6:
      ret = -_rows*_cols-_rows+_rows*layer-col*layer;
      break;
    case 7:
      ret = -_rows*_cols+layer+(1-layer)*(col-1); //-1-layer*(col-2)+col
      break;
    case 8:
      ret = -_rows*_cols+_rows*layer+(1-layer)*col;
      break;
    case 9:
      ret = _rows*_cols-_rows+layer*_rows-col*layer;
      break;
    case 10:
      ret = _rows*_cols+(1-col)*(2*layer-1);
      break;
    case 11:
      ret = _rows*_cols+_rows*layer+(1-layer)*col;
      break;
    }
  if(ret < curr || ret >= _voxs-curr)
    {
      return curr;
    }
  return ret+curr;
}



