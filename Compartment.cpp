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

#include <math.h>
#include <Compartment.hpp>
#include <Species.hpp>


Compartment::Compartment(const double voxRadius, const double lenX,
                         const double lenY, const double lenZ):
  _hcpX(voxRadius*sqrt(3)),
  _hcpO(voxRadius/sqrt(3)), //protruding lenX at an odd numbered lay
  _hcpZ(voxRadius*sqrt(8.0/3)),
  _cols(rint(lenX/_hcpX)),
  _lays(rint(lenZ/_hcpZ)),
  _rows(rint(lenY/voxRadius/2)),
  /*
  _cols(227),
  _lays(227),
  _rows(227),
  */
  _ecol(_cols-1),
  _erow(_rows-1),
  _layVoxs(_cols*_rows),
  _voxs(_cols*_lays*_rows),
  _lenX(lenX),
  _lenY(lenY),
  _lenZ(lenZ),
  _voxRadius(voxRadius),
  _center(lenX/2, lenY/2, lenZ/2),
  _boundary(0, 0, *this)
{
  std::cout << "rows:" << _rows << " cols:" << _cols << " lays:" <<
    _lays << std::endl;
  _lattice.resize(ceil(double(_voxs)/WORD), 0);
  setBoundary();
  //_lattice.resize(_voxs, 0);
}

void Compartment::setBoundary()
{
  std::vector<unsigned>& mols(_boundary.getMols());
  for(unsigned i(0); i != _layVoxs; ++i)
    {
      unsigned coord(i);
      _lattice[coord/WORD] |= 1 << coord%WORD;
      mols.push_back(coord);
    }
}

unsigned Compartment::getTar(const unsigned curr, const unsigned aRand) const
{
  const bool col((curr%_layVoxs/_rows)&1);
  const bool lay((curr/_layVoxs)&1);
  int ret(-1);
  switch(aRand)
    {
    case 1:
      ret = 1;
      break;
    case 2:
      ret = (col^lay)-_rows-1 ;
      break;
    case 3:
      ret = (col^lay)-_rows;
      break;
    case 4:
      ret = (col^lay)+_rows-1;
      break;
    case 5:
      ret = (col^lay)+_rows;
      break;
    case 6:
      ret = _rows*(lay-_cols-1)-(col&lay);
      break;
    case 7:
      ret = !col*(lay-!lay)-_rows*_cols;
      break;
    case 8:
      ret = _rows*(lay-_cols)+(col&!lay);
      break;
    case 9:
      ret = _rows*(_cols-!lay)-(col&lay);
      break;
    case 10:
      ret = _rows*_cols+!col*(lay-!lay);
      break;
    case 11:
      ret = _rows*(_cols+lay)+(col&!lay);
      break;
    }
  if(long(curr)+ret < 0 || ret+long(curr) >= _voxs)
    {
      return curr;
    }
  return ret+curr;
}


