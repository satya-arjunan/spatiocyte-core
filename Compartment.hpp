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


#ifndef __Compartment_hpp
#define __Compartment_hpp

#include <Common.hpp>

class Compartment
{ 
public: 
  Compartment(const double voxRadius, const double lenX,
              const double lenY, const double lenZ);
  ~Compartment() {}
  void populate();
  std::vector<unsigned>& getLattice()
    {
      return _lattice;
    }
  int getCols()
    {
      return _cols;
    }
  int getLays()
    {
      return _lays;
    }
  int getRows()
    {
      return _rows;
    }
  unsigned getVoxs() const
    {
      return _voxs;
    }
  double getVoxRadius() const
    {
      return _voxRadius;
    }
  Vector getCenter() const
    {
      return _center;
    }
  unsigned getTar(const unsigned, const unsigned) const;
  unsigned getTar2(const unsigned, const unsigned) const;
private:
  void setOffsets();
private:
  const double _hcpX;
  const double _hcpO;
  const double _hcpZ;
  const int _cols;
  const int _lays;
  const int _rows;
  const unsigned _ecol;
  const unsigned _erow;
  const unsigned _layVoxs;
  const unsigned _voxs;
  const double _lenX;
  const double _lenY;
  const double _lenZ;
  const double _voxRadius;
  const Vector _center;
  std::vector<unsigned> _lattice;
  std::vector<int> _offsets;
};

#endif /* __Compartment_hpp */

