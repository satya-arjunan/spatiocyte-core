//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of the Spatiocyte package
//
//        Copyright (C) 2006-2009 Keio University
//        Copyright (C) 2010-2014 RIKEN
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

#include <cstring>
#include <Lattice.hpp>

Lattice::Lattice(const unsigned num_voxel,
    const Vector<unsigned>& dimensions, const unsigned num_box)
  : num_voxel_(num_voxel),
    num_box_(num_box),
    dimensions_(dimensions),
    box_dimensions_(dimensions) {}

void Lattice::initialize() {
  voxels_ = new voxel_t[num_voxel_];
  memset(voxels_, 0, sizeof(voxel_t)*num_voxel_);
  std::cout << "num_x:" << dimensions_.x << " num_y:" << dimensions_.y <<
    " num_z:" << dimensions_.z  << " num_voxel:" << num_voxel_ << " memory:"
    << num_voxel_*sizeof(voxel_t)/(1024*1024.0) << " MB" << std::endl;
}

unsigned Lattice::get_num_voxel() const {
  return num_voxel_;
}

const Vector<unsigned>& Lattice::get_dimensions() const {
  return dimensions_;
}

const Vector<unsigned>& Lattice::get_box_dimensions() const {
  return box_dimensions_;
}

voxel_t* Lattice::get_voxels() {
  return voxels_;
}
