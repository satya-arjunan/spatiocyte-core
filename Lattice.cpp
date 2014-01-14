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

#include <cmath>
#include <cstring>
#include <Lattice.hpp>

Lattice::Lattice(const unsigned num_voxel,
    const Vector<unsigned>& dimensions, const unsigned num_box)
  : num_box_(num_box),
    box_dimensions_(Vector<unsigned>(unsigned(ceil(pow(num_box,1/3.0))),
        unsigned(pow(num_box,1/3.0)), unsigned(pow(num_box,1/3.0)))),
    box_voxel_dimensions_(dimensions),
    dimensions_(Vector<unsigned>(box_voxel_dimensions_.x*box_dimensions_.x,
        box_voxel_dimensions_.y*box_dimensions_.y,
        box_voxel_dimensions_.z*box_dimensions_.z)),
    num_box_voxel_(box_voxel_dimensions_.x*box_voxel_dimensions_.y*
        box_voxel_dimensions_.z) {}

void Lattice::initialize() {
  voxels_ = new voxel_t[get_num_voxel()];
  memset(get_voxels(), 0, sizeof(voxel_t)*get_num_voxel());
  box_voxels_ = new voxel_t*[get_num_box()];
  for(unsigned box(0); box != get_num_box(); ++box) {
    box_voxels_[box] = voxels_+box*get_num_box_voxel();
  }
  std::cout << "num_x:" << get_dimensions().x << " num_y:" << 
    get_dimensions().y << " num_z:" << get_dimensions().z  << " num_voxel:" <<
    get_num_voxel() << " memory:" << get_num_voxel()*sizeof(voxel_t)/
    (1024*1024.0) << " MB" << std::endl;
}

unsigned Lattice::get_num_box_voxel() const {
  return num_box_voxel_;
}

unsigned Lattice::get_num_box() const {
  return num_box_;
}

unsigned Lattice::get_num_voxel() const {
  return get_num_box_voxel()*get_num_box();
}

unsigned Lattice::box_mol_to_mol(const unsigned box, const umol_t box_mol)
    const {
  Vector<unsigned> mol_coord(box_mol_to_coord(box_mol));
  const Vector<unsigned> box_coord(box_to_coord(box));
  mol_coord.x += box_coord.x*box_voxel_dimensions_.x;
  mol_coord.y += box_coord.y*box_voxel_dimensions_.y;
  mol_coord.z += box_coord.z*box_voxel_dimensions_.z;
  return coord_to_mol(mol_coord);
}

unsigned Lattice::coord_to_mol(const Vector<unsigned>& coord) const {
  return coord.x*dimensions_.y+coord.y+coord.z*dimensions_.x*dimensions_.y;
}

Vector<unsigned> Lattice::box_mol_to_coord(const umol_t box_mol) const {
  const unsigned xy(box_voxel_dimensions_.x*box_voxel_dimensions_.y);
  return Vector<unsigned>(box_mol%(xy)/box_voxel_dimensions_.y, box_mol%(xy),
      box_mol/xy);
}

Vector<unsigned> Lattice::box_to_coord(const unsigned box) const {
  const unsigned xy(box_dimensions_.x*box_dimensions_.y);
  return Vector<unsigned>(box%(xy)/box_dimensions_.y, box%(xy), box/xy);
}

const Vector<unsigned>& Lattice::get_dimensions() const {
  return dimensions_;
}

const Vector<unsigned>& Lattice::get_box_dimensions() const {
  return box_dimensions_;
}

const Vector<unsigned>& Lattice::get_box_voxel_dimensions() const {
  return box_voxel_dimensions_;
}

voxel_t* Lattice::get_voxels() {
  return voxels_;
}

voxel_t** Lattice::get_box_voxels() {
  return box_voxels_;
}
