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
#include <Model.hpp>


Compartment::Compartment(std::string name, const double len_x,
                         const double len_y, const double len_z, Model& model)
  : name_(name),
    length_(len_x, len_y, len_z),
    center_(len_x/2, len_y/2, len_z/2),
    model_(model),
    volume_species_("volume", 0, 0, model, *this, volume_species_, true),
    surface_species_("surface", 0, 0, model, *this, surface_species_, true) {}

void Compartment::initialize() {
  nbit_ = model_.get_nbit();
  sur_xor_ = surface_species_.get_id()^volume_species_.get_id();
  lattice_.resize(ceil(double(NUM_VOXEL)*nbit_/WORD), 0);
  set_surface();
  std::cout << "nrow:" << NUM_ROW << " ncol:" << NUM_COL << " nlay:" <<
    NUM_LAY << " latticeSize:" << lattice_.size() << " memory:" << 
    lattice_.size()*sizeof(unsigned)/(1024*1024.0) << " MB" << std::endl;
}

unsigned Compartment::get_num_col() const {
  return NUM_COL;
}

unsigned Compartment::get_num_lay() const {
  return NUM_LAY;
}

unsigned Compartment::get_num_row() const {
  return NUM_ROW;
}

unsigned Compartment::get_num_voxel() const {
  return NUM_VOXEL;
}

/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  const bool odd_lay((vdx/NUM_COLROW)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-NUM_ROW-1 ;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx+NUM_ROW*(odd_lay-NUM_COL-1)-(odd_col&odd_lay);
    case 7:
      return vdx+!odd_col*(odd_lay-!odd_lay)-NUM_COLROW;
    case 8:
      return vdx+NUM_ROW*(odd_lay-NUM_COL)+(odd_col&!odd_lay);
    case 9:
      return vdx+NUM_ROW*(NUM_COL-!odd_lay)-(odd_col&odd_lay);
    case 10:
      return vdx+NUM_COLROW+!odd_col*(odd_lay-!odd_lay);
    case 11:
      return vdx+NUM_ROW*(NUM_COL+odd_lay)+(odd_col&!odd_lay);
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  unsigned value(0);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      value = -1;
    case 3:
      return vdx+(odd_col^odd_lay)-202+value;
    case 4:
      value = -1;
    case 5:
      return vdx+(odd_col^odd_lay)+202+value;
    case 6:
      value = 94132;
    case 7:
      return vdx+value-47066-(odd_col&odd_lay)-202*(!odd_lay);
    case 8:
      value = 94132;
    case 9:
      return vdx+value-47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 10:
      value = 94132;
    case 11:
      return vdx+value-47066+(odd_col&!odd_lay)+202*odd_lay;
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-nrowm_;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+nrowm_;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx-NUM_COLROW-(odd_col&odd_lay)-NUM_ROW*(!odd_lay);
    case 7:
      return vdx-NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-NUM_COLROW+(odd_col&!odd_lay)+NUM_ROW*odd_lay;
    case 9:
      return vdx+NUM_COLROW-(odd_col&odd_lay)+NUM_ROW*(!odd_lay);
    case 10:
      return vdx+NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+NUM_COLROW+(odd_col&!odd_lay)+NUM_ROW*odd_lay;
    }
  return vdx-1;
}
*/
/*
unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%47066/202)&1);
  const bool odd_lay((vdx/47066)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-201;
    case 3:
      return vdx+(odd_col^odd_lay)-202;
    case 4:
      return vdx+(odd_col^odd_lay)+201;
    case 5:
      return vdx+(odd_col^odd_lay)+202;
    case 6:
      return vdx-47066-(odd_col&odd_lay)-(202&(-!odd_lay));
    case 7:
      return vdx-47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-47066+(odd_col&!odd_lay)+(202&(-odd_lay));
    case 9:
      return vdx+47066-(odd_col&odd_lay)-(202&(-!odd_lay));
    case 10:
      return vdx+47066-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+47066+(odd_col&!odd_lay)+(202&(-odd_lay));
    }
  return vdx-1;
}
*/

unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const {
  const bool odd_col((vdx%NUM_COLROW/NUM_ROW)&1);
  const bool odd_lay((vdx/NUM_COLROW)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-NUM_ROW-1;
    case 3:
      return vdx+(odd_col^odd_lay)-NUM_ROW;
    case 4:
      return vdx+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
      return vdx+(odd_col^odd_lay)+NUM_ROW;
    case 6:
      return vdx-NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 7:
      return vdx-NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 8:
      return vdx-NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    case 9:
      return vdx+NUM_COLROW-(odd_col&odd_lay)-(NUM_ROW&(-!odd_lay));
    case 10:
      return vdx+NUM_COLROW-!(odd_col|odd_lay)+(!odd_col&odd_lay);
    case 11:
      return vdx+NUM_COLROW+(odd_col&!odd_lay)+(NUM_ROW&(-odd_lay));
    }
  return vdx-1;
}

const Vector& Compartment::get_center() const {
  return center_;
}

Species& Compartment::get_surface_species() {
  return surface_species_;
}

Species& Compartment::get_volume_species() {
  return volume_species_;
}

Model& Compartment::get_model() {
  return model_;
}

const std::string& Compartment::get_name() const {
  return name_;
}

std::vector<unsigned>& Compartment::get_lattice() {
  return lattice_;
}

void Compartment::set_surface() {
  //row_col xy-plane
  for (unsigned i(0); i != NUM_COLROW; ++i) {
      populate_mol(i);
      populate_mol(NUM_VOXEL-1-i);
    }
  for (unsigned i(1); i != NUM_LAY-1; ++i) {
      //layer_row yz-plane
      for (unsigned j(0); j != NUM_ROW; ++j) {
          populate_mol(i*NUM_COLROW+j);
          populate_mol(i*NUM_COLROW+j+NUM_ROW*(NUM_COL-1));
        }
      //layer_col xz-plane
      for (unsigned j(1); j != NUM_COL-1; ++j) {
          populate_mol(i*NUM_COLROW+j*NUM_ROW);
          populate_mol(i*NUM_COLROW+j*NUM_ROW+NUM_ROW-1);
        }
    }
  //std::cout << "surface size:" << mols.size() << " actual size:" <<
  //NUM_COLROW*2 + NUM_ROW*(NUM_LAY-2)*2 + (NUM_COL-2)*(NUM_LAY-2)*2 <<
  //std::endl;
}

void Compartment::populate_mol(const unsigned vdx) {
  lattice_[vdx*nbit_/WORD] ^= sur_xor_ << vdx*nbit_%WORD;
  surface_species_.get_mols().push_back(vdx);
}

