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

#include <Spatiocyte.hpp>
#include <Common.hpp>
#include <Species.hpp>

#define HCP_X double(VOXEL_RADIUS*1.7320508075688772)
#define HCP_Z double(VOXEL_RADIUS*1.632993161855452)
#define NUM_COL umol_t(LENGTH_X/HCP_X+3)
#define NUM_LAY umol_t(LENGTH_Z/HCP_Z+3)
//#define NUM_ROW unsigned(LENGTH_Y/VOXEL_RADIUS/2+3) //correct version
#define NUM_ROW umol_t(LENGTH_Y/VOXEL_RADIUS/2+2)
#define NUM_COLROW umol_t(NUM_COL*NUM_ROW)
#define NUM_COLROWROW umol_t(NUM_COLROW*NUM_ROW)
#define NUM_VOXEL umol_t(NUM_COLROW*NUM_LAY)

class Compartment { 
 public: 
  Compartment(std::string, const double, const double, const double, Model&);
  ~Compartment() {}
  void initialize();
  __m256i set_tars(const __m256i, __m256i) const;
  umol_t get_num_col() const;
  umol_t get_num_lay() const;
  umol_t get_num_row() const;
  umol_t get_num_voxel() const;
  umol_t get_lattice_size() const;
  umol_t get_tar(const umol_t, const unsigned) const;
  const Vector& get_center() const;
  Species& get_surface_species();
  Species& get_volume_species();
  Model& get_model();
  const std::string& get_name() const;
  int* get_lattice();
 private:
  void set_surface();
  void populate_mol(const umol_t);
  void setOffsets();
 private:
  const std::string name_;
  const Vector length_;
  const Vector center_;
  umol_t lattice_size_;
  int* lattice_;
  Model& model_;
  Species volume_species_;
  Species surface_species_;
  voxel_t nbit_;
  voxel_t sur_xor_;
  mol2_t* offsets_;
  const __m256i m256i_1_;
  const __m256i m256i_m1_;
  const __m256i m256i_12_;
  const __m256i m256i_24_;
  const __m256i m256i_24_12_;
  const __m256i m256i_num_colrow_;
  __m256i multiplier_colrow_;
  __m256i multiplier_row_;
  __m128i nshift_colrow_;
  __m128i nshift_row_;
};

#endif /* __Compartment_hpp */

