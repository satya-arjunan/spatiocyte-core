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


Compartment::Compartment(const double vox_radius, const double len_x,
                         const double len_y, const double len_z,
                         Model& model):
  hcpx_(vox_radius*sqrt(3)),
  hcpo_(vox_radius/sqrt(3)), //protruding length_x at an odd numbered layer
  hcpz_(vox_radius*sqrt(8.0/3)),
  vox_radius_(vox_radius),
  ncol_(rint(len_x/hcpx_)+2),
  nlay_(rint(len_z/hcpz_)+2),
  nrow_(rint(len_y/vox_radius_/2)+2),
  ncolrow_(ncol_*nrow_),
  nvox_(ncolrow_*nlay_),
  length_(len_x, len_y, len_z),
  center_(len_x/2, len_y/2, len_z/2),
  model_(model),
  volume_(0, 0, model, *this, true),
  surface_(0, 0, model, *this, true) {}

void Compartment::initialize()
{
  nbit_ = model_.get_nbit();
  sur_xor_ = surface_.get_id()^volume_.get_id();
  lattice_.resize(ceil(double(nvox_)*nbit_/WORD), 0);
  set_surface();
  std::cout << "nrow:" << nrow_ << " ncol:" << ncol_ << " nlay:" << nlay_ <<
    std::endl;
}

unsigned Compartment::get_ncol() const
{
  return ncol_;
}

unsigned Compartment::get_nlay() const
{
  return nlay_;
}

unsigned Compartment::get_nrow() const
{
  return nrow_;
}

unsigned Compartment::get_nvox() const
{
  return nvox_;
}

unsigned Compartment::get_tar(const unsigned vdx, const unsigned nrand) const
{
  const bool odd_col((vdx%ncolrow_/nrow_)&1);
  const bool odd_lay((vdx/ncolrow_)&1);
  switch(nrand)
    {
    case 1:
      return vdx+1;
    case 2:
      return vdx+(odd_col^odd_lay)-nrow_-1 ;
    case 3:
      return vdx+(odd_col^odd_lay)-nrow_;
    case 4:
      return vdx+(odd_col^odd_lay)+nrow_-1;
    case 5:
      return vdx+(odd_col^odd_lay)+nrow_;
    case 6:
      return vdx+nrow_*(odd_lay-ncol_-1)-(odd_col&odd_lay);
    case 7:
      return vdx+!odd_col*(odd_lay-!odd_lay)-ncolrow_;
    case 8:
      return vdx+nrow_*(odd_lay-ncol_)+(odd_col&!odd_lay);
    case 9:
      return vdx+nrow_*(ncol_-!odd_lay)-(odd_col&odd_lay);
    case 10:
      return vdx+ncolrow_+!odd_col*(odd_lay-!odd_lay);
    case 11:
      return vdx+nrow_*(ncol_+odd_lay)+(odd_col&!odd_lay);
    }
  return vdx-1;
}

double Compartment::get_vox_radius() const
{
  return vox_radius_;
}

const Vector& Compartment::get_center() const
{
  return center_;
}

Species& Compartment::get_surface()
{
  return surface_;
}

std::vector<unsigned>& Compartment::get_lattice()
{
  return lattice_;
}

void Compartment::set_surface()
{
  //row_col xy-plane
  for(unsigned i(0); i != ncolrow_; ++i)
    {
      populate_mol(i);
      populate_mol(nvox_-1-i);
    }
  for(unsigned i(1); i != nlay_-1; ++i)
    {
      //layer_row yz-plane
      for(unsigned j(0); j != nrow_; ++j)
        {
          populate_mol(i*ncolrow_+j);
          populate_mol(i*ncolrow_+j+nrow_*(ncol_-1));
        }
      //layer_col xz-plane
      for(unsigned j(1); j != ncol_-1; ++j)
        {
          populate_mol(i*ncolrow_+j*nrow_);
          populate_mol(i*ncolrow_+j*nrow_+nrow_-1);
        }
    }
  //std::cout << "surface size:" << mols.size() << " actual size:" <<
  //ncolrow_*2 + nrow_*(nlay_-2)*2 + (ncol_-2)*(nlay_-2)*2 << std::endl;
}

void Compartment::populate_mol(const unsigned vdx)
{
  lattice_[vdx*nbit_/WORD] ^= sur_xor_ << vdx*nbit_%WORD;
  surface_.get_mols().push_back(vdx);
}

