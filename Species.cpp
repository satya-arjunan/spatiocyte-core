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

#include <Species.hpp>
#include <Compartment.hpp>
#include <Model.hpp>
#include <iostream>
#include <string>

Species::Species(const std::string name, const unsigned nmols, const double D,
                 Model& model, Compartment& comp, Species& vacant,
                 const bool is_comp_vacant):
  name_(get_init_name(name, comp, vacant, is_comp_vacant)),
  comp_(comp),
  vacant_(vacant),
  is_comp_vacant_(is_comp_vacant),
  lattice_(comp_.get_lattice()),
  id_(model.push_species(*this)),
  vac_id_(vacant_.get_id()),
  vac_xor_(vac_id_^id_),
  diffuser_(D, *this)
{
  std::cout << get_name_id() << std::endl;
  mols_.resize(nmols);
}

void Species::initialize()
{
  nbit_ = comp_.get_model().get_nbit();
  diffuser_.initialize();
}

void Species::populate()
{
  for(unsigned i(0), j(mols_.size()); i != j; ++i)
    {
      mol_t vdx(rng_.IntegerC(lattice_.size()*WORD/nbit_-1));
      while(lattice_[vdx*nbit_/WORD] & (1 << vdx*nbit_%WORD))
        {
          vdx = rng_.IntegerC(lattice_.size()*WORD/nbit_-1);
        }
      mols_[i] = vdx;
      lattice_[vdx*nbit_/WORD] ^= vac_xor_ << vdx*nbit_%WORD;
    }
}

bool Species::is_comp_vacant() const
{
  return is_comp_vacant_;
}

unsigned Species::get_id() const
{
  return id_;
}

unsigned Species::get_vac_id() const
{
  return vac_id_;
}

unsigned Species::get_vac_xor() const
{
  return vac_xor_;
}

Diffuser& Species::get_diffuser()
{
  return diffuser_;
}

Compartment& Species::get_comp()
{
  return comp_;
}

const std::string& Species::get_name() const
{
  return name_;
}

const std::string Species::get_name_id() const
{
  std::stringstream sid;
  sid << get_id();
  return std::string(get_name()+":"+sid.str());
  //return std::string(get_name()+":id:"+std::to_string(get_id())); c++11
}

const std::string Species::get_init_name(const std::string name,
                                           const Compartment& comp, 
                                           const Species& vacant,
                                           const bool is_comp_vacant) const
{
  if(is_comp_vacant)
    {
      return std::string(comp.get_name()+"/"+name);
    }
  return std::string(vacant.get_name()+"/"+name);
}

std::vector<mol_t>& Species::get_mols()
{
  return mols_;
}
