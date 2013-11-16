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
#include <Species.hpp>

class Compartment
{ 
public: 
  Compartment(std::string, const double, const double, const double,
              const double, Model&);
  ~Compartment() {}
  void initialize();
  unsigned get_ncol() const;
  unsigned get_nlay() const;
  unsigned get_nrow() const;
  unsigned get_nvox() const;
  unsigned get_tar(const unsigned, const unsigned) const;
  double get_vox_radius() const;
  const Vector& get_center() const;
  Species& get_surface();
  Species& get_volume();
  Model& get_model();
  const std::string& get_name() const;
  std::vector<unsigned>& get_lattice();
private:
  void set_surface();
  void populate_mol(const unsigned);
private:
  const std::string name_;
  const double hcpx_;
  const double hcpo_;
  const double hcpz_;
  const double vox_radius_;
  const static unsigned ncol_;
  const static unsigned nlay_;
  const static unsigned nrow_;
  const static unsigned nrowm_;
  const static unsigned ncolrow_;
  const unsigned nvox_;
  const Vector length_;
  const Vector center_;
  std::vector<unsigned> lattice_;
  Model& model_;
  Species volume_;
  Species surface_;
  unsigned nbit_;
  unsigned sur_xor_;
};

#endif /* __Compartment_hpp */

