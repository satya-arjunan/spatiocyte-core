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


#ifndef __Species_hpp
#define __Species_hpp

#include <RandomLib/Random.hpp>
#include <Common.hpp>
#include <Diffuser.hpp>

class Species
{ 
public: 
  Species(const unsigned, const double, Model&, Compartment&,
          const bool is_comp_vacant=false);
  ~Species() {}
  void populate();
  bool is_comp_vacant() const;
  unsigned get_id() const;
  Diffuser& get_diffuser();
  Compartment& get_comp();
  std::vector<unsigned>& get_mols();
private:
  const bool is_comp_vacant_;
  Compartment& comp_;
  RandomLib::Random rng_;
  std::vector<unsigned>& lattice_;
  std::vector<unsigned> mols_;
  const unsigned id_;
  Diffuser diffuser_;
  const unsigned nbit_;
  const unsigned vac_id_;
  const unsigned vac_xor_;
};

#endif /* __Species_hpp */

