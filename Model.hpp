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


#ifndef __Model_hpp
#define __Model_hpp

#include <Common.hpp>
#include <Compartment.hpp>
#include <Stepper.hpp>

class Model
{ 
public: 
  Model(const double, const double, const double, const double);
  ~Model() {}
  void run(const double);
  unsigned push_species(Species&);
  Compartment& get_comp();
  Stepper& get_stepper();
private:
  std::vector<Species*> species_;
  Stepper stepper_;
  Compartment comp_; //must declare this at the end after initializing others
};

#endif /* __Model_hpp */

