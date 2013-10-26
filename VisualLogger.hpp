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


#ifndef __VisualLogger_hpp
#define __VisualLogger_hpp

#include <fstream>
#include <Species.hpp>
#include <Compartment.hpp>
#include <Stepper.hpp>
#include <climits>

class Stepper;

class VisualLogger
{
public:
  VisualLogger(Compartment& comp, Stepper& stepper):
    _marker(UINT_MAX),
    _fileName("VisualLog.dat"),
    _comp(comp),
    _stepper(stepper) {}
  virtual ~VisualLogger() {}
  void initialize()
    {
      std::ostringstream fileName;
      fileName << _fileName << std::ends;
      _logFile.open(fileName.str().c_str(), std::ios::binary | std::ios::trunc);
      initializeLog();
      logCompVacant();
      logSpecies();
      _logFile.flush();
    }
  void fire()
    {
      logSpecies();
      _logFile.flush();
    }
  void addSpecies(Species& species)
    {
      _species.push_back(&species);
    }
protected:
  virtual void initializeLog();
  virtual void logCompVacant();
  void logSpecies();
  void logMolecules(const unsigned);
protected:
  unsigned _marker;
  std::string _fileName;
  std::ofstream _logFile;
  Compartment& _comp;
  Stepper& _stepper;
  std::vector<Species*> _species;
};

#endif /* __VisualLogger_hpp */
