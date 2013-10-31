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


#include <iostream> 
#include <boost/date_time/posix_time/posix_time.hpp>
#include <Compartment.hpp>
#include <Stepper.hpp>
#include <Model.hpp>
#include <VisualLogger.hpp>

int main()
{
  const double aVoxelRadius(2.5e-9);
  const double aLength(1e-6);
  Compartment aRootComp(aVoxelRadius, aLength, aLength, aLength);
  Species A(10000, 1e-12, aRootComp);
  A.populate();
  Stepper aStepper(aRootComp, A);
  Model aModel(aStepper);
  VisualLogger visualLogger(aRootComp, aStepper);
  aStepper.setLogger(visualLogger);
  visualLogger.addSpecies(A);
  visualLogger.initialize();

  boost::posix_time::ptime start(
                 boost::posix_time::microsec_clock::universal_time()); 
  //aModel.run(0.001);
  aModel.run(0.1);
  boost::posix_time::ptime end(
                 boost::posix_time::microsec_clock::universal_time());
  std::cout << "duration:" << end-start << std::endl;
}
