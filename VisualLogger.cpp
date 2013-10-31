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

#include <VisualLogger.hpp>

void VisualLogger::initializeLog()
{
  const unsigned latticeType(0); //HCP
  _logFile.write((char*)(&latticeType), sizeof(latticeType));
  const unsigned meanCount(0);
  _logFile.write((char*)(&meanCount), sizeof(meanCount));
  const unsigned startCoord(0);
  _logFile.write((char*)(&startCoord), sizeof(startCoord));
  const unsigned colSize(_comp.getCols());
  _logFile.write((char*)(&colSize), sizeof(colSize));
  const unsigned layerSize(_comp.getLays());
  _logFile.write((char*)(&layerSize), sizeof(layerSize));
  const unsigned rowSize(_comp.getRows());
  _logFile.write((char*)(&rowSize), sizeof(rowSize));
  const double voxRadius(_comp.getVoxRadius());
  Vector center(_comp.getCenter());
  const double realColSize(center.x*2/(voxRadius*2));
  _logFile.write((char*)(&realColSize), sizeof(realColSize));
  const double realLayerSize(center.y*2/(voxRadius*2));
  _logFile.write((char*)(&realLayerSize), sizeof(realLayerSize));
  const double realRowSize(center.z*2/(voxRadius*2));
  _logFile.write((char*)(&realRowSize), sizeof(realRowSize));
  const unsigned latticeSpSize(_species.size());
  _logFile.write((char*)(&latticeSpSize), sizeof(latticeSpSize));
  const unsigned polymerSize(0);
  _logFile.write((char*)(&polymerSize), sizeof(polymerSize));
  const unsigned reservedSize(0);
  _logFile.write((char*)(&reservedSize), sizeof(reservedSize));
  const unsigned offLatticeSpSize(0);
  _logFile.write((char*)(&offLatticeSpSize), sizeof(offLatticeSpSize));
  _logFile.write((char*)(&_marker), sizeof(_marker));
  _logFile.write((char*)(&voxRadius), sizeof(voxRadius));
  for(unsigned i(0); i != _species.size(); ++i)
    {
      std::string string("Hello");
      const unsigned stringSize(string.size());
      _logFile.write((char*)(&stringSize), sizeof(stringSize));
      _logFile.write(string.c_str(), stringSize);
      _logFile.write((char*)(&voxRadius), sizeof(voxRadius));
    }
}

void VisualLogger::logMolecules(const unsigned index)
{
  Species& species(*_species[index]);
  //No need to log lipid or non diffusing vacant molecules since we have
  //already logged them once during initialization:
  if(species.getIsCompVacant())
    {
      return;
    }
  _logFile.write((char*)(&index), sizeof(index));
  const std::vector<unsigned>& mols(species.getMols());
  const unsigned size(mols.size());
  _logFile.write((char*)(&size), sizeof(size)); 
  for(unsigned i(0); i != mols.size(); ++i)
    {
      _logFile.write((char*)(&mols[i]), sizeof(mols[i]));
    }
}  

void VisualLogger::logSpecies()
{
  const double currentTime(_stepper.getCurrentTime());
  _logFile.write((char*)(&currentTime), sizeof(currentTime));
  for(unsigned i(0); i != _species.size(); ++i)
    {
      logMolecules(i);
    }
  _logFile.write((char*)(&_marker), sizeof(_marker));
  _logFile.write((char*)(&_marker), sizeof(_marker));
}

void VisualLogger::logCompVacant()
{
  const double currentTime(_stepper.getCurrentTime());
  _logFile.write((char*)(&currentTime), sizeof(currentTime));
  for(unsigned i(0); i != _species.size(); ++i)
    {
      if(_species[i]->getIsCompVacant())
        {
          Species& species(*_species[i]);
          //The species index in the process:
          _logFile.write((char*)(&i), sizeof(i)); 
          const std::vector<unsigned>& mols(species.getMols());
          const unsigned size(mols.size());
          _logFile.write((char*)(&size), sizeof(size)); 
          for(unsigned i(0); i != mols.size(); ++i)
            {
              _logFile.write((char*)(&mols[i]), sizeof(mols[i]));
            }
        }
    }
  _logFile.write((char*)(&_marker), sizeof(_marker));
  _logFile.write((char*)(&_marker), sizeof(_marker));
}
