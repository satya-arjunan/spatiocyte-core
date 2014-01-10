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

#include <sstream>
#include <climits>
#include <Spatiocyte.hpp>
#include <VisualLogger.hpp>
#include <Model.hpp>
#include <Stepper.hpp>
#include <Compartment.hpp>


VisualLogger::VisualLogger(Model& model):
  marker_(UINT_MAX),
  filename_("VisualLog.dat"),
  comp_(model.get_comp()),
  stepper_(model.get_stepper()) {}

void VisualLogger::fire()
{
  log_species();
  logfile_.flush();
}

void VisualLogger::initialize()
{
  std::ostringstream fileName;
  fileName << filename_ << std::ends;
  logfile_.open(fileName.str().c_str(), std::ios::binary | std::ios::trunc);
  initialize_log();
  log_structure_species();
  log_species();
  logfile_.flush();
}

void VisualLogger::push_species(Species& species)
{
  species_.push_back(&species);
}

void VisualLogger::initialize_log()
{
  const unsigned latticeType(0); //HCP
  logfile_.write((char*)(&latticeType), sizeof(latticeType));
  const unsigned meanCount(0);
  logfile_.write((char*)(&meanCount), sizeof(meanCount));
  const unsigned startCoord(0);
  logfile_.write((char*)(&startCoord), sizeof(startCoord));
  const unsigned colSize(comp_.get_num_col());
  logfile_.write((char*)(&colSize), sizeof(colSize));
  const unsigned layerSize(comp_.get_num_lay());
  logfile_.write((char*)(&layerSize), sizeof(layerSize));
  const unsigned rowSize(comp_.get_num_row());
  logfile_.write((char*)(&rowSize), sizeof(rowSize));
  const double voxRadius(VOXEL_RADIUS);
  Vector center(comp_.get_center());
  const double realColSize(center.x*2/(voxRadius*2));
  logfile_.write((char*)(&realColSize), sizeof(realColSize));
  const double realLayerSize(center.y*2/(voxRadius*2));
  logfile_.write((char*)(&realLayerSize), sizeof(realLayerSize));
  const double realRowSize(center.z*2/(voxRadius*2));
  logfile_.write((char*)(&realRowSize), sizeof(realRowSize));
  const unsigned latticeSpSize(species_.size());
  logfile_.write((char*)(&latticeSpSize), sizeof(latticeSpSize));
  const unsigned polymerSize(0);
  logfile_.write((char*)(&polymerSize), sizeof(polymerSize));
  const unsigned reservedSize(0);
  logfile_.write((char*)(&reservedSize), sizeof(reservedSize));
  const unsigned offLatticeSpSize(0);
  logfile_.write((char*)(&offLatticeSpSize), sizeof(offLatticeSpSize));
  logfile_.write((char*)(&marker_), sizeof(marker_));
  logfile_.write((char*)(&voxRadius), sizeof(voxRadius));
  for(unsigned i(0); i != species_.size(); ++i)
    {
      const unsigned stringSize(species_[i]->get_name_id().size());
      logfile_.write((char*)(&stringSize), sizeof(stringSize));
      logfile_.write(species_[i]->get_name_id().c_str(), stringSize);
      logfile_.write((char*)(&voxRadius), sizeof(voxRadius));
    }
}

void VisualLogger::log_structure_species()
{
  const double currentTime(stepper_.get_current_time());
  logfile_.write((char*)(&currentTime), sizeof(currentTime));
  for(unsigned i(0); i != species_.size(); ++i)
    {
      if(species_[i]->is_structure_species())
        {
          Species& species(*species_[i]);
          //The species index in the process:
          logfile_.write((char*)(&i), sizeof(i)); 
          const std::vector<umol_t>& mols(species.get_mols());
          const unsigned size(mols.size());
          logfile_.write((char*)(&size), sizeof(size)); 
          for(unsigned j(0); j != mols.size(); ++j)
            {
              unsigned mol(mols[j]);
              logfile_.write((char*)(&mol), sizeof(mol));
            }
        }
    }
  logfile_.write((char*)(&marker_), sizeof(marker_));
  logfile_.write((char*)(&marker_), sizeof(marker_));
}

void VisualLogger::log_species()
{
  const double currentTime(stepper_.get_current_time());
  logfile_.write((char*)(&currentTime), sizeof(currentTime));
  for(unsigned i(0); i != species_.size(); ++i)
    {
      log_mols(i);
    }
  logfile_.write((char*)(&marker_), sizeof(marker_));
  logfile_.write((char*)(&marker_), sizeof(marker_));
}

void VisualLogger::log_mols(const unsigned index)
{
  Species& species(*species_[index]);
  //No need to log lipid or non diffusing vacant molecules since we have
  //already logged them once during initialization:
  if(species.is_structure_species())
    {
      return;
    }
  logfile_.write((char*)(&index), sizeof(index));
  const std::vector<umol_t>& mols(species.get_mols());
  const unsigned size(mols.size());
  logfile_.write((char*)(&size), sizeof(size)); 
  for(unsigned i(0); i != mols.size(); ++i)
    {
      unsigned mol(mols[i]);
      logfile_.write((char*)(&mol), sizeof(mol));
    }
}  

