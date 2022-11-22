// $Id: LundWithSecondary.cc 1289 2021-11-09 11:53:53Z scyboz $
//
// Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
// Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "LundWithSecondary.hh"
#include <sstream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
/// return LundDeclustering sequence of primary plane
std::vector<LundDeclustering> LundWithSecondary::primary(const PseudoJet& jet) const {
  return lund_gen_(jet);
}

//----------------------------------------------------------------------
/// return LundDeclustering sequence of secondary plane (slow version)
std::vector<LundDeclustering> LundWithSecondary::secondary(const PseudoJet& jet) const {
  // this is not optimal as one is computing the primary plane twice.
  std::vector<LundDeclustering> declusts = lund_gen_(jet);
  return secondary(declusts);
}
//----------------------------------------------------------------------
/// return LundDeclustering sequence of secondary plane with primary sequence as input
std::vector<LundDeclustering> LundWithSecondary::secondary(
				  const std::vector<LundDeclustering> & declusts) const {

  int sec_index = secondary_index(declusts);
  // if we found the index of secondary emission, return its declustering sequence
  if (sec_index >= 0) {
    return lund_gen_(declusts[sec_index].softer());
  } else  {
    return std::vector<LundDeclustering>();
  }
}

//----------------------------------------------------------------------
/// return the index associated of the primary declustering that is to be
/// used for the secondary plane.
int LundWithSecondary::secondary_index(const std::vector<LundDeclustering> & declusts) const {
  if (secondary_def_ == 0) {
    throw Error("secondary class is a null pointer, cannot identify element to use for secondary plane");
  }
  return (*secondary_def_)(declusts);
}
  
//----------------------------------------------------------------------
/// description
std::string LundWithSecondary::description() const {
  std::ostringstream oss;
  oss << "LundWithSecondary using " << secondary_def_->description()
      << " and " << lund_gen_.description();
  return oss.str();
}

} // namespace contrib

FASTJET_END_NAMESPACE
