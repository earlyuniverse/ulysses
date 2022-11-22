// $Id: SecondaryLund.cc 1305 2021-12-06 11:36:03Z gsoyez $
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

#include "SecondaryLund.hh"
#include <math.h>
#include <sstream>
#include <limits>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
/// retrieve the vector of declusterings of the secondary plane of a jet
int SecondaryLund_mMDT::result(const std::vector<LundDeclustering>& declusts) const {
  // iterate through primary branchings
  for (unsigned int i=0; i < declusts.size(); ++i) {
    // mMDTZ: find the first emission passing z>zcut
    if (declusts[i].z() > zcut_) return i;
  }
  return -1;
}

int SecondaryLund_dotmMDT::result(const std::vector<LundDeclustering>& declusts) const {
  // set up leading emission pointer and bookkeeping variables
  int secondary_index  = -1;
  double dot_prod_max   = 0.0;
  
  // iterate through primary branchings
  for (unsigned int i=0; i < declusts.size(); ++i) {
    // dotmMDTZ: find emission passing z>zcut with largest dot product
    if (declusts[i].z() > zcut_) {
      double dot_prod = declusts[i].harder().pt()*declusts[i].softer().pt()
	*declusts[i].Delta()*declusts[i].Delta();
      if (dot_prod > dot_prod_max) {
	dot_prod_max    = dot_prod;
	secondary_index = i;
      }
    }
  }
  
  // return index of secondary emission
  return secondary_index;
}


int SecondaryLund_Mass::result(const std::vector<LundDeclustering>& declusts) const {
  // set up leading emission pointer and bookkeeping variables
  int secondary_index  = -1;
  double mass_diff     = std::numeric_limits<double>::max();
  
  // iterate through primary branchings
  for (unsigned int i=0; i < declusts.size(); ++i) {
    // Mass: find emission that minimises the distance to reference mass
    double dist =
      std::abs(log(declusts[i].harder().pt()*declusts[i].softer().pt()
		   * declusts[i].Delta()*declusts[i].Delta() / mref2_)
	       * log(1.0/declusts[i].z()));
    if (dist < mass_diff) {
      mass_diff     = dist;
      secondary_index = i;
    }
  }
  
  // return index of secondary emission
  return secondary_index;
}


//----------------------------------------------------------------------
/// description
std::string SecondaryLund::description() const {
  std::ostringstream oss;
  oss << "SecondaryLund";
  return oss.str();
}

//----------------------------------------------------------------------
/// description
std::string SecondaryLund_mMDT::description() const {
  std::ostringstream oss;
  oss << "SecondaryLund (mMDT selection of leading emission, zcut=" << zcut_<<")";
  return oss.str();
}

//----------------------------------------------------------------------
/// description
std::string SecondaryLund_dotmMDT::description() const {
  std::ostringstream oss;
  oss << "SecondaryLund (dotmMDT selection of leading emission, zcut=" << zcut_<<")";
  return oss.str();
}

//----------------------------------------------------------------------
/// description
std::string SecondaryLund_Mass::description() const {
  std::ostringstream oss;
  oss << " (Mass selection of leading emission, m="<<sqrt(mref2_)<<")";
  return oss.str();
}

} // namespace contrib

FASTJET_END_NAMESPACE
