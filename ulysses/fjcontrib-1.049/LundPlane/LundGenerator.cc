// $Id: LundGenerator.cc 1288 2021-11-09 11:53:11Z scyboz $
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

#include <sstream>
#include "LundGenerator.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

LundDeclustering::LundDeclustering(const PseudoJet& pair,
				   const PseudoJet& j1, const PseudoJet& j2)
  : m_(pair.m()), Delta_(j1.delta_R(j2)), pair_(pair) {

  // establish which of j1 and j2 is softer
  if (j1.pt2() > j2.pt2()) {
    harder_ = j1;
    softer_ = j2;
  } else {
    harder_ = j2;
    softer_ = j1;
  }

  // now work out the various Lund declustering variables
  double softer_pt = softer_.pt();
  z_   = softer_pt / (softer_pt + harder_.pt());
  kt_  = softer_pt * Delta_;
  psi_ = atan2(softer_.rap()-harder_.rap(), harder_.delta_phi_to(softer_));
  kappa_ = z_ * Delta_;
}

//----------------------------------------------------------------------
/// retrieve the vector of declusterings of the primary plane of a jet
std::vector<LundDeclustering> LundGenerator::result(const PseudoJet& jet) const {
  std::vector<LundDeclustering> result;
  PseudoJet j = recluster_(jet);
  PseudoJet pair, j1, j2;
  pair = j;
  while (pair.has_parents(j1, j2)) {
    if (j1.pt2() < j2.pt2()) std::swap(j1,j2);
    result.push_back(LundDeclustering(pair, j1, j2));
    pair = j1;
  }
  return result;
}

//----------------------------------------------------------------------
/// description
std::string LundGenerator::description() const {
  std::ostringstream oss;
  oss << "LundGenerator with " << recluster_.description();
  return oss.str();
}

} // namespace contrib

FASTJET_END_NAMESPACE
