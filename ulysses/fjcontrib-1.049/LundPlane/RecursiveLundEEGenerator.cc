// $Id: RecursiveLundEEGenerator.cc 1290 2021-11-09 11:55:12Z scyboz $
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
#include "RecursiveLundEEGenerator.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE

namespace contrib{
  
LundEEDeclustering::LundEEDeclustering(const PseudoJet& pair,
		const PseudoJet& j1, const PseudoJet& j2,
		int iplane, double psi, double psibar,
		int depth, int leaf_iplane, int sign_s)
  : iplane_(iplane), psi_(psi), psibar_(psibar), m_(pair.m()), pair_(pair), depth_(depth), leaf_iplane_(leaf_iplane), sign_s_(sign_s) {

  double omc = one_minus_costheta(j1,j2);
  
  if (omc > sqrt(numeric_limits<double>::epsilon())) {
    double cos_theta = 1.0 - omc;
    double theta     = acos(cos_theta);
    sin_theta_ = sin(theta);
    eta_       = -log(tan(theta/2.0));
  } else {
    // we are at small angles, so use small-angle formulas
    double theta = sqrt(2. * omc);
    sin_theta_ = theta;
    eta_       = -log(theta/2);
  }

  // establish which of j1 and j2 is softer
  if (j1.modp2() > j2.modp2()) {
    harder_ = j1;
    softer_ = j2;
  } else {
    harder_ = j2;
    softer_ = j1;
  }

  // now work out the various LundEE declustering variables
  double softer_modp = softer_.modp();
  z_   = softer_modp / (softer_modp + harder_.modp());
  kt_  = softer_modp * sin_theta_;
  lnkt_ = log(kt_);
  kappa_ = z_ * sin_theta_;
}

} // namespace contrib


FASTJET_END_NAMESPACE

std::ostream & operator<<(std::ostream & ostr, const fastjet::contrib::LundEEDeclustering & d) {
  ostr << "kt = "   << d.kt()
       << " z = "   << d.z()
       << " eta = " << d.eta()
       << " psi = " << d.psi()
       << " psibar = " << d.psibar()
       << " m = "   << d.m()
       << " iplane = " << d.iplane()
       << " depth = " << d.depth()
       << " leaf_iplane = " << d.leaf_iplane();
  return ostr;
}

