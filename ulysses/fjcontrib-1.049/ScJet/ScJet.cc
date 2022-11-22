// $Id$
//
// Copyright (c) 2013, Oxford University.
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib ScJet.
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
#include <limits>
#include <cmath>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/NNH.hh>

#include "ScJet.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

namespace contrib{

//----------------------------------------------------------------------
/// class to help run a semi-classical jet clustering
class ScBriefJet {
public:
  void init(const PseudoJet & jet, const ScJet* plugin) {
    et  = sqrt(energy2_value(jet, plugin->energyMode()));
    rap = jet.rap();
    phi = jet.phi();
    R = plugin->R();
    Rexp = plugin->Rexp();

    et4 = et * et * et * et;
    RR2 = 1.0 / (R * R);
  }

  double distance(const ScBriefJet * jet) const {

    double sumet = et + jet->et;

    double rho = dr2(jet) * RR2;
    double dij = 0.25 * 0.25 * sumet * sumet * sumet * sumet;
    for (int n = 0; n < Rexp; ++n) dij *= rho;

    return dij;
  }

  double dr2(const ScBriefJet* jet) const {
    double drap = rap - jet->rap;
    double dphi = abs(phi - jet->phi);
    if (dphi > M_PI) dphi = 2.0 * M_PI - dphi;
    return drap*drap + dphi*dphi;
  }

  double beam_distance() const {
    return et4;
  }

  static double beam_distance(const PseudoJet& jet, ScJet::energyModeType mode) {
    double d = energy2_value(jet, mode);
    return d;
  }

  static double energy2_value(const PseudoJet& jet, ScJet::energyModeType mode) {
    double x;
    switch (mode) {
      case ScJet::use_et : x = jet.Et2();
                           break;
      case ScJet::use_pt : x = jet.pt2();
                           break;
      default :            x = jet.mt2();
    }
    return x;
  }

private:
  double et, rap, phi, R;
  int Rexp;
  double et4, RR2;
};

//----------------------------------------------------------------------
/// implementation of ScJet algorithm
string ScJet::description () const {
  ostringstream desc;
  desc << "ScJet plugin using " << energyModeString()
       << " with R = " << R() << " and exponent " << Rexp() ;
  return desc.str();
}

void ScJet::run_clustering(ClusterSequence & cs) const {
  int njets = cs.jets().size();
  NNH<ScBriefJet,const ScJet> nnh(cs.jets(), this);

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j); // i,j are return values...

    if (j >= 0) {
      cs.plugin_record_ij_recombination(i, j, dij, k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
    } else {
      double diB = ScBriefJet::beam_distance(cs.jets()[i], energyMode()); // get new diB
      cs.plugin_record_iB_recombination(i, diB*diB);
      nnh.remove_jet(i);
    }
    njets--;
  }
    
}

//----------------------------------------------------------------------

} // namespace contrib

FASTJET_END_NAMESPACE
