//  ClusteringVetoPlugin Package
//  Questions/Comments? liew@hep-th.phys.s.u-tokyo.ac.jp
//                      stoll@hep-th.phys.s.u-tokyo.ac.jp
//
//  Copyright (c) 2014-2015
//  Seng Pei Liew, Martin Stoll
//
//  $Id: ClusteringVetoPlugin.cc 792 2015-05-04 03:42:26Z martinstoll $
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

#include "ClusteringVetoPlugin.hh"
#include <fastjet/NNH.hh>

#include <cstdio>
#include "math.h"
#include <iomanip>
#include <cmath>
#include <sstream>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
      
  // Default constructor.
  ClusteringVetoPlugin::ClusteringVetoPlugin(double mu, double theta, double max_r, ClusterType clust_type)
    :_max_r2(max_r*max_r),  _mu(mu), _theta(theta),  _clust_type(clust_type),
     _veto_function(NULL)
  {
    // Added errors for user input.
    if (mu < 0.0) throw Error("ClusteringVetoPlugin: mu must be positive.");
    if (theta > 1.0 || theta < 0.0) throw Error("ClusteringVetoPlugin: theta must be in [0.0,1.0].");
    if (max_r < 0.0) throw Error("ClusteringVetoPlugin: Maximum radius must be positive.");
  }

  // Implements MJ clustering algorithm (virtual function from JetDefinition::Plugin)
  void ClusteringVetoPlugin::run_clustering(ClusterSequence & cs) const {
    vector<int> pass_jets;

    // set up NNH
    ClusteringVetoJetInfo myinfo;
    myinfo.clust_type = _clust_type;
    myinfo.max_r2 = _max_r2;
    NNH<ClusteringVetoJet,ClusteringVetoJetInfo> nnh(cs.jets(), &myinfo); 
         
    int njets = cs.jets().size();
    while (njets > 0) {

      int i(-1), j(-1);
      double dij = nnh.dij_min(i, j);

      // If closest distance is to beam, then merge with beam
      if(j < 0) {
	// add to passive jet queue and merge with beam
	pass_jets.push_back(i);
	cs.plugin_record_iB_recombination(i,dij);
	nnh.remove_jet(i);
	njets--;
	continue;
      }

      // Else check veto
      switch ( CheckVeto ( cs.jets()[i],cs.jets()[j] ) ) {
      case CLUSTER: { // below mass threshold mu
	int k=-1;
	cs.plugin_record_ij_recombination(i, j, dij, k);
	nnh.merge_jets(i, j, cs.jets()[k], k);
	njets--;
	break;
      }
      case VETO: // MJ veto called
	// add to passive jet queue and merge with beam
	pass_jets.push_back(i);
	pass_jets.push_back(j);
	cs.plugin_record_iB_recombination(i,dij);
	cs.plugin_record_iB_recombination(j,dij);
	nnh.remove_jet(i);
	nnh.remove_jet(j);
	njets=njets-2;
	break;
      case NOVETO: // check active-passive veto
	int pass_1(-1), pass_2(-1);
	double d_1(dij), d_2(dij);

	// find closest passive veto candidates
	for ( unsigned jj=0; jj < pass_jets.size(); ++jj ) {
	  double dd_1 =
	    GetJJDistanceMeasure(cs.jets()[pass_jets[jj]],cs.jets()[i]),
	    dd_2 =
	    GetJJDistanceMeasure(cs.jets()[pass_jets[jj]],cs.jets()[j]),
	    db = GetJBDistanceMeasure(cs.jets()[pass_jets[jj]]);

	  // check if they could recombine in terms of the distances
	  if ( dd_1 < d_1 && dd_1 < db )
	    { d_1 = dd_1; pass_1 = (int)pass_jets[jj]; }
	  if ( dd_2 < d_2 && dd_2 < db )
	    { d_2 = dd_2; pass_2 = (int)pass_jets[jj]; }
	}

	// check veto
	bool had_active_passive_veto = false;
	if ( pass_1 >= 0 &&
	     CheckVeto(cs.jets()[i], cs.jets()[pass_1])==VETO) {

	  // add to passive jet queue and merge with beam
	  pass_jets.push_back(i);
	  cs.plugin_record_iB_recombination(i,GetJJDistanceMeasure(cs.jets()[i],cs.jets()[pass_1]));
	  nnh.remove_jet(i);
	  had_active_passive_veto = true;
	  njets--;
	}
	if ( pass_2 >= 0 &&
	     CheckVeto(cs.jets()[j], cs.jets()[pass_2])==VETO) {

	  // add to passive jet queue and merge with beam
	  pass_jets.push_back(j);
	  cs.plugin_record_iB_recombination(j,GetJJDistanceMeasure(cs.jets()[j],cs.jets()[pass_2]));
	  nnh.remove_jet(j);
	  had_active_passive_veto = true;
	  njets--;
	}

	// no veto has been called: merge jets
	if ( !had_active_passive_veto ) {
	  int k=-1;
	  cs.plugin_record_ij_recombination(i, j, dij, k);
	  nnh.merge_jets(i, j, cs.jets()[k], k);
	  njets--;
	}
      }
    }

  }
      
  // Description of algorithm, including parameters
  string ClusteringVetoPlugin::description() const{
    stringstream sstream("");

    sstream << "Clustering Veto (1410.4637), ";
         
    switch (_clust_type) {
    case AKTLIKE:
      sstream << "AKT";
      break;
    case CALIKE:
      sstream << "CA";
      break;
    case KTLIKE:
      sstream << "KT";
      break;
    }
    sstream << "-like";
         
    sstream << fixed << setprecision(1) << ", theta=" << _theta;
    sstream << ", mu=" << _mu;
    sstream << ", max_r=" << sqrt(_max_r2);
    if ( _veto_function != NULL )
      sstream << ", have user-defined veto function";
         
    return sstream.str();
  }

  // get the dij between two jets
  double ClusteringVetoPlugin::GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2) const {
    double ret;
    switch(_clust_type) {
    case AKTLIKE:
      ret = min(1./j1.perp2(), 1./j2.perp2());
      break;
    case CALIKE :
      ret = 1.;
      break;
    case KTLIKE:
      ret = min(j1.perp2(), j2.perp2());
      break;
    default:
      assert(false);
    }
         
    ret *= j1.squared_distance(j2) / _max_r2;
    return ret;
  }
      
  // jet diB between jet and beam
  double ClusteringVetoPlugin::GetJBDistanceMeasure(const PseudoJet& jet) const{
    switch(_clust_type) {
    case AKTLIKE:
      return 1./jet.perp2();
    case CALIKE :
      return 1.;
    case KTLIKE:
      return jet.perp2();
    default:
      assert(false);
    }
  }

  // check the terminating veto
  ClusteringVetoPlugin::VetoResult ClusteringVetoPlugin::CheckVeto(const PseudoJet& j1, const PseudoJet& j2) const {
    if ( _veto_function == NULL )
      return CheckVeto_MJ(j1,j2);
    return _veto_function(j1,j2);
  }

  // check the MJ terminating veto
  ClusteringVetoPlugin::VetoResult ClusteringVetoPlugin::CheckVeto_MJ(const PseudoJet& j1, const PseudoJet& j2) const {

    PseudoJet combj = j1+j2;

    double mj1 = abs (j1.m());
    double mj2 = abs (j2.m());      
    double mcombj = abs (combj.m());

    if (mcombj < _mu)  return CLUSTER; // recombine
    else if (_theta*mcombj > max(mj1,mj2)) return VETO; // label passive
    else return NOVETO; // mass jump step 3 (check active-passive veto)
  }
      
} // namespace contrib

FASTJET_END_NAMESPACE
