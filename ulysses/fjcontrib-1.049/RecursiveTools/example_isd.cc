//----------------------------------------------------------------------
/// \file example_isd.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::IteratedSoftDrop class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_isd < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_isd.cc 1115 2018-04-21 13:37:04Z jthaler $
//
// Copyright (c) 2017, Jesse Thaler, Kevin Zhou
// based on arXiv:1704.06266 by Christopher Frye, Andrew J. Larkoski,
// Jesse Thaler, Kevin Zhou
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

#include <iostream>
#include <sstream>

#include <sstream>
#include <iomanip>
#include <cmath>
#include "fastjet/ClusterSequence.hh"

#include "IteratedSoftDrop.hh" // In external code, this should be fastjet/contrib/IteratedSoftDrop.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  // first get some anti-kt jets
  double R = 0.5;
  JetDefinition jet_def(antikt_algorithm, R);
  double ptmin = 200.0;
  Selector pt_min_selector = SelectorPtMin(ptmin);

  ClusterSequence cs(event,jet_def);
  vector<PseudoJet> jets = pt_min_selector(sorted_by_pt(cs.inclusive_jets()));
  
  // Determine optimal scale from 1704.06266 for beta = -1
  double expected_jet_pt = 1000; // in GeV
  double NP_scale = 1; // approximately Lambda_QCD in GeV
  double optimal_z_cut = (NP_scale / expected_jet_pt / R);
  
  // set up iterated soft drop objects
  double z_cut = optimal_z_cut;
  double beta  = -1.0;
  double theta_cut = 0.0;

  contrib::IteratedSoftDrop isd(beta, z_cut, theta_cut, R);
  contrib::IteratedSoftDrop isd_ee(beta, z_cut, contrib::RecursiveSoftDrop::theta_E,
                                   theta_cut, R, std::numeric_limits<double>::infinity(),
                                   contrib::RecursiveSoftDrop::larger_E);
  
  cout << "---------------------------------------------------" << endl;
  cout << "Iterated Soft Drop" << endl;
  cout << "---------------------------------------------------" << endl;
  
  cout << endl;
  cout << "Computing with:" << endl;
  cout << "     " << isd.description() << endl;
  cout << "     " << isd_ee.description() << endl;
  cout << endl;

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    
    cout << "---------------------------------------------------" << endl;
    cout << "Processing Jet " << ijet << endl;
    cout << "---------------------------------------------------" << endl;

    
    cout << endl;
    cout << "Jet pT: " << jets[ijet].pt() << " GeV" << endl;
    cout << endl;
    
    
    // Run full iterated soft drop
    //
    // Instead of asking for an IteratedSoftDropInfo which contains
    // all the information to calculate physics observables, one can
    // use
    //   isd.all_zg_thetag(jet)    for the list of symmetry factors and angles
    //   isd.multiplicity(jet)     for the ISD multiplicity
    //   isd.angularity(jet,alpha) for the angularity obtained from the (zg, thetag)
    //
    contrib::IteratedSoftDropInfo syms = isd(jets[ijet]);
    
    cout << "Soft Drop Multiplicity (pt_R measure, beta=" << beta << ", z_cut=" << z_cut << "):" << endl;
    cout << syms.multiplicity() << endl;
    cout << endl;
    
    cout << "Symmetry Factors (pt_R measure, beta=" << beta << ", z_cut=" << z_cut << "):" << endl;
    for (unsigned i = 0; i < syms.size(); i++){
      cout << syms[i].first << " ";
    }
    cout << endl;
    cout << endl;
    
    cout << "Soft Drop Angularities (pt_R measure, beta=" << beta << ", z_cut=" << z_cut << "):" << endl;
    cout << "  alpha = 0,   kappa = 0 :  " << syms.angularity(0.0, 0.0) << endl;
    cout << "  alpha = 0,   kappa = 2 :  " << syms.angularity(0.0, 2.0) << endl;
    cout << "  alpha = 0.5, kappa = 1 :  " << syms.angularity(0.5) << endl;
    cout << "  alpha = 1,   kappa = 1 :  " << syms.angularity(1.0) << endl;
    cout << "  alpha = 2,   kappa = 1 :  " << syms.angularity(2.0) << endl;
    cout << endl;
    
    // Alternative version with e+e- measure
    contrib::IteratedSoftDropInfo syms_ee = isd_ee(jets[ijet]);
    cout << "Symmetry Factors (E_theta measure, beta=" << beta << ", z_cut=" << z_cut << "):" << endl;
    for (unsigned i = 0; i < syms_ee.size(); i++){
      cout << syms_ee[i].first << " ";
    }
    cout << endl;
    cout << endl;
    
    
  }

  return 0;
}

//----------------------------------------------------------------------
/// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}
