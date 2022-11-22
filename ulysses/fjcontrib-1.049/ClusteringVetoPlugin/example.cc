//  ClusteringVetoPlugin Package
//  Questions/Comments? liew@hep-th.phys.s.u-tokyo.ac.jp
//                      stoll@hep-th.phys.s.u-tokyo.ac.jp
//
//  Copyright (c) 2014-2015
//  Seng Pei Liew, Martin Stoll
//
//  Example showing usage of ClusteringVetoPlugin
//
//  Compile with "make example" and run with
//    ./example < ../data/single-event.dat
//
//  $Id: example.cc 792 2015-05-04 03:42:26Z martinstoll $
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
#include <stdio.h>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "ClusteringVetoPlugin.hh" // In external code, this should be fastjet/contrib/ClusteringVetoPlugin.hh

using namespace std;
using namespace fastjet;
using namespace contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

void print_jets(const vector<PseudoJet> &jets,
		const ClusterSequence &clust_seq);
ClusteringVetoPlugin::VetoResult user_veto_function
(const PseudoJet& j1, const PseudoJet& j2);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how this ClusteringVetoPlugin contrib works
  // basic Cambridge-Aachen-like mass-jump
  //----------------------------------------------------------

  {
    // clustering parameters
    double mu(30.), theta(0.7), max_r(1.0), ptmin(5.);

    ClusteringVetoPlugin mj_plugin(mu, theta, max_r, ClusteringVetoPlugin::CALIKE);
    JetDefinition jet_def(&mj_plugin);
    ClusterSequence clust_seq(event, jet_def);

    // print setup
    cout << endl << "Run " << jet_def.description() << endl;

    // get inclusive jets and print
    vector<PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  
    cout << "Inclusive jets with pT > " << ptmin << " GeV" << endl;
    print_jets(jets, clust_seq);
  }

  //----------------------------------------------------------
  // kT-like and anti-kT-like
  //----------------------------------------------------------

  {
    // clustering parameters
    double mu(30.), theta(0.7), max_r(1.0), ptmin(5.);

    // kT
    ClusteringVetoPlugin mj_plugin_kt(mu, theta, max_r, ClusteringVetoPlugin::KTLIKE);
    JetDefinition jet_def_kt(&mj_plugin_kt);
    ClusterSequence clust_seq_kt(event, jet_def_kt);

    cout << endl << "Run " << jet_def_kt.description() << endl;

    vector<PseudoJet> jets_kt =
      sorted_by_pt(clust_seq_kt.inclusive_jets(ptmin));
  
    cout << "Inclusive jets with pT > " << ptmin << " GeV" << endl
	 << " number of jets: " << jets_kt.size() << endl;

    // anti-kT
    ClusteringVetoPlugin mj_plugin_akt(mu, theta, max_r, ClusteringVetoPlugin::AKTLIKE);
    JetDefinition jet_def_akt(&mj_plugin_akt);
    ClusterSequence clust_seq_akt(event, jet_def_akt);

    cout << endl << "Run " << jet_def_akt.description() << endl;

    vector<PseudoJet> jets_akt =
      sorted_by_pt(clust_seq_akt.inclusive_jets(ptmin));
  
    cout << "Inclusive jets with pT > " << ptmin << " GeV" << endl
	 << " number of jets: " << jets_akt.size() << endl;
  }

  //----------------------------------------------------------
  // Cambridge-Aachen-like with user-defined veto function
  //----------------------------------------------------------

  {
    // clustering parameters
    double mu(30.), theta(0.7), max_r(1.0), ptmin(5.);

    ClusteringVetoPlugin mj_plugin(mu, theta, max_r, ClusteringVetoPlugin::CALIKE);

    // set veto function, renders mu and theta irrelevant
    mj_plugin.set_veto_function(user_veto_function);

    JetDefinition jet_def(&mj_plugin);
    ClusterSequence clust_seq(event, jet_def);

    // print setup
    cout << endl << "Run " << jet_def.description() << endl;

    // get inclusive jets and print
    vector<PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  
    cout << "Inclusive jets with pT > " << ptmin << " GeV" << endl;
    print_jets(jets, clust_seq);
  }



  return 0;
}

// read in input particles
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

// prints a vector of jets
void print_jets(const vector<PseudoJet> &jets,
		const ClusterSequence &clust_seq){

  // columns labels
  printf("%5s %10s %10s %10s %10s %10s\n",
	 "jet #", "pt", "rap", "phi", "m", "last d_ij");
   
  // print out the jets
  for (unsigned i=0; i<jets.size(); ++i) {
    printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	   i, jets[i].pt(), jets[i].rap(), jets[i].phi(), jets[i].m(),
	   clust_seq.exclusive_subdmerge ( jets[i], 1 ));
  }
}

// example of user-defined veto function
ClusteringVetoPlugin::VetoResult user_veto_function
( const PseudoJet& j1, const PseudoJet& j2 ) {
  
  // DeltaR instead of jet-mass (mu) threshold
  if ( j1.delta_R(j2) < 0.3 )
    return ClusteringVetoPlugin::CLUSTER;
  // mass-jump veto with theta=0.5
  else if ( 0.5 * (j1+j2).m() > max( j1.m(), j2.m() ) )
    return ClusteringVetoPlugin::VETO;
  // no veto, but algorithm may need to check active-passive veto
  return ClusteringVetoPlugin::NOVETO;
}
