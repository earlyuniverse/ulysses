// JetCleanser Package
// Questions/Comments? dkrohn@physics.harvard.edu mattlow@uchicago.edu schwartz@physics.harvard.edu liantaow@uchicago.edu
//
// Copyright (c) 2013
// David Krohn, Matthew Low, Matthew Schwartz, and Lian-Tao Wang
//
// Compile it with "make example" and run it with
//
//   ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//
// $Id$
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
#include <cstdlib>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "JetCleanser.hh" // In external code, this should be fastjet/contrib/JetCleanser.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &hard_event_charged, vector<PseudoJet> &hard_event_neutral, 
                vector<PseudoJet> &pileup_charged,     vector<PseudoJet> &pileup_neutral);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> hard_event_charged, hard_event_neutral, pileup_charged, pileup_neutral;
  read_event(hard_event_charged, hard_event_neutral, pileup_charged, pileup_neutral);

  //----------------------------------------------------------
  // construct measureable quantities
  vector<PseudoJet> hard_event, full_event;  // via "calo cells"
  vector<PseudoJet> full_event_neutral;      // via "particle flow"

  for (unsigned i=0; i<hard_event_charged.size(); i++) {
    hard_event.push_back( hard_event_charged[i] ); 
    full_event.push_back( hard_event_charged[i] ); }

  for (unsigned i=0; i<hard_event_neutral.size(); i++) {
    hard_event.push_back( hard_event_neutral[i] ); 
    full_event.push_back( hard_event_neutral[i] );
    full_event_neutral.push_back( hard_event_neutral[i] ); }

  for (unsigned i=0; i<pileup_charged.size(); i++) { 
    full_event.push_back( pileup_charged[i] ); }

  for (unsigned i=0; i<pileup_neutral.size(); i++) { 
    full_event.push_back( pileup_neutral[i] ); 
    full_event_neutral.push_back( pileup_neutral[i] ); }

  //----------------------------------------------------------
  // illustrate how this JetCleanser contrib works

  // find jets
  JetDefinition jet_def( antikt_algorithm, 1.0 );
  vector< vector<PseudoJet> > sets;
  sets.push_back( full_event );           // calorimeter cells
  sets.push_back( hard_event_charged );   // tracks from primary interaction
  sets.push_back( pileup_charged );       // tracks from pileup
  sets.push_back( full_event_neutral );   // neutral particles
  sets.push_back( hard_event );           // for truth comparison

  // collect jets
  vector< vector<PseudoJet> > jet_sets = ClusterSets(jet_def, full_event, sets, 25.0);
  vector<PseudoJet> jets_plain     = jet_sets[0];
  vector<PseudoJet> jets_tracks_LV = jet_sets[1];
  vector<PseudoJet> jets_tracks_PU = jet_sets[2];
  vector<PseudoJet> jets_neutrals  = jet_sets[3];
  vector<PseudoJet> jets_truth     = jet_sets[4];

  PseudoJet plain_dijet, truth_dijet, jvf_cleansed_dijet, lin_cleansed_dijet, gau_cleansed_dijet;
  unsigned n_jets;

  //----------------------------------------------------------
  // ATLAS-like: cleansers
  cout << "ATLAS-like cleansing:" << endl << endl;

  JetDefinition subjet_def_A(kt_algorithm, 0.3);
  JetCleanser jvf_cleanser_A(subjet_def_A, JetCleanser::jvf_cleansing, JetCleanser::input_nc_together);
  jvf_cleanser_A.SetTrimming(0.01);

  JetCleanser linear_cleanser_A(0.25, JetCleanser::linear_cleansing, JetCleanser::input_nc_together);
  linear_cleanser_A.SetLinearParameters(0.65);

  JetCleanser gaussian_cleanser_A(0.3, JetCleanser::gaussian_cleansing, JetCleanser::input_nc_together);
  gaussian_cleanser_A.SetGaussianParameters(0.67,0.62,0.20,0.25);

  // print info about cleansers
  cout << jvf_cleanser_A.description() << endl;
  cout << linear_cleanser_A.description() << endl;
  cout << gaussian_cleanser_A.description() << endl;
 
  // ATLAS-like: cleanse jets
  n_jets = min((int) jets_plain.size(),3);
  for (unsigned i=0; i<n_jets; i++){
    PseudoJet plain_jet = jets_plain[i];
    PseudoJet truth_jet = jets_truth[i];
    PseudoJet jvf_cleansed_jet = jvf_cleanser_A( jets_plain[i], jets_tracks_LV[i].constituents(), jets_tracks_PU[i].constituents() );
    PseudoJet lin_cleansed_jet = linear_cleanser_A( jets_plain[i], jets_tracks_LV[i].constituents(), jets_tracks_PU[i].constituents() );
    PseudoJet gau_cleansed_jet = gaussian_cleanser_A( jets_plain[i], jets_tracks_LV[i].constituents(), jets_tracks_PU[i].constituents() );
    
    cout << "                  no pileup: pt = " << truth_jet.pt()
                                    << " eta = " << truth_jet.eta()
                                    << " phi = " << truth_jet.phi()
                                    << "   m = " << truth_jet.m()
                                    << endl;

    cout << "                with pileup: pt = " << plain_jet.pt()
                                    << " eta = " << plain_jet.eta()
                                    << " phi = " << plain_jet.phi()
                                    << "   m = " << plain_jet.m()
                                    << endl;

    cout << " with pileup + jvf cleansed: pt = " << jvf_cleansed_jet.pt()
                                    << " eta = " << jvf_cleansed_jet.eta()
                                    << " phi = " << jvf_cleansed_jet.phi()
                                    << "   m = " << jvf_cleansed_jet.m()
                                    << endl;

    cout << " with pileup + lin cleansed: pt = " << lin_cleansed_jet.pt()
                                    << " eta = " << lin_cleansed_jet.eta()
                                    << " phi = " << lin_cleansed_jet.phi()
                                    << "   m = " << lin_cleansed_jet.m()
                                    << endl;

    cout << " with pileup + gau cleansed: pt = " << gau_cleansed_jet.pt()
                                    << " eta = " << gau_cleansed_jet.eta()
                                    << " phi = " << gau_cleansed_jet.phi()
                                    << "   m = " << gau_cleansed_jet.m()
                                    << endl
                                    << endl;

    if ( i<2 ) { 
      plain_dijet = plain_dijet + plain_jet;
      truth_dijet = truth_dijet + truth_jet;
      jvf_cleansed_dijet = jvf_cleansed_dijet + jvf_cleansed_jet;
      lin_cleansed_dijet = lin_cleansed_dijet + lin_cleansed_jet;
      gau_cleansed_dijet = gau_cleansed_dijet + gau_cleansed_jet;
    }
  }

  cout << "Dijet Masses: " << endl
       << " plain = " << plain_dijet.m() << endl
       << " truth = " << truth_dijet.m() << endl
       << " jvf   = " << jvf_cleansed_dijet.m() << endl
       << " lin   = " << lin_cleansed_dijet.m() << endl
       << " gau   = " << gau_cleansed_dijet.m() << endl << endl;

  plain_dijet = PseudoJet();
  truth_dijet = PseudoJet();
  jvf_cleansed_dijet = PseudoJet();
  lin_cleansed_dijet = PseudoJet();
  gau_cleansed_dijet = PseudoJet();

  //----------------------------------------------------------
  // CMS-like: cleansers
  cout << "CMS-like cleansing:" << endl << endl;

  JetDefinition subjet_def_B(kt_algorithm, 0.3);
  JetCleanser jvf_cleanser_B(subjet_def_B, JetCleanser::jvf_cleansing, JetCleanser::input_nc_separate);
  jvf_cleanser_B.SetTrimming(0.01);

  JetCleanser linear_cleanser_B(0.25, JetCleanser::linear_cleansing, JetCleanser::input_nc_separate);
  linear_cleanser_B.SetLinearParameters(0.65);

  JetCleanser gaussian_cleanser_B(0.3, JetCleanser::gaussian_cleansing, JetCleanser::input_nc_separate);
  gaussian_cleanser_B.SetGaussianParameters(0.67,0.62,0.20,0.25);

  // print info about cleansers
  cout << jvf_cleanser_B.description() << endl;
  cout << linear_cleanser_B.description() << endl;
  cout << gaussian_cleanser_B.description() << endl;

  // CMS-like: cleanse jets
  n_jets = min((int) jets_plain.size(),3);
  for (unsigned i=0; i<n_jets; i++){
    PseudoJet plain_jet = jets_plain[i];
    PseudoJet truth_jet = jets_truth[i];
    PseudoJet jvf_cleansed_jet = jvf_cleanser_B( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(), 
                                                 jets_tracks_PU[i].constituents() );
    PseudoJet lin_cleansed_jet = linear_cleanser_B( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(), 
                                                    jets_tracks_PU[i].constituents() );
    PseudoJet gau_cleansed_jet = gaussian_cleanser_B( jets_neutrals[i].constituents(), jets_tracks_LV[i].constituents(), 
                                                      jets_tracks_PU[i].constituents() );
    
    cout << "                  no pileup: pt = " << truth_jet.pt()
                                    << " eta = " << truth_jet.eta()
                                    << " phi = " << truth_jet.phi()
                                    << "   m = " << truth_jet.m()
                                    << endl;

    cout << "                with pileup: pt = " << plain_jet.pt()
                                    << " eta = " << plain_jet.eta()
                                    << " phi = " << plain_jet.phi()
                                    << "   m = " << plain_jet.m()
                                    << endl;

    cout << " with pileup + jvf cleansed: pt = " << jvf_cleansed_jet.pt()
                                    << " eta = " << jvf_cleansed_jet.eta()
                                    << " phi = " << jvf_cleansed_jet.phi()
                                    << "   m = " << jvf_cleansed_jet.m()
                                    << endl;

    cout << " with pileup + lin cleansed: pt = " << lin_cleansed_jet.pt()
                                    << " eta = " << lin_cleansed_jet.eta()
                                    << " phi = " << lin_cleansed_jet.phi()
                                    << "   m = " << lin_cleansed_jet.m()
                                    << endl;

    cout << " with pileup + gau cleansed: pt = " << gau_cleansed_jet.pt()
                                    << " eta = " << gau_cleansed_jet.eta()
                                    << " phi = " << gau_cleansed_jet.phi()
                                    << "   m = " << gau_cleansed_jet.m()
                                    << endl
                                    << endl;

    if ( i<2 ) { 
      plain_dijet = plain_dijet + plain_jet;
      truth_dijet = truth_dijet + truth_jet;
      jvf_cleansed_dijet = jvf_cleansed_dijet + jvf_cleansed_jet;
      lin_cleansed_dijet = lin_cleansed_dijet + lin_cleansed_jet;
      gau_cleansed_dijet = gau_cleansed_dijet + gau_cleansed_jet;
    }
  }

  cout << "Dijet Masses: " << endl
       << " plain = " << plain_dijet.m() << endl
       << " truth = " << truth_dijet.m() << endl
       << " jvf   = " << jvf_cleansed_dijet.m() << endl
       << " lin   = " << lin_cleansed_dijet.m() << endl
       << " gau   = " << gau_cleansed_dijet.m() << endl << endl;

  return 0;
}

//------------------------------------------------------------------------
// read the event with and without pileup
void read_event(vector<PseudoJet> &hard_event_charged, 
                vector<PseudoJet> &hard_event_neutral, 
                vector<PseudoJet> &pileup_charged,
                vector<PseudoJet> &pileup_neutral){
  string line;
  int  nsub  = 0; // counter to keep track of which sub-event we're reading
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {break;}
    if (line.substr(0,9) == "#SUBSTART") {
      // if more sub events follow, make copy of first one (the hard one) here
      //if (nsub == 1) {hard_event_charged = full_event_charged; hard_event_neutral = full_event_neutral;}
      nsub += 1;
    }
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    int pid,charge;
    linestream >> px >> py >> pz >> E >> pid >> charge;
    PseudoJet particle(px,py,pz,E);

    if ( nsub <= 1 ) {
      if ( charge != 0 ) hard_event_charged.push_back(particle);
      else hard_event_neutral.push_back(particle);
    } else {
      if ( charge != 0 ) pileup_charged.push_back(particle);
      else pileup_neutral.push_back(particle);
    }
  }

  // if there was nothing in the event 
  if (nsub == 0) {
    cerr << "Error: read empty event\n";
    exit(-1);
  }

  cout << "# " << nsub-1 << " pileup events on top of the hard event" << endl;
}
