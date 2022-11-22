//----------------------------------------------------------------------
/// \file example_bottomup_softdrop.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::BottomUpSoftDrop class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_bottomup_softdrop < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_bottomup_softdrop.cc 1075 2017-09-18 15:27:09Z gsoyez $
//
// Copyright (c) 2017-, Gavin P. Salam, Gregory Soyez, Jesse Thaler,
// Kevin Zhou, Frederic Dreyer
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

#include <iomanip>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "BottomUpSoftDrop.hh" // In external code, this should be fastjet/contrib/BottomUpSoftDrop.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void print_prongs(const PseudoJet &jet, const string &pprefix);
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  // first get some anti-kt jets
  double R = 1.0, ptmin = 100.0;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

  //----------------------------------------------------------------------
  // give the soft drop groomer a short name
  // Use a symmetry cut z > z_cut R^beta
  // By default, there is no mass-drop requirement
  double z_cut = 0.2;
  double beta  = 1.0;
  contrib::BottomUpSoftDrop busd(beta, z_cut);
  
  //----------------------------------------------------------------------
  cout << "BottomUpSoftDrop groomer is: " << busd.description() << endl;

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    // Run SoftDrop and examine the output
    PseudoJet busd_jet = busd(jets[ijet]);
    cout << endl;
    cout << "original            jet: " << jets[ijet] << endl;
    cout << "BottomUpSoftDropped jet: " << busd_jet << endl;
    
    assert(busd_jet != 0); //because bottom-up soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

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

//----------------------------------------------------------------------
/// overloaded jet info output
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi();
  }
  return ostr;
}
