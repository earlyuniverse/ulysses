//----------------------------------------------------------------------
/// \file example_softdrop.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::SoftDrop class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_softdrop < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_softdrop.cc 705 2014-07-07 14:37:03Z gsoyez $
//
// Copyright (c) 2014, Gavin P. Salam
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
#include "SoftDrop.hh" // In external code, this should be fastjet/contrib/SoftDrop.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  // first get some anti-kt jets
  double R = 1.0, ptmin = 20.0;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
  
  // give the soft drop groomer a short name
  // Use a symmetry cut z > z_cut R^beta
  // By default, there is no mass-drop requirement
  double z_cut = 0.10;
  double beta  = 2.0;
  contrib::SoftDrop sd(beta, z_cut);
  cout << "SoftDrop groomer is: " << sd.description() << endl;

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    // Run SoftDrop and examine the output
    PseudoJet sd_jet = sd(jets[ijet]);
    cout << endl;
    cout << "original    jet: " << jets[ijet] << endl;
    cout << "SoftDropped jet: " << sd_jet << endl;
    
    assert(sd_jet != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
     
    cout << "  delta_R between subjets: " << sd_jet.structure_of<contrib::SoftDrop>().delta_R() << endl;
    cout << "  symmetry measure(z):     " << sd_jet.structure_of<contrib::SoftDrop>().symmetry() << endl;
    cout << "  mass drop(mu):           " << sd_jet.structure_of<contrib::SoftDrop>().mu() << endl;
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
