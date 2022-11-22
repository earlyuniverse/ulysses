//----------------------------------------------------------------------
/// \file example_mmdt.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::ModifiedMassDropTagger class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_mmdt < ../data/single-event.dat
/// \endverbatim
///
/// It also shows operation in conjunction with a Filter.
//----------------------------------------------------------------------

// $Id: example_mmdt.cc 686 2014-06-14 03:25:09Z jthaler $
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
#include "fastjet/tools/Filter.hh"
#include "ModifiedMassDropTagger.hh" // In external code, this should be fastjet/contrib/ModifiedMassDropTagger.hh

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

  // first get some Cambridge/Aachen jets
  double R = 1.0, ptmin = 20.0;
  JetDefinition jet_def(cambridge_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
  
  // give the tagger a short name
  typedef contrib::ModifiedMassDropTagger MMDT;
  // use just a symmetry cut, with no mass-drop requirement
  double z_cut = 0.10;
  MMDT tagger(z_cut);
  cout << "tagger is: " << tagger.description() << endl;

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    // first run MMDT and examine the output
    PseudoJet tagged_jet = tagger(jets[ijet]);
    cout << endl;
    cout << "original jet: " << jets[ijet] << endl;
    cout << "tagged   jet: " << tagged_jet << endl;
    if (tagged_jet == 0) continue;  // If symmetry condition not satisified, jet is not tagged
    cout << "  delta_R between subjets: " << tagged_jet.structure_of<MMDT>().delta_R() << endl;
    cout << "  symmetry measure(z):     " << tagged_jet.structure_of<MMDT>().symmetry() << endl;
    cout << "  mass drop(mu):           " << tagged_jet.structure_of<MMDT>().mu() << endl;

    // then filter the jet (useful for studies at moderate pt)
    // with a dynamic Rfilt choice (as in arXiv:0802.2470)
    double Rfilt = min(0.3, tagged_jet.structure_of<MMDT>().delta_R()*0.5);
    int    nfilt = 3;
    Filter filter(Rfilt, SelectorNHardest(nfilt));
    PseudoJet filtered_jet = filter(tagged_jet);
    cout << "filtered jet: " << filtered_jet << endl;
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
