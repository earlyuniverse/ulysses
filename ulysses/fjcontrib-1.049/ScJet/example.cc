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

#include <iostream>
#include <sstream>
#include <cstdio>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <sstream>
#include "ScJet.hh" // In external code, this should be fastjet/contrib/ScJet.hh

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

  //----------------------------------------------------------
  // illustrate how this ScJet contrib works
  JetDefinition::Plugin* scjet = new fastjet::contrib::ScJet(1.0);
  JetDefinition jet_def(scjet);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  //----------------------------------------------------------
  // print out jets (based on VariableR contrib code)
  printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity", 
         "phi", "pt","m","e", "n constituents");
  for (unsigned int i = 0; i < jets.size(); ++i) {
    PseudoJet j = jets[i];
    //int nc = cs.constituents(j).size();
    int nc = j.constituents().size();
    printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %8u\n",
           i, j.rap(), j.phi(), j.perp(), j.m(), j.e(), nc);
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
