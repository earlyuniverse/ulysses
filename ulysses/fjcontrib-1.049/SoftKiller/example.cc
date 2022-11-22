//----------------------------------------------------------------------
// Example on how to use the SoftKiller
//
// run it with
//  ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//----------------------------------------------------------------------

// $Id$
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin. P. Salam and Gregory Soyez
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
#include <iomanip>

#include "fastjet/ClusterSequence.hh"
#include "SoftKiller.hh" // In external code, this should be fastjet/contrib/SoftKiller.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event);

//----------------------------------------------------------------------
int main(){
  // read in input particles
  //
  // since we use here simulated data we can split the hard event
  // from the full (i.e. with pileup added) one
  //
  // (see also example 07 in FastJet)
  //----------------------------------------------------------
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  // keep the particles up to 5 units in rapidity
  double rapmax = 5.0;
  hard_event = SelectorAbsRapMax(rapmax)(hard_event);
  full_event = SelectorAbsRapMax(rapmax)(full_event);
  
  // do the clustering and get the jets from the hard event
  // as well as for the full event without applying any
  // subtraction
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.4);

  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);

  Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);

  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // now apply the soft killerto the full event
  // Then cluster the resulting event
  //----------------------------------------------------------
  double grid_size = 0.4;
  contrib::SoftKiller soft_killer(rapmax, grid_size);
  //contrib::SoftKiller soft_killer(-0.0, 4.0, 0.5, 0.3);


  double pt_threshold;
  vector<PseudoJet> soft_killed_event;
  soft_killer.apply(full_event, soft_killed_event, pt_threshold);
  cout << "# Ran the following soft killer: " << soft_killer.description() << endl;

  // alternative, more compact invocation
  //vector<PseudoJet> soft_killed_event = soft_killer(full_event);

  ClusterSequence clust_seq_kill(soft_killed_event, jet_def);  

  vector<PseudoJet> kill_jets = sel_jets(clust_seq_kill.inclusive_jets());

  cout << setprecision(4);
  cout << "Soft Killer applied a pt threshold of " << pt_threshold << endl;

  // run things and print the result
  //----------------------------------------------------------
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
  }
  cout << endl;

  cout << "# original full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
  }
  cout << endl;

  cout << "# jets after applying the soft killer" << endl;
  for (unsigned int i=0; i<kill_jets.size(); i++){
    const PseudoJet &jet = kill_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m() << endl;
  }
  cout << endl;

  return 0;
}

//------------------------------------------------------------------------
// read the event with and without pileup
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event){
  string line;
  int  nsub  = 0; // counter to keep track of which sub-event we're reading
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {break;}
    if (line.substr(0,9) == "#SUBSTART") {
      // if more sub events follow, make copy of first one (the hard one) here
      if (nsub == 1) hard_event = full_event;
      nsub += 1;
    }
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    full_event.push_back(particle);
  }

  // if we have read in only one event, copy it across here...
  if (nsub == 1) hard_event = full_event;

  // if there was nothing in the event 
  if (nsub == 0) {
    cerr << "Error: read empty event\n";
    exit(-1);
  }

  cout << "# " << nsub-1 << " pileup events on top of the hard event" << endl;
}
