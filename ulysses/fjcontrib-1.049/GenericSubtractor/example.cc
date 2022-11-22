//----------------------------------------------------------------------
// Example on how to use this contribution
//
// run it with
//  ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//----------------------------------------------------------------------

// $Id: example.cc 859 2015-09-21 10:11:32Z gsalam $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include "ExampleShapes.hh"
#include "GenericSubtractor.hh"

using namespace std;
using namespace fastjet;

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

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsRapMax(4.0)(hard_event);
  full_event = SelectorAbsRapMax(4.0)(full_event);
  
  // do the clustering and get the jets
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.7);
  AreaDefinition area_def(active_area_explicit_ghosts,
                          GhostedAreaSpec(SelectorAbsRapMax(4.0)));

  ClusterSequenceArea clust_seq_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea clust_seq_full(full_event, jet_def, area_def);

  Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);

  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // create what we need for the background estimation
  //----------------------------------------------------------
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(3.0);
  ClusterSequenceArea clust_seq_rho(full_event, jet_def, area_def);  

  JetMedianBackgroundEstimator bge_rho(rho_range, jet_def_for_rho, area_def);
  bge_rho.set_particles(full_event);

  Subtractor subtractor(&bge_rho);

  // the shape part
  //----------------------------------------------------------
  contrib::Angularity shape(1.0); // angularity with alpha=1.0
  //contrib::KtDij shape; // for a test of ShapeWithPartition
  contrib::GenericSubtractor gen_sub(&bge_rho);
  contrib::GenericSubtractorInfo info;

  cout << gen_sub.description() << endl;
  cout << setprecision(4);

  // uncomment this if you also want rho_m to be estimated (using the
  // same background estimator)
  //gen_sub.use_common_bge_for_rho_and_rhom(true);

  // run things and print the result
  //----------------------------------------------------------
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", angularity = " << shape(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", angularity = " << shape(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    PseudoJet subtracted_jet = subtractor(jet);

    // compute the subtracted shape, and retrieve additional information
    double subtracted_shape = gen_sub(shape, jet, info);

    // print out the result and some extra information
    //
    // Note that, the last line is commented so that it is not taken
    // into account by 'make check' since its result may be platform
    // dependent
    cout << "pt = " << subtracted_jet.pt()
	 << ", rap = " << subtracted_jet.rap()
	 << ", angularity = " << subtracted_shape << endl;
    cout << "  rho  = " << info.rho() << endl;
    cout << "  rhom = " << info.rhom() << endl;
    cout << "  1st derivative: " << info.first_derivative() << endl;
    cout << "  2nd derivative: " << info.second_derivative() << endl;
    cout << "  unsubtracted: " << info.unsubtracted() << endl;
    cout << "  1st order: " << info.first_order_subtracted() << endl;
    cout << "# step used: " << info.ghost_scale_used() << endl;
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
