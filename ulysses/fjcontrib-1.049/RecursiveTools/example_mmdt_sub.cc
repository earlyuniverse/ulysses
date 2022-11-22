// $Id: example_mmdt_sub.cc 686 2014-06-14 03:25:09Z jthaler $
//
// Copyright (c) 2014, Gavin Salam
//
/// \file example_mmdt_sub.cc
///
/// An example to illustrate how to use the ModifiedMassDropTagger
/// together with pileup subtraction.
///
/// Usage:
///
/// \verbatim
///    ./example_mmdt_sub < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
/// \endverbatim
///
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
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/config.h"
#include "ModifiedMassDropTagger.hh" // In external code, this should be fastjet/contrib/ModifiedMassDropTagger.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event);
void do_analysis(const vector<PseudoJet> & jets, const Subtractor * subtractor);
ostream & operator<<(ostream &, const PseudoJet &);

// give the tagger a short name
typedef contrib::ModifiedMassDropTagger MMDT;

//----------------------------------------------------------------------
int main(){

  // Specify our basic jet tools: 
  // jet definition & area definition
  double R = 1.0, rapmax = 5.0, ghost_area = 0.01;
  int    repeat = 1;
  JetDefinition jet_def(cambridge_algorithm, R);
  AreaDefinition area_def(active_area_explicit_ghosts,
                          GhostedAreaSpec(SelectorAbsRapMax(rapmax),repeat,ghost_area));
  cout << "# " << jet_def.description() << endl;
  cout << "# " << area_def.description() << endl;

  // then our background estimator: use the (fast) grid-median method,
  // and manually include reasonable rapidity dependence for
  // particle-level 14 TeV events.
  double grid_size = 0.55;
  GridMedianBackgroundEstimator bge(rapmax, grid_size);
  BackgroundRescalingYPolynomial rescaling (1.1685397, 0, -0.0246807, 0, 5.94119e-05);
  bge.set_rescaling_class(&rescaling);
  // define a corresponding subtractor
  Subtractor subtractor(&bge);


  //----------------------------------------------------------
  // next read in input particles and get corresponding jets
  // for the hard event (no pileup) and full event (with pileup)
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);
  cout << "# read a hard event with " << hard_event.size() << " particles" ;
#if (FASTJET_VERSION_NUMBER >= 30100) // Selector.sum(..) works only for FJ >= 3.1
  cout << ", centre-of-mass energy = " << SelectorIdentity().sum(hard_event).m();
#endif
  cout << endl;
  cout << "# read a full event with " << full_event.size() << " particles" << endl;

  // then get the CS and jets for both hard and full events
  ClusterSequenceArea csa_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea csa_full(full_event, jet_def, area_def);
  vector<PseudoJet> hard_jets = SelectorNHardest(2)(csa_hard.inclusive_jets());
  vector<PseudoJet> full_jets = SelectorNHardest(2)(csa_full.inclusive_jets());
  hard_jets = sorted_by_rapidity(hard_jets);
  full_jets = sorted_by_rapidity(full_jets);

  // estimate the background (the subtractor is automatically tied to this)
  bge.set_particles(full_event); 

  // then do analyses with and without PU, and with and without subtraction
  cout << endl << "-----------------------------------------" << endl
       << "No pileup, no subtraction" << endl;
  do_analysis(hard_jets, 0);
  cout << endl << "-----------------------------------------" << endl
       << "Pileup, no subtraction" << endl;
  do_analysis(full_jets, 0);
  cout << endl << "-----------------------------------------" << endl
       << "Pileup, with subtraction" << endl;
  do_analysis(full_jets, &subtractor);
  
  return 0;
}


//----------------------------------------------------------------------
/// do a simple MMDT + filter analysis, optionally with a subtractor
void do_analysis(const vector<PseudoJet> & jets, const Subtractor * subtractor) {
  
  // use just a symmetry cut for the tagger, with no mass-drop requirement
  double z_cut = 0.10;
  MMDT tagger(z_cut);
  cout << "tagger is: " << tagger.description() << endl;
  // tell the tagger that we will subtract the input jet ourselves
  tagger.set_subtractor(subtractor);
  tagger.set_input_jet_is_subtracted(true); 


  PseudoJet jet;
  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    if (subtractor) {
      jet = (*subtractor)(jets[ijet]);
    } else {
      jet = jets[ijet];
    }
    PseudoJet tagged_jet = tagger(jet);
    cout << endl;
    cout << "original jet" << jet << endl;
    cout << "tagged   jet" << tagged_jet << endl;
    if (tagged_jet != 0) {  // get additional informaition about satisified symmetry condition
      cout << "  delta_R between subjets: " << tagged_jet.structure_of<MMDT>().delta_R() << endl;
      cout << "  symmetry measure(z):     " << tagged_jet.structure_of<MMDT>().symmetry() << endl;
      cout << "  mass drop(mu):           " << tagged_jet.structure_of<MMDT>().mu() << endl;
    }

    // then filter the jet (useful for studies at moderate pt)
    // with a dynamic Rfilt choice (as in arXiv:0802.2470)
    double Rfilt = min(0.3, tagged_jet.structure_of<MMDT>().delta_R()*0.5);
    int    nfilt = 3;
    Filter filter(Rfilt, SelectorNHardest(nfilt));
    filter.set_subtractor(subtractor);
    PseudoJet filtered_jet = filter(tagged_jet);
    cout << "filtered jet: " << filtered_jet << endl;
    cout << endl;
  }


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
