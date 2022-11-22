//----------------------------------------------------------------------
// This example program is meant to illustrate how the
// JetFFMoments class can be used when processing events. 
// In this particular case, the calculation of the
// fragmentation function moments is done for a 'hard event' and also
// for the same event with superimposed 20 pileup events ('full
// event'), which are then subtracted from the moments.
//
// Run this example with
//  ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//
// For the hard event, the unsubtracted and subtracted moments are given in output.
// For the full event, the subtracted moments improved by the unfolding correction
// are also given.
//----------------------------------------------------------------------

// $Id: example.cc 3602 2012-09-25 13:03:36Z salam $
//
// Copyright (c) 2012-, Matteo Cacciari, Paloma Quiroga-Arias, Gavin P. Salam and Gregory Soyez
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

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <iostream>

#include "JetFFMoments.hh"

using namespace std;
using namespace fastjet;


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

  // keep the particles up to 5 units in rapidity
  hard_event = SelectorAbsRapMax(5.0)(hard_event);
  full_event = SelectorAbsRapMax(5.0)(full_event);
  
  // create what we need for the clustering and the background estimation and subtraction
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.4);
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  AreaDefinition area_def(active_area,
                          GhostedAreaSpec(SelectorAbsRapMax(5.0)));
// NB explicit ghosts do not work for moments with N < 0 with FastJet versions < 3.1			  
//  AreaDefinition area_def(active_area_explicit_ghosts,
//                          GhostedAreaSpec(SelectorAbsRapMax(5.0)));
  Selector rho_range =  SelectorDoughnut(0.4,1.2);
  JetMedianBackgroundEstimator bge(rho_range, jet_def_for_rho, area_def);
  Subtractor subtractor(&bge);

  // do the clustering
  //----------------------------------------------------------
  ClusterSequenceArea clust_seq_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea clust_seq_full(full_event, jet_def, area_def);

  // use only the two hardest jets in |y| < 4
  Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(4.0);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // fragmentation function moments
  //----------------------------------------------------------
  contrib::JetFFMoments ffms_unsubtracted(-0.5, 6.0, 14);
  contrib::JetFFMoments ffms_subtracted  (-0.5, 6.0, 14, &bge);
  contrib::JetFFMoments ffms_improved    (-0.5, 6.0, 14, &bge);
  double mu = 25.0; // typical value in the 150-200 GeV range
  ffms_improved.set_improved_subtraction(mu, rho_range,full_event,jet_def_for_rho,area_def);

  // With FastJet-3.1, this could be done this way:
  // ffms_improved.set_improved_subtraction(mu);

  // process and output hard event
  bge.set_particles(hard_event); // tell the background estimator what to use
  for (unsigned int ijet = 0; ijet<hard_jets.size(); ijet++){
    PseudoJet jet = hard_jets[ijet];
    PseudoJet subtracted_jet = subtractor(jet); // not needed for fragmentation functions, just as a reference
    vector<double> ffm_unsubtracted = ffms_unsubtracted(jet);
    vector<double> ffm_subtracted   = ffms_subtracted(jet);
    cout << "# Fragmentation function moments for hard jet " << ijet+1 
         << ": (pt,y,phi) = (" << jet.pt() << ", " 
         << jet.rap() << ", " << jet.phi() << "), #constituents="
         << jet.constituents().size() << endl;
    cout << "#                        [subtracted hard jet " << ijet+1 
         << "]:(pt,y,phi) = (" << subtracted_jet.pt() << ", " 
         << subtracted_jet.rap() << ", " << subtracted_jet.phi() << ")" << endl;
    cout << "# N  M_N(unsubtracted)  M_N(subtracted)" << endl;
    for (unsigned int in=0; in<ffm_subtracted.size(); in++){
      cout << ffms_subtracted.N(in) << " "
           << ffm_unsubtracted[in] << " "
           << ffm_subtracted[in] << endl;
    }
    cout << endl << endl;
  }

  // process and output full event
  bge.set_particles(full_event); // tell the background estimator what to use
  for (unsigned int ijet = 0; ijet<full_jets.size(); ijet++){
    PseudoJet jet = full_jets[ijet];
    PseudoJet subtracted_jet = subtractor(jet); // not needed for fragmentation functions, just as a reference
    vector<double> ffm_unsubtracted = ffms_unsubtracted(jet);
    vector<double> ffm_subtracted   = ffms_subtracted(jet);
    vector<double> ffm_improved     = ffms_improved(jet);
    cout << "# Fragmentation function moments for full jet " << ijet+1 
         << ": (pt,y,phi) = (" << jet.pt() << ", " 
         << jet.rap() << ", " << jet.phi() << "), #constituents="
         << jet.constituents().size() << endl;
    cout << "#                        [subtracted full jet " << ijet+1 
         << "]:(pt,y,phi) = (" << subtracted_jet.pt() << ", " 
         << subtracted_jet.rap() << ", " << subtracted_jet.phi() << ")" << endl;
    cout << "# N  M_N(unsubtracted)  M_N(subtracted)  M_N(improved)" << endl;
    for (unsigned int in=0; in<ffm_subtracted.size(); in++){
      cout << ffms_subtracted.N(in) << " "
           << ffm_unsubtracted[in] << " "
           << ffm_subtracted[in] << " "
           << ffm_improved[in] << endl;
    }
    cout << endl << endl;
  }
  
  cout << "# " << ffms_improved.description() << endl;

  return 0;
}
