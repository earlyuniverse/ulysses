//----------------------------------------------------------------------
/// \file
/// \page Example10 10 - extracting subjets
///
/// fastjet example program to show how to access subjets;
/// 
/// See also 12-boosted_higgs.cc to see the use of subjets for
/// identifying boosted higgs (and other objects)
///
/// run it with    : ./10-subjets < data/single-event.dat
///
/// Source code: 10-subjets.cc
//----------------------------------------------------------------------

//STARTHEADER
// $Id$
//
// Copyright (c) 2005-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include "fastjet/ClusterSequence.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io

using namespace std;
using namespace fastjet;

int main(){
  
  // read in input particles
  //----------------------------------------------------------
  vector<PseudoJet> input_particles;
  
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(PseudoJet(px,py,pz,E)); 
  }
  

  // create a jet definition: 
  // for subjet studies, Cambridge/Aachen is the natural algorithm
  //----------------------------------------------------------
  double R = 1.0;
  JetDefinition jet_def(cambridge_algorithm, R);


  // run the jet clustering with the above jet definition
  // and get the jets above 5 GeV
  //----------------------------------------------------------
  ClusterSequence clust_seq(input_particles, jet_def);
  double ptmin = 6.0;
  vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

  // extract the subjets at a smaller angular scale (Rsub=0.5)
  //
  // This is done by using ClusterSequence::exclusive_subjets(dcut):
  // for the Cambridge/Aachen algorithm, running with R and then
  // asking for exclusive subjets with dcut should give the same
  // subjets as rerunning the algorithm with R'=R*sqrt(dcut) on the
  // jet's constituents.
  //
  // At the same time we output a summary of what has been done and the 
  // resulting subjets
  //----------------------------------------------------------
  double Rsub = 0.5;
  double dcut = pow(Rsub/R,2);

  // a "header" for the output
  cout << "Ran " << jet_def.description() << endl;
  cout << "Showing the jets above " << ptmin << " GeV" << endl;
  cout << "And their subjets for Rsub = " << Rsub << endl;
  printf("%10s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "n constituents");

  // show the jets and their subjets
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    // get the subjets
    vector<PseudoJet> subjets = sorted_by_pt(inclusive_jets[i].exclusive_subjets(dcut));

    cout << endl;
    // print the jet and its subjets
    printf("%5u      %15.8f %15.8f %15.8f %8d\n", i, 
	   inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp(), int(inclusive_jets[i].constituents().size()));

    for (unsigned int j=0; j<subjets.size(); j++)
      printf("   sub%4u %15.8f %15.8f %15.8f %8u\n", j, 
	     subjets[j].rap(), subjets[j].phi(),
	     subjets[j].perp(), 
	     (unsigned int) subjets[j].constituents().size());
  }

  return 0;
}
