//----------------------------------------------------------------------
/// \file
/// \page Example08 08 - using the Selector tool
///
/// fastjet sample program to illustrate the use of fastjet::Selector
///
/// run it with    : ./08-selector < data/single-event.dat
///
/// Source code: 08-selector.cc
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
#include "fastjet/Selector.hh" 
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
    // the back of the input_particles vector
    input_particles.push_back(PseudoJet(px,py,pz,E)); 
  }

  // Selector application #1: keep particles within a given acceptance
  // e.g. all particles with 1<|y|<2.5, and all particles with pt>1 for |y|<1
  // (we include the (redundant) fastjet:: prefix just as a reminder that this
  // is the namespace where Selector is to be found).
  //----------------------------------------------------------
  fastjet::Selector particle_selector = fastjet::SelectorAbsRapRange(1.0,2.5)
    || (fastjet::SelectorAbsRapMax(1.0) && fastjet::SelectorPtMin(1.0));
  cout << input_particles.size() << " particles before selector" << endl;
  input_particles = particle_selector(input_particles);
  cout << input_particles.size() << " particles after selector" << endl;

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
  JetDefinition jet_def(kt_algorithm, R);


  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  ClusterSequence clust_seq(input_particles, jet_def);


  // get the 5 hardest jets within |y|<2 
  //
  // Note that this nicely illustrates that you should watch out that
  // Selectors do not necessarily commute.
  //
  // The && operator behaves like a logical and, i.e. it keeps objects
  // that satisfy both criteria (independently). It does commute.
  //
  // The * operator applies Selectors successively (starting from the
  // rightmost one as in a usual operator product). Here, order may matter.
  //----------------------------------------------------------
  Selector jet_selector = SelectorNHardest(5) * SelectorAbsRapMax(2.0);
  vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));


  // tell the user what was done
  //  - the description of the algorithm used
  //  - the description of teh selectors used
  //  - the jets
  //    show the output as 
  //      {index, rap, phi, pt}
  //----------------------------------------------------------
  double rapmin, rapmax;

  cout << "Ran " << jet_def.description() << endl;

  cout << "Selected particles: " << particle_selector.description() << endl;
  particle_selector.get_rapidity_extent(rapmin, rapmax);
  cout << "  with a total rapidity range of [" << rapmin << ", " << rapmax << "]" << endl;

  cout << "Selected jets: " << jet_selector.description() << endl;
  jet_selector.get_rapidity_extent(rapmin, rapmax);
  cout << "  with a total rapidity range of [" << rapmin << ", " << rapmax << "]" << endl;

  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp());
  }

  return 0;
}
