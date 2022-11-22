//----------------------------------------------------------------------
/// \file
/// \page Example01 01 - basic usage example
///
/// fastjet basic example program:
///   simplest illustration of the usage of the basic classes:
///   fastjet::PseudoJet, fastjet::JetDefinition and 
///   fastjet::ClusterSequence
///
/// run it with    : ./01-basic < data/single-event.dat
///
/// Source code: 01-basic.cc
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

/// an example program showing how to use fastjet
int main(){
  
  // read in input particles
  //----------------------------------------------------------
  vector<fastjet::PseudoJet> input_particles;
  
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fastjet::PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
  }
  

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);


  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);


  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));


  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt}
  //----------------------------------------------------------
  cout << "Ran " << jet_def.description() << endl;

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
