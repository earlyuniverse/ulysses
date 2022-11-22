//----------------------------------------------------------------------
/// \file
/// \page Example04 04 - accessing clustering information in a PseudoJet
///
/// illustrate how a jet can carry information about its clustering
///
/// We do it by associating a user index to each of the input particles
/// and show what particles are in each jets (above 5 GeV)
///
/// We also illustrate a few other features about how a fastjet::PseudoJet
/// can access its underlying structure.
///
/// run it with    : ./04-constituents < data/single-event.dat
///
/// Source code: 04-constituents.cc
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
  
  valarray<double> fourvec(4);
  int index=0;
  while (cin >> fourvec[0] >> fourvec[1] >> fourvec[2] >> fourvec[3]) {
    // create a particle with the approprite 4-momentum and 
    // set its user index to keep track of its index.
    // you can construct a PseudoJet from any object that allows subscripts
    // from [0] .. [3] (the last one must be the energy)
    fastjet::PseudoJet particle(fourvec);

    particle.set_user_index(index);
    input_particles.push_back(particle); 

    index++;
  }
  

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, R);


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
  //      {index, rap, phi, pt, number of constituents}
  //----------------------------------------------------------
  cout << "Ran " << jet_def.description() << endl << endl;

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "n constituents");
  printf("        indices of constituents\n\n");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    // get the constituents of the jet
    vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();

    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp(), (unsigned int) constituents.size());

    printf("       ");
    for (unsigned int j=0; j<constituents.size(); j++){
      printf("%4u ", constituents[j].user_index());
      if (j%10==9) printf("\n       ");
    }
    printf("\n\n");
  }

  return 0;
}
