//----------------------------------------------------------------------
/// \file
/// \page Example05 05 - using e+e- algorithms
///
/// illustrate the use of e^+ e^- algorithms
///
/// They mostly differ from the pp algorithm by the fact that rather
/// than using a radius parameter and inclusive jets, they use
/// exclusive jets in one of the following ways:
///  - a fixed number of them
///  - with a dcut
///  - with a ycut
///
/// Note that natively, FastJet includes the kt (ee_kt_algorithm) and
/// genkt (ee_genkt_algorithm) algorithms. Others (like Cambridge for
/// e+ e-, Jade or SISCone in spherical coordinates) are available as
/// plugins (see 03-plugin.cc)
///
/// run it with    : ./05-eplus_eminus < data/single-ee-event.dat
///
/// Source code: 05-eplus_eminus.cc
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
  

  // create a jet definition for the kt algorithm (note that one
  // should not specify an R value here)
  //----------------------------------------------------------
  fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm);

  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  // get 3 exclusive jets
  //----------------------------------------------------------
  int n = 3;
  vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(n);


  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt, number of constituents}
  //----------------------------------------------------------
  cout << "Ran " << jet_def.description() << endl;

  // label the columns
  printf("%5s %15s\n","jet #", "E");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < exclusive_jets.size(); i++) {
    printf("%5u %15.8f\n",
	   i, exclusive_jets[i].perp());
  }

  return 0;
}
