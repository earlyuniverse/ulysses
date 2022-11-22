//----------------------------------------------------------------------
/// \file
/// \page Example02 02 - changing the jet definition
///
/// fastjet basic example program:
///   illustration of the usage of how to change the jet definition
///   used for the clustering (see also fastjet::JetDefinition)
///
/// run it with    : ./02-jetdef < data/single-event.dat
///
/// Source code: 02-jetdef.cc
//----------------------------------------------------------------------
//
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

  // select a jet algorithm to use
  //
  // this could be one of
  //   {kt_algorithm, cambridge_algorithm, antikt_algorithm,
  //    genkt_algorithm, ee_kt_algorithm, ee_genkt_algorithm}
  // see example 03-plugin.cc for extra options using plugins
  // instead of naive algorithms)
  fastjet::JetAlgorithm jet_alg = fastjet::kt_algorithm;

  // when appropriate select a radius to use 
  //
  // this would not be mandatory for e^+ e^- algorithms (see 05-eplus_eminus.cc)
  double R = 0.6;

  // select an __optional__ recombination scheme
  //
  // this could be one of
  //   {E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, 
  //    BIpt2_scheme, WTA_pt_scheme, WTA_E_scheme, WTA_modp_scheme, 
  //    external_sheme}
  // 
  // Notes:
  //  - for the usage of a user-defined recombination scheme
  //    (external_scheme), see 11-boosted_higgs.cc
  //  - WTA_E_scheme, WTA_modp_scheme are meant for e+e- clusterings
  //
  // By default, the E_scheme is used 
  fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;

  // select an __optional__ strategy
  //
  // this could be chosen among
  //   {N2MinHeapTiled, N2Tiled, N2PoorTiled, N2Plain, N3Dumb, 
  //    Best, 
  //    NlnN, NlnN3pi, NlnN4pi, NlnNCam4pi, NlnNCam2pi2R, NlnNCam}
  //
  // By default, the Best strategy is chosen and we advise to keep
  // that default unless you are targeting a very specific usage. Note
  // also that the N log (N) strategies for algorithms other than
  // Cambridge/Aachen need CGAL support.
  fastjet::Strategy strategy = fastjet::Best;

  // create the JetDefinition from the above information
  fastjet::JetDefinition jet_def(jet_alg, R, recomb_scheme, strategy);


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
