//----------------------------------------------------------------------
/// \file
/// \page Example06 06 - using jet areas
///
/// fastjet example program for jet areas
/// It mostly illustrates the usage of the 
/// fastjet::AreaDefinition and fastjet::ClusterSequenceArea classes
///
/// run it with    : ./06-area < data/single-event.dat
///
/// Source code: 06-area.cc
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

#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support
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
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, R);


  // Now we also need an AreaDefinition to define the properties of the 
  // area we want
  //
  // This is made of 2 building blocks:
  //  - the area type:
  //    passive, active, active with explicit ghosts, or Voronoi area
  //  - the specifications:
  //    a VoronoiSpec or a GhostedAreaSpec for the 3 ghost-bases ones
  // 
  //---------------------------------------------------------- For
  // GhostedAreaSpec (as below), the minimal info you have to provide
  // is up to what rapidity ghosts are placed. 
  // Other commonm parameters (that mostly have an impact on the
  // precision on the area) include the number of repetitions
  // (i.e. the number of different sets of ghosts that are used) and
  // the ghost density (controlled through the ghost_area).
  // Other, more exotic, parameters (not shown here) control how ghosts
  // are placed.
  //
  // The ghost rapidity interval should be large enough to cover the
  // jets for which you want to calculate. E.g. if you want to
  // calculate the area of jets up to |y|=4, you need to put ghosts up
  // to at least 4+R (or, optionally, up to the largest particle
  // rapidity if this is smaller).
  double maxrap = 5.0;
  unsigned int n_repeat = 3; // default is 1
  double ghost_area = 0.01; // this is the default
  fastjet::GhostedAreaSpec area_spec(maxrap, n_repeat, ghost_area);

  fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);

  // run the jet clustering with the above jet and area definitions
  //
  // The only change is the usage of a ClusterSequenceArea rather than
  //a ClusterSequence
  //----------------------------------------------------------
  fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def, area_def);


  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));


  // tell the user what was done
  //  - the description of the algorithm and area used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt, number of constituents}
  //----------------------------------------------------------
  cout << endl;
  cout << "Ran " << jet_def.description() << endl;
  cout << "Area: " << area_def.description() << endl << endl;

  // label the columns
  printf("%5s %15s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "area", "area error");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n", i,
	   inclusive_jets[i].rap(), inclusive_jets[i].phi(), inclusive_jets[i].perp(),
	   inclusive_jets[i].area(), inclusive_jets[i].area_error());
  }

  return 0;
}
