//----------------------------------------------------------------------
/// \file
/// \page Example03 03 - using plugins
///
/// fastjet plugins example program:
///   we illustrate the plugin usage
///   here, we use the SISCone plugin though different choices are possible
///   see the output of 'fastjet-config --list-plugins' for more details
///
/// Note that when using plugins, the code needs to be linked against
/// the libfastjetplugins library (with the default monolithic
/// build. For non-monolithic build, individual libraries have to be
/// used for each plugin). 
/// This is ensured in practice by calling
///   fastjet-config --libs --plugins
///
/// run it with    : ./03-plugin < data/single-event.dat
///
/// Source code: 03-plugin.cc
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

// include the SISCone plugin header if enabled
#include "fastjet/config.h"
#ifdef FASTJET_ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#else
#warning "SISCone plugin not enabled. Skipping the example"
#endif // FASTJET_ENABLE_PLUGIN_SISCONE


using namespace std;

int main(){

#ifdef FASTJET_ENABLE_PLUGIN_SISCONE

  // read in input particles
  //----------------------------------------------------------
  vector<fastjet::PseudoJet> input_particles;
  
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fastjet::PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
  }
  

  // create a jet definition fron a plugin.
  // It basically requires declaring a JetDefinition from a pointer to
  // the plugin
  //
  // we will use the SISCone plugin here. Its (mandatory) parameters
  // are a cone radius and an overlap threshold, plus other optional
  // parameters
  //
  // for other plugin, see individual documentations for a description
  // of their parameters
  //
  // the list of available plugins for a given build of FastJet can be
  // obtained using
  //   fastjet-config --list-plugins 
  // from the command line.
  //----------------------------------------------------------
  double cone_radius = 0.7;
  double overlap_threshold = 0.75;
  fastjet::SISConePlugin siscone(cone_radius, overlap_threshold);
  fastjet::JetDefinition jet_def(& siscone);


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

#endif

  return 0;

}
