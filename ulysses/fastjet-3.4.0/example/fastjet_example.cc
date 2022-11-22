/// \file
/// \page Examples FastJet examples
///
/// The FastJet examples have been organised by order of complexity,
/// starting by the simplest case and introducing features one after
/// another.
///   - \subpage Example01
///   - \subpage Example02
///   - \subpage Example03
///   - \subpage Example04
///   - \subpage Example05
///   - \subpage Example06
///   - \subpage Example07 (\subpage Example07old "old version")
///   - \subpage Example08
///   - \subpage Example09
///   - \subpage Example10
///   - \subpage Example11
///   - \subpage Example12 (\subpage Example12old "old version")
///   - \subpage Example13
///   - \subpage Example14

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


//----------------------------------------------------------------------
// fastjet example program. 
//
// Compile it with: make fastjet_example
// run it with    : ./fastjet_example < data/single-event.dat
//
// People who are familiar with the ktjet package are encouraged to
// compare this file to the ktjet_example.cc program which does the
// same thing in the ktjet framework.
//----------------------------------------------------------------------
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 
#include <cstdio>

using namespace std;

// a declaration of a function that pretty prints a list of jets
void print_jets (const vector<fastjet::PseudoJet> &);

/// an example program showing how to use fastjet
int main () {
  
  vector<fastjet::PseudoJet> input_particles;
  
  // Read in input particles
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fastjet::PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
  }
  
  // create an object that represents your choice of jet algorithm and 
  // the associated parameters
  double Rparam = 1.0;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, Rparam, recomb_scheme, strategy);

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  // tell the user what was done
  cout << "Ran " << jet_def.description() << endl;
  cout << "Strategy adopted by FastJet was "<<
       clust_seq.strategy_string()<<endl<<endl;

  // extract the inclusive jets with pt > 5 GeV
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);

  // print them out
  cout << "Printing inclusive jets with pt > "<< ptmin<<" GeV\n";
  cout << "---------------------------------------\n";
  print_jets(inclusive_jets);
  cout << endl;

  // extract the exclusive jets with dcut = 25 GeV^2 
  double dcut = 25.0;
  vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(dcut);

  // print them out
  cout << "Printing exclusive jets with dcut = "<< dcut<<" GeV^2\n";
  cout << "--------------------------------------------\n";
  print_jets(exclusive_jets);


}


//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const vector<fastjet::PseudoJet> & jets) {

  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", 
	 "phi", "pt", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    // the following is not super efficient since it creates an
    // intermediate constituents vector
    int n_constituents = sorted_jets[i].constituents().size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].perp(), n_constituents);
  }

}
