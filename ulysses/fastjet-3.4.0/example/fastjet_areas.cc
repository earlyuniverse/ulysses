
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
// fastjet areas example program. 
//
// Compile it with: make fastjet_areas
// run it with    : ./fastjet_areas < data/single-event.dat
//
//----------------------------------------------------------------------
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"

// get info on how fastjet was configured
#include "fastjet/config.h"

#ifdef ENABLE_PLUGIN_SISCONE
#include "fastjet/SISConePlugin.hh"
#endif

#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 

using namespace std;

// a declaration of a function that pretty prints a list of jets
void print_jets (const vector<fastjet::PseudoJet> &);

/// an example program showing how to use fastjet
int main () {
  
  vector<fastjet::PseudoJet> input_particles;
  
  // read in input particles
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a fastjet::PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
  }

  // create an object that represents your choice of jet algorithm, and 
  // the associated parameters
  double Rparam = 1.0;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::kt_algorithm, Rparam, fastjet::E_scheme, strategy);
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, Rparam, strategy);
  //fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, strategy);
  //fastjet::JetDefinition jet_def(new fastjet::SISConePlugin(Rparam,0.75));

  // create an object that specifies how we to define the area
  fastjet::AreaDefinition area_def;
  bool use_voronoi = false;
  if (!use_voronoi) {
    double ghost_etamax = 6.0;
    double ghost_area    = 0.01;
    int    active_area_repeats = 1;

    // now create the object that holds info about ghosts, and from that
    // get an area definition
    fastjet::GhostedAreaSpec ghost_spec(ghost_etamax, active_area_repeats, 
                                        ghost_area);
    area_def = fastjet::AreaDefinition(fastjet::active_area,ghost_spec);
    //area_def = fastjet::AreaDefinition(fastjet::passive_area,ghost_spec);
  } else {
    double effective_Rfact = 1.0;
    area_def = fastjet::VoronoiAreaSpec(effective_Rfact);
  }

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequenceArea clust_seq(input_particles, 
                                             jet_def, area_def);
  // you can also run the individual area classes directly
  //fastjet::ClusterSequencePassiveArea clust_seq(input_particles, jet_def, 
  //                                              area_def.ghost_spec());

  // you may want to find out how much area in a given range (|y|<range)
  // is empty of real jets (or corresponds to pure "ghost" jets).
  //double range = 4.0;
  //cout << clust_seq.empty_area(range) << endl;
  //cout << clust_seq.n_empty_jets(range) << endl;

  // tell the user what was done
  cout << "Jet definition was: " << jet_def.description() << endl;
  cout << "Area definition was: " << area_def.description() << endl;
  cout << "Strategy adopted by FastJet was "<<
       clust_seq.strategy_string()<<endl<<endl;

  // extract the inclusive jets with pt > 5 GeV, sorted by pt
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);

  // print them out
  cout << "Printing inclusive jets with pt > "<< ptmin<<" GeV\n";
  cout << "---------------------------------------\n";
  print_jets(inclusive_jets);
  cout << endl;

  
  cout << "Number of unclustered particles: " 
       << clust_seq.unclustered_particles().size() << endl;


}


//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const vector<fastjet::PseudoJet> & unsorted_jets) {

  // sort jets into increasing pt
  vector<fastjet::PseudoJet> jets = sorted_by_pt(unsorted_jets);  

  printf(" ijet   rap      phi        Pt         area  +-   err\n");
  for (unsigned int j = 0; j < jets.size(); j++) {

    double area       = jets[j].area();
    double area_error = jets[j].area_error();

    printf("%5u %9.5f %8.5f %10.3f %8.3f +- %6.3f\n",j,jets[j].rap(),
	   jets[j].phi(),jets[j].perp(), area, area_error);
  }


}







