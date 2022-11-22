//----------------------------------------------------------------------
/// \file
/// \page Example13 13 - boosted top tagging
///
/// fastjet example program, illustration of carrying out boosted
/// top subjet ID analysis using the Johns Hopkins top tagger as
/// introduced in arXiv:0806.0848 (Kaplan, Rehermann, Schwartz
/// and Tweedie)
///
/// run it with    : ./13-boosted_top < data/boosted_top_event.dat
///
/// Source code: 13-boosted_top.cc
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

#include <iostream> // needed for io
#include <sstream>  // needed for internal io
#include <iomanip>  
#include <cmath>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>

using namespace std;
using namespace fastjet;


//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
// core of the program
//----------------------------------------------------------------------
int main(){

  vector<PseudoJet> particles;

  // read in data in format px py pz E b-tag [last of these is optional]
  // lines starting with "#" are considered as comments and discarded
  //----------------------------------------------------------
  string line;
  while (getline(cin,line)) {
    if (line.substr(0,1) == "#") {continue;}
    istringstream linestream(line);
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;

    // construct the particle
    particles.push_back(PseudoJet(px,py,pz,E));
  }

  // compute the parameters to be used through the analysis
  // ----------------------------------------------------------
  double Et=0;
  for (unsigned int i=0; i<particles.size(); i++)
    Et += particles[i].perp();

  double R, delta_p, delta_r;
  if      (Et>2600){ R=0.4; delta_p=0.05; delta_r=0.19;}
  else if (Et>1600){ R=0.6; delta_p=0.05; delta_r=0.19;}
  else if (Et>1000){ R=0.8; delta_p=0.10; delta_r=0.19;}
  else{ cerr << "Et has to be at least 1 TeV"<< endl; return 1;}

  double ptmin = min(500.0, 0.7*Et/2);

  // find the jets
  // ----------------------------------------------------------
  JetDefinition jet_def(cambridge_algorithm, R);
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  cout << "Ran: " << jet_def.description() << endl << endl;
  cout << "2 Hardest jets: " << jets[0] << endl
       << "                " << jets[1] << endl << endl;

  if (jets[0].perp()<ptmin){
    cout << "No jet above the ptmin threshold" << endl;
    return 2;
  }

  // now do jet tagging using the Johns Hopkins top tagger
  // For simplicity, we just apply it to the hardest jet.
  //
  // In addition to delta_p and delta_r, note that there are two
  // further parameters to the JH top tagger that here are implicitly
  // set to their defaults:
  //
  // - cos_theta_W_max (defaults to 0.7) 
  // - mW (defaults to 80.4). 
  //
  // The value for mW implicitly assumes that momenta are passed in
  // GeV.
  // ----------------------------------------------------------
  JHTopTagger top_tagger(delta_p, delta_r);
  top_tagger.set_top_selector(SelectorMassRange(150,200));
  top_tagger.set_W_selector  (SelectorMassRange( 65, 95));

  PseudoJet tagged = top_tagger(jets[0]);

  cout << "Ran the following top tagger: " << top_tagger.description() << endl;

  if (tagged == 0){
    cout << "No top substructure found" << endl;
    return 0;
  }

  cout << "Found top substructure from the hardest jet:" << endl;
  cout << "  top candidate:     " << tagged << endl;
  cout << "  |_ W   candidate:  " << tagged.structure_of<JHTopTagger>().W() << endl;
  cout << "  |  |_  W subjet 1: " << tagged.structure_of<JHTopTagger>().W1() << endl;
  cout << "  |  |_  W subjet 2: " << tagged.structure_of<JHTopTagger>().W2() << endl;
  cout << "  |  cos(theta_W) =  " << tagged.structure_of<JHTopTagger>().cos_theta_W() << endl;
  cout << "  |_ non-W subjet:   " << tagged.structure_of<JHTopTagger>().non_W() << endl;
}


//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "pt, y, phi =" << setprecision(6)
       << " " << setw(9) << jet.perp() 
       << " " << setw(9)  <<  jet.rap()  
       << " " << setw(9)  <<  jet.phi()  
       << ", mass = " << setw(9) << jet.m();
  return ostr;
}
