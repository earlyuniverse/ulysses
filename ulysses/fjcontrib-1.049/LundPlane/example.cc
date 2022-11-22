//----------------------------------------------------------------------
/// \file example.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::LundGenerator class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example < ../data/single-event.dat
/// \endverbatim

//----------------------------------------------------------------------
// $Id: example.cc 1289 2021-11-09 11:53:53Z scyboz $
//
// Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
// Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include <string>
#include "LundGenerator.hh" // In external code, this should be fastjet/contrib/LundGenerator.hh
#include "LundJSON.hh"  // In external code, this should be fastjet/contrib/LundJSON.hh
using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  string filename = "jets.json";
  cout << "# read an event with " << event.size() << " particles" << endl;
  cout << "# writing declusterings of primary and secondary plane to file "
       << filename << endl;

  ofstream outfile;
  outfile.open(filename.c_str());
  
  // first get some anti-kt jets
  double R = 1.0, ptmin = 100.0;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));
  
  //----------------------------------------------------------
  // create an instance of LundGenerator, with default options
  contrib::LundGenerator lund;

  cout << lund.description() << endl;
  
  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    cout << endl << "Lund coordinates ( ln 1/Delta, ln kt ) of declusterings of jet " 
         << ijet << " are:" << endl;
    vector<contrib::LundDeclustering> declusts = lund(jets[ijet]);

    for (int idecl = 0; idecl < declusts.size(); idecl++) {
      pair<double,double> coords = declusts[idecl].lund_coordinates();
      cout << "(" << coords.first << ", " << coords.second << ")";
      if (idecl < declusts.size() - 1) cout << "; ";
    }

    cout << endl;
    
    // outputs the primary Lund plane
    lund_to_json(outfile, declusts); outfile << endl;
    // outputs the full Lund tree
    //to_json(cout, lund_gen, jets[ijet]); cout << endl;
  }

  cout << endl << "File " << filename << " written." << endl;



  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there is extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}
