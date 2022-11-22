//----------------------------------------------------------------------
/// \file example.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::FlavorConePlugin class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example.cc 1045 2017-08-29 15:07:57Z philten $
//
// Copyright (c) 2017, Philip Ilten, Nicholas Rodd, Jesse Thaler,
// Michael Williams 
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
#include <sstream>
#include <iomanip>

#include "fastjet/PseudoJet.hh"
#include <sstream>
// In external code, this should be fastjet/contrib/FlavorCone.hh
#include "FlavorCone.hh"

using namespace std;
using namespace fastjet;

// Forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // Read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# Read an event with " << event.size() << " particles" << endl;

  // setting user indices to keep track of things, and then sorting by pT
  for (unsigned int ijet = 0; ijet < event.size(); ++ijet) {
    event[ijet].set_user_index(ijet);
  }
  event = sorted_by_pt(event);

  // Create the seeds. (These are not real seeds, which would be
  // heavy-flavor-tagged objects, but that information is not contained
  // in the sample event.  So we just use the two hardest pT objects)
  if (event.size() < 2) return 1;
  vector<PseudoJet> seeds;
  cout << "# Setting two seeds for clustering" << endl;
  seeds.push_back(event[0]);
  seeds.push_back(event[1]);


  // Set the jet radius
  double rcut(0.5);

  // Defining FlavorCone as a FastJet plugin
  // Note that you have to create a new jet definition for each seed choice
  contrib::FlavorConePlugin jdf(seeds, rcut);

  // Running flavor cone
  // Also note that there is no minimum pT cut applied
  cout << "# Running " << jdf.description() << endl;
  ClusterSequence jcs(event, &jdf);
  vector<PseudoJet> jets = jcs.inclusive_jets(0);
  
  // Print the FlavorCone results
  cout << "#" << setw(9) << "Seed ID"
    << setw(10) << "px"
    << setw(10) << "py"
    << setw(10) << "pz"
    << setw(10) << "e"
    << setw(10) << "const." << endl;
  
  for (unsigned int ijet = 0; ijet < jets.size(); ++ijet) {
    // Find the seed used to cluster the jet
    const contrib::FlavorConePlugin::Extras* extras = 
      dynamic_cast<const contrib::FlavorConePlugin::Extras*>(jcs.extras());
    cout << setw(10) << extras->seed(jets[ijet]).user_index();
    // Print out the jet four-vector
    for (int icmp = 0; icmp < 4; ++icmp) {
      cout << setw(10) << jets[ijet][icmp];
    }
    cout << setw(10) << jets[ijet].constituents().size();
    cout << endl;
  }

  return 0;
}

// Read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // Take substrings to avoid problems when there is extra "pollution"
    // characters (e.g. line-feed)
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // Push event onto back of full_event vector
    event.push_back(particle);
  }
}
