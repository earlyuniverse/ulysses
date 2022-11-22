//----------------------------------------------------------------------
/// \file example_dpsi_collinear.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::RecursiveLundEEGenerator class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_dpsi_collinear < ../data/single-ee-event.dat
/// \endverbatim

//----------------------------------------------------------------------
// $Id: example_dpsi_collinear.cc 1292 2021-11-09 11:55:44Z scyboz $
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
#include "fastjet/EECambridgePlugin.hh"
#include <string>
#include "RecursiveLundEEGenerator.hh" // In external code, this should be fastjet/contrib/RecursiveLundEEGenerator.hh
using namespace std;
using namespace fastjet;

// Definitions for the collinear observable
double z1_cut = 0.1;
double z2_cut = 0.1;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);

  cout << "# read an event with " << event.size() << " particles" << endl;
  //----------------------------------------------------------
  // create an instance of RecursiveLundEEGenerator, with default options
  int depth = -1;
  bool dynamic_psi_reference = true;
  fastjet::contrib::RecursiveLundEEGenerator lund(depth, dynamic_psi_reference);
  
  // first get some C/A jets
  double y3_cut = 1.0;
  JetDefinition::Plugin* ee_plugin = new EECambridgePlugin(y3_cut);
  JetDefinition jet_def(ee_plugin);
  ClusterSequence cs(event, jet_def);

  // Get the list of primary declusterings
  const vector<contrib::LundEEDeclustering> declusts = lund.result(cs);
  
  // Find the highest-kt primary declustering (ordered in kt by default)
  int i_primary = -1;
  double psi_1;
  for (int i=0; i<declusts.size(); i++) {
    if (declusts[i].depth() == 0 && declusts[i].z() > z1_cut) {
        i_primary = i;
        psi_1 = declusts[i].psibar();
        break;
    }
  }
  if (i_primary < 0) return 0;

  // Find the highest-kt secondary associated to that Lund leaf
  int iplane_to_follow = declusts[i_primary].leaf_iplane();
  vector<const contrib::LundEEDeclustering *> secondaries;
  for (const auto & declust: declusts){
      if (declust.iplane() == iplane_to_follow) secondaries.push_back(&declust);
  }
  if(secondaries.size() < 1) return 0;

  int i_secondary = -1;
  double dpsi;
  for (int i=0; i<secondaries.size(); i++) {
      if (secondaries[i]->z() > z2_cut) {
          i_secondary = i;
          double psi_2 = secondaries[i]->psibar();
          dpsi = contrib::map_to_pi(psi_2-psi_1);
          break;
      }
  }
  if (i_secondary < 0) return 0;

  cout << "Primary that passes the zcut of "
       << z1_cut << ", with Lund coordinates  ( ln 1/Delta, ln kt, psibar ):" << endl;
  pair<double,double> coords = declusts[i_primary].lund_coordinates();
  cout << "index [" << i_primary << "](" << coords.first << ", " << coords.second << ", "
       << declusts[i_primary].psibar() << ")" << endl;
  cout << endl << "with Lund coordinates for the secondary plane that passes the zcut of "
       << z2_cut << endl;
  coords = secondaries[i_secondary]->lund_coordinates();
  cout << "index [" << i_secondary << "](" << coords.first << ", " << coords.second << ", "
       << secondaries[i_secondary]->psibar() << ")";

  cout << " --> delta_psi = " << dpsi << endl;

  // Delete the EECambridge plugin
  delete ee_plugin;

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
