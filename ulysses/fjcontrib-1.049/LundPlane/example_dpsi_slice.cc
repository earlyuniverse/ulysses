//----------------------------------------------------------------------
/// \file example_dpsi_slice.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::RecursiveLundEEGenerator class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_dpsi_slice < ../data/single-ee-event.dat
/// \endverbatim

//----------------------------------------------------------------------
// $Id: example_dpsi_slice.cc 1292 2021-11-09 11:55:44Z scyboz $
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

// Definitions for the slice observable (ymax = 1, zcut = 0.1)
double abs_rap_slice = 1;
double z2_cut = 0.1;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

// returns true if PseudoJet p is in the central slice of half-width abs_rap_slice
// w.r.t. the event axis pref
bool in_slice(const PseudoJet & p, const PseudoJet & pref) {
  
  // Rotate p, and use the rapidity to determine if p is in the slice
  contrib::Matrix3 rotmat = contrib::Matrix3::from_direction(pref).transpose();
  const PseudoJet p_rot = rotmat*p;

  return fabs(p_rot.rap()) < abs_rap_slice;
}

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

  // Get the event axis (to be used for the slice)
  std::vector<PseudoJet> excl = cs.exclusive_jets(2);
  assert(excl.size() == 2);
  // order the two jets according to momentum along z axis
  if (excl[0].pz() < excl[1].pz()) {
    std::swap(excl[0],excl[1]);
  }
  const PseudoJet ev_axis = excl[0]-excl[1];

  // Get the list of primary declusterings
  const vector<contrib::LundEEDeclustering> declusts = lund.result(cs);

  // find the declustering that throws something into a fixed slice
  // (from a parent that was not in the slice) and that throws the
  // largest pt (relative to the z axis) into that slice

  int index_of_max_pt_in_slice = -1;
  double max_pt_in_slice = 0.0;
  double psi_1;
  for (unsigned i = 0; i < declusts.size(); i++) {
    const auto & decl = declusts[i];
    if (!in_slice(decl.harder(), ev_axis) && in_slice(decl.softer(), ev_axis)) {
      if (decl.softer().pt() > max_pt_in_slice) {
        index_of_max_pt_in_slice = i;
        max_pt_in_slice = decl.softer().pt();
        psi_1 = decl.psibar();
      }
    }
  }
  if (index_of_max_pt_in_slice < 0) return 0;

  // establish what we need to follow and the reference psi_1
  int iplane_to_follow = declusts[index_of_max_pt_in_slice].leaf_iplane();
  
  vector<const contrib::LundEEDeclustering *> secondaries;
  for (const auto & declust: declusts){
    if (declust.iplane() == iplane_to_follow) secondaries.push_back(&declust);
  }

  int index_of_max_kt_secondary = -1;
  double dpsi;
  for (uint i_secondary=0; i_secondary<secondaries.size(); i_secondary++) {
    if (secondaries[i_secondary]->z() > z2_cut) {

      index_of_max_kt_secondary = i_secondary;
      double psi_2 = secondaries[i_secondary]->psibar();
      dpsi = contrib::map_to_pi(psi_2 - psi_1);

      break;
    }
  }
  if (index_of_max_kt_secondary < 0) return 0;

  cout << "Primary in the central slice, with Lund coordinates ( ln 1/Delta, ln kt, psibar ):" << endl;
  pair<double,double> coords = declusts[index_of_max_pt_in_slice].lund_coordinates();
  cout << "index [" << index_of_max_pt_in_slice << "](" << coords.first << ", " << coords.second << ", "
       << declusts[index_of_max_pt_in_slice].psibar() << ")" << endl;
  cout << endl << "with Lund coordinates for the (highest-kT) secondary plane that passes the zcut of "
       << z2_cut << endl;
  coords = secondaries[index_of_max_kt_secondary]->lund_coordinates();
  cout << "index [" << index_of_max_kt_secondary << "](" << coords.first << ", " << coords.second << ", "
       << secondaries[index_of_max_kt_secondary]->psibar() << ")";

  cout << " --> delta_psi,slice = " << dpsi << endl;

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
