//----------------------------------------------------------------------
/// \file example_recursive_softdrop.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::RecursiveSoftDrop class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_recursive_softdrop < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_recursive_softdrop.cc 1074 2017-09-18 15:15:20Z gsoyez $
//
// Copyright (c) 2017-, Gavin P. Salam, Gregory Soyez, Jesse Thaler,
// Kevin Zhou, Frederic Dreyer
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
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "RecursiveSoftDrop.hh" // In external code, this should be fastjet/contrib/RecursiveSoftDrop.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

void print_prongs_with_clustering_info(const PseudoJet &jet, const string &pprefix);
void print_raw_prongs(const PseudoJet &jet);

ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  // first get some anti-kt jets
  double R = 1.0, ptmin = 100.0;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

  //----------------------------------------------------------------------
  // give the soft drop groomer a short name
  // Use a symmetry cut z > z_cut R^beta
  // By default, there is no mass-drop requirement
  double z_cut = 0.2;
  double beta  = 0.5;
  int n=4; // number of layers (-1 <> infinite)
  contrib::RecursiveSoftDrop rsd(beta, z_cut, n, R);

  // keep addittional structure info (used below)
  rsd.set_verbose_structure(true);

  // (optionally) use the same-depth variant
  //
  // instead of recursing into the largest Delta R branch until "n+1"
  // branches hav ebeen found, the same-depth variant recurses n times
  // into all the branches found in the previous iteration
  //
  //rsd.set_fixed_depth_mode();

  // (optionally) use a dynamical R0
  //
  // Instead of being normalised by the initial jet radios R0, angles
  // are notrmalised by the delta R of the previous iteration
  //
  rsd.set_dynamical_R0();

  // (optionally) recurse only in the hardest branch
  //
  // Instead of recursing into both branches found by the previous
  // iteration, only keep recursing into the hardest one
  //
  //rsd.set_hardest_branch_only();

  
  //----------------------------------------------------------------------
  cout << "RecursiveSoftDrop groomer is: " << rsd.description() << endl;

  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
    // Run SoftDrop and examine the output
    PseudoJet rsd_jet = rsd(jets[ijet]);
    cout << endl;
    cout << "original             jet: " << jets[ijet] << endl;
    cout << "RecursiveSoftDropped jet: " << rsd_jet << endl;
    
    assert(rsd_jet != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet

    // print the prong structure of the jet
    //
    // This can be done in 2 ways:
    //
    //  - either keeping the clustering information and get the
    //    branches as a succession of 2->1 recombinations (this is
    //    done calling "pieces" recursively)
    cout << endl
         << "Prongs with clustering information" << endl
         << "----------------------------------" << endl;
    print_prongs_with_clustering_info(rsd_jet, " ");
    //
    //  - or getting all the branches in a single go (done directly
    //    through the jet associated structure)
    cout << endl
         << "Prongs without clustering information" << endl
         << "-------------------------------------" << endl;
    print_raw_prongs(rsd_jet);

    cout << "Groomed prongs information:" << endl;
    cout << "index            zg        thetag" << endl;
    vector<pair<double, double> > ztg = rsd_jet.structure_of<contrib::RecursiveSoftDrop>().sorted_zg_and_thetag();
    for (unsigned int i=0; i<ztg.size();++i)
      cout << setw(5) << i+1
           << setw(14) << ztg[i].first << setw(14) << ztg[i].second << endl;
    
  }

  return 0;
}

//----------------------------------------------------------------------
// print the prongs inside the jet, showing the clustering info
void print_prongs_with_clustering_info(const PseudoJet &jet, const string &prefix){
  if (prefix.size() == 1){
    cout << " " << setw(14) << " "
         << setw(8) << "branch" << setw(14) << "branch"
         << setw(10) << "N_groomed"
         << setw(11) << "max loc"
         << setw(22) << "substructure" << endl;
    cout << " " << setw(14) << " "
         << setw(8) << "pt" << setw(14) << "mass"
         << setw(5) << "loc"
         << setw(5) << "tot"
         << setw(11) << "zdrop"
         << setw(11) << "zg"
         << setw(11) << "thetag"<< endl;
  }
  const contrib::RecursiveSoftDrop::StructureType &structure = jet.structure_of<contrib::RecursiveSoftDrop>();
  double dR = structure.delta_R();
  cout << " " << left << setw(14) << (prefix.substr(0, prefix.size()-1)+"+--> ") << right
       << setw(8) << jet.pt() << setw(14) << jet.m()
       << setw(5) << structure.dropped_count(false)
       << setw(5) << structure.dropped_count()
       << setw(11) << structure.max_dropped_symmetry(false);
  
  if (structure.has_substructure()){
    cout << setw(11) << structure.symmetry()
         << setw(11) << structure.delta_R();
  }
  cout << endl;
  
  if (dR>=0){
    vector<PseudoJet> pieces = jet.pieces();
    assert(pieces.size()==2);
    print_prongs_with_clustering_info(pieces[0], prefix+" |");
    print_prongs_with_clustering_info(pieces[1], prefix+"  ");
  }    
}

//----------------------------------------------------------------------
// print all the prongs inside the jet (no clustering info)
void print_raw_prongs(const PseudoJet &jet){
  cout << "(Raw) list of prongs:" << endl;
  if (!jet.has_structure_of<contrib::RecursiveSoftDrop>()){
    cout << "  None (bad structure)" << endl;
    return;
  }
  
  cout << setw(5) << " " << setw(11) << "pt" << setw(14) << "mass" << endl;

  vector<PseudoJet> prongs = contrib::recursive_soft_drop_prongs(jet);
  for (unsigned int iprong=0; iprong<prongs.size(); ++iprong){
    const PseudoJet & prong = prongs[iprong];
    const contrib::RecursiveSoftDrop::StructureType &structure = prong.structure_of<contrib::RecursiveSoftDrop>();
    cout << setw(5) << iprong << setw(11) << prong.pt() << setw(14) << prong.m() << endl;
  
    assert(!structure.has_substructure());
  }
  cout << endl;
}

//----------------------------------------------------------------------
/// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
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

//----------------------------------------------------------------------
/// overloaded jet info output
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi();
  }
  return ostr;
}
