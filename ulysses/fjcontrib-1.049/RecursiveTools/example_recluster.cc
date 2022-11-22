//----------------------------------------------------------------------
/// \file example_recluster.cc
///
/// This example program is meant to illustrate how the
/// fastjet::contrib::ModifiedMassDropTagger class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_recluster < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_recluster.cc 705 2014-07-07 14:37:03Z gsoyez $
//
// Copyright (c) 2014, Gavin P. Salam
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

#include <sstream>
#include <iomanip>
#include <cmath>
#include "fastjet/ClusterSequence.hh"

#include "Recluster.hh" // In external code, this should be fastjet/contrib/Recluster.hh

using namespace std;
using namespace fastjet;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  double R     = 1.0;
  double ptmin = 20.0;
  double Rsub  = 0.3;

  //----------------------------------------------------------
  // start with an example from anti-kt jets
  cout << "--------------------------------------------------" << endl;
  JetDefinition jet_def_akt(antikt_algorithm, R);
  ClusterSequence cs_akt(event, jet_def_akt);
  vector<PseudoJet> jets_akt = sorted_by_pt(cs_akt.inclusive_jets(ptmin));
  PseudoJet jet_akt = jets_akt[0];
  cout << "Starting from a jet obtained from: " << jet_def_akt.description() << endl
       << "  " << jet_akt << endl << endl;

  // recluster with C/A ("infinite" radius)
  contrib::Recluster recluster_ca_inf(cambridge_algorithm, JetDefinition::max_allowable_R);
  PseudoJet rec_jet_ca_inf = recluster_ca_inf(jet_akt);
  cout << "Reclustering with: " << recluster_ca_inf.description() << endl
       << "  " << rec_jet_ca_inf << endl << endl;;

  // recluster with C/A (small radius), keeping all subjets
  contrib::Recluster recluster_ca_sub(cambridge_algorithm, Rsub, false);
  PseudoJet rec_jet_ca_sub = recluster_ca_sub(jet_akt);
  cout << "Reclustering with: " << recluster_ca_sub.description() << endl
       << "  " << rec_jet_ca_sub << endl;
  vector<PseudoJet> pieces = rec_jet_ca_sub.pieces();
  cout << "   subjets: " << endl;
  for (unsigned int i=0;i<pieces.size();i++)
    cout << "    " << pieces[i] << endl;
  cout << endl;

  // recluster with kt (small radius), keeping all subjets
  contrib::Recluster recluster_kt_sub(kt_algorithm, Rsub, false);
  PseudoJet rec_jet_kt_sub = recluster_kt_sub(jet_akt);
  cout << "Reclustering with: " << recluster_kt_sub.description() << endl
       << "  " << rec_jet_kt_sub << endl;
  pieces = rec_jet_kt_sub.pieces();
  cout << "   subjets: " << endl;
  for (unsigned int i=0;i<pieces.size();i++)
    cout << "    " << pieces[i] << endl;
  cout << endl;
  
  //----------------------------------------------------------
  // now an example starting from C/A jets
  cout << "--------------------------------------------------" << endl;
  JetDefinition jet_def_ca(cambridge_algorithm, R);
  ClusterSequence cs_ca(event, jet_def_ca);
  vector<PseudoJet> jets_ca = sorted_by_pt(cs_ca.inclusive_jets(ptmin));
  PseudoJet jet_ca = jets_ca[0];
  cout << "Starting from a jet obtained from: " << jet_def_ca.description() << endl
       << "  " << jet_ca << endl << endl;

  // recluster with C/A ("infinite" radius)
  rec_jet_ca_inf = recluster_ca_inf(jet_ca);
  cout << "Reclustering with: " << recluster_ca_inf.description() << endl
       << "  " << rec_jet_ca_inf << endl << endl;

  // recluster with C/A (small radius), keeping all subjets
  rec_jet_ca_sub = recluster_ca_sub(jet_ca);
  cout << "Reclustering with: " << recluster_ca_sub.description() << endl
       << "  " << rec_jet_ca_sub << endl;
  pieces = rec_jet_ca_sub.pieces();
  cout << "   subjets: " << endl;
  for (unsigned int i=0;i<pieces.size();i++)
    cout << "    " << pieces[i] << endl;
  cout << endl;

  // recluster with kt (small radius), keeping all subjets
  rec_jet_kt_sub = recluster_kt_sub(jet_ca);
  cout << "Reclustering with: " << recluster_kt_sub.description() << endl
       << "  " << rec_jet_kt_sub << endl;
  pieces = rec_jet_kt_sub.pieces();
  cout << "   subjets: " << endl;
  for (unsigned int i=0;i<pieces.size();i++)
    cout << "    " << pieces[i] << endl;
  cout << endl;

  return 0;
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
         << " phi = " << jet.phi()
         << " ClusSeq = " << (jet.has_associated_cs() ? "yes" : "no");
  }
  return ostr;
}
