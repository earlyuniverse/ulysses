// $Id$
//
// Copyright (c) -, 
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

///////////////////////////////////////////////////////////////////////////
//
// example code illustrating the use of the SubjetCounting routines
// 
// usage: ./example < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//        ./example < ../data/single-event.dat
//
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include <cstdio>
#include "SubjetCounting.hh" // In external code, this should be fastjet/contrib/SubjetCounting.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  JetAlgorithm algorithm = antikt_algorithm;
  double jet_rad = 1.50; // jet radius for anti-Kt algorithm
  JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
  ClusterSequence clust_seq(event,jetDef);
  vector<PseudoJet> antiKt_jets  = sorted_by_pt(clust_seq.inclusive_jets());

  cout << "Event has a total of " << (int)antiKt_jets.size() << " anti-Kt jets with jet radius 1.50" << endl;
  cout << "Going to compute n_Kt and n_CA for each such jet (with two different sets of parameters";
  cout << " for each observable)" << endl;

  //----------------------------------------------------------
  // illustrate how this SubjetCounting contrib works

  double f_Kt_1 = 0.06; // parameters defining n_Kt and n_CA
  double f_Kt_2 = 0.12;
  double pt_cut_1 = 40.0;
  double pt_cut_2 = 50.0;
  double mass_cut_off_1 = 30.0;
  double mass_cut_off_2 = 50.0;
  double ycut_1 = 0.10;
  double ycut_2 = 0.15;
  double R_min = 0.15;

  SubjetCountingKt scKt1(f_Kt_1, pt_cut_1);
  SubjetCountingKt scKt2(f_Kt_2, pt_cut_2);
  SubjetCountingCA scca1(mass_cut_off_1,ycut_1,R_min, pt_cut_1);
  SubjetCountingCA scca2(mass_cut_off_2,ycut_2,R_min, pt_cut_2);

  for (int k=0; k<(int)antiKt_jets.size(); k++)
  {
  printf("n_Kt(jet %d; f_Kt = %.2f, pt_cut = %.1f GeV) = %d\n", k+1, f_Kt_1, pt_cut_1, \
         scKt1(antiKt_jets[k]));
	 //std::cout << scKt1.description() << std::endl;
  printf("n_Kt(jet %d; f_Kt = %.2f, pt_cut = %.1f GeV) = %d\n", k+1, f_Kt_2, pt_cut_2, \
         scKt2(antiKt_jets[k]));
  printf("n_CA(jet %d; mass_cut_off = %.1f, ycut = %.2f, R_min = %.2f, pt_cut = %.1f GeV) = %d\n", \
         k+1, mass_cut_off_1, ycut_1, R_min, pt_cut_1, \
         scca1(antiKt_jets[k]));
	 //std::cout << scca1.description() << std::endl;
  printf("n_CA(jet %d; mass_cut_off = %.1f, ycut = %.2f, R_min = %.2f, pt_cut = %.1f GeV) = %d\n", \
         k+1, mass_cut_off_2, ycut_2, R_min, pt_cut_2, \
         scca2(antiKt_jets[k]));
  }

  return 0;
}

// read in input particles
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
