// $Id: example.cc 887 2015-10-08 08:27:29Z cspollard $
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

#include <iostream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "QCDAwarePlugin.hh" // In external code, this should be fastjet/contrib/QCDAwarePlugin.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib::QCDAwarePlugin;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

// return a flavor label by the pseudojet index
// just to make sure things work
// 40% gluons
// 20% u/ubar
// 20% d/dbar
// 20% s/sbar
int get_flavor_label(unsigned int idx) {
    int l = (idx % 10) - 3;

    // quark flavors
    // l = -3, -2, -1, 1, 2, 3
    if (abs(l) <= 3 && l != 0)
        return l;
    // gluons
    else
        return 21;
}

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how this QCDAwarePlugin contrib works

  AntiKtMeasure *akt04dm = new AntiKtMeasure(0.4);
  QCDAwarePlugin *qcdawareakt04 = new QCDAwarePlugin(akt04dm);

  cout << "using " << qcdawareakt04->description() << endl;
  ClusterSequence qcdawareakt04cs(event, qcdawareakt04);

  const vector<PseudoJet> akt04PartonJets =
      sorted_by_pt(qcdawareakt04cs.inclusive_jets());

  for (unsigned int iPJ = 0; iPJ < akt04PartonJets.size(); iPJ++) {
      PseudoJet pj = akt04PartonJets[iPJ];
      cout << "parton jet " << iPJ << " pt eta phi e label:" << endl
          << pj.pt() << " "
          << pj.eta() << " "
          << pj.phi() << " "
          << pj.e() << " "
          << pj.user_index() << endl;
  }

  delete akt04dm;
  delete qcdawareakt04;

  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  unsigned int idx = 0;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;

    PseudoJet particle(px,py,pz,E);
    particle.set_user_index(get_flavor_label(idx));

    // push event onto back of full_event vector
    event.push_back(particle);
    idx++;
  }

  return;
}
