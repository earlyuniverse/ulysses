// run it with
//  ./example < single-epDIS-event.dat
//----------------------------------------------------------------------
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

#include <iostream>
#include <iomanip>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "Centauro.hh"
#include "fastjet/EECambridgePlugin.hh"


using namespace std;
using namespace fastjet;

void read_event(vector<PseudoJet> &full_event);

//----------------------------------------------------------------------
int main(){
  // example for the creation and execution of the fastjet plugin
  // for the Centauro jet algorithm
  //
  // read in input particles
  //----------------------------------------------------------
  vector<PseudoJet> full_event;
  read_event(full_event);

  // create the jet definition using the plugin mechanism
  //----------------------------------------------------------
  //The Centauro jet algorithm is described in arXiv:2006.10751, "Asymmetric jet clustering in deep-inelastic scattering", Miguel Arratia, Yiannis Makris, Duff Neill, Felix Ringer, Nobuo Sato.
  fastjet::contrib::CentauroPlugin * centauro_plugin = new fastjet::contrib::CentauroPlugin(1.0);
  fastjet::JetDefinition jet_def(centauro_plugin);
  ClusterSequence clust_seq(full_event, jet_def);

  vector<PseudoJet> jets = clust_seq.inclusive_jets(0);

  vector<PseudoJet> sortedJets = sorted_by_E(jets);
  cout << "jets in inclusive clustering " << endl;
  for (unsigned int i=0; i<sortedJets.size(); i++){
    const PseudoJet &jet = sortedJets[i];
    vector<fastjet::PseudoJet> constituents = jet.constituents();
    cout << " rap = " << jet.rap() << " e " << jet.e() << " n constituents " << constituents.size() << endl;
  }
  cout << endl;

  return 0;
}


//------------------------------------------------------------------------
// read the event
void read_event(vector<PseudoJet> &full_event){
  string line;

  std::cout << "Particles that are input" << std::endl;
  while (getline(cin, line)) {                                                                                                                                                                                         istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;                                                                                                                                                                                 PseudoJet particle(px,py,pz,E);
    std::cout << px << " " << py << " " << pz << " " << E <<std::endl;
    // push event onto back of full_event vector
    full_event.push_back(particle);                                                                                                                                                                                  }

}
