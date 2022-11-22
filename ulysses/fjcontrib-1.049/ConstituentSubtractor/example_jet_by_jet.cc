// $Id: example_jet_by_jet.cc 1240 2020-02-23 13:51:05Z peter.berta $
//
//----------------------------------------------------------------------
// Example on how to use this contribution
//
// run it with
//  ./example_jet_by_jet < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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
//----------------------------------------------------------------------


#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "ConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/ConstituentSubtractor.hh

#include "functions.hh"


//----------------------------------------------------------------------
int main(){

  // It is highly recommended to use the ICS method instead of the jet-by-jet CS method since the ICS method has much better performance. The jet-by-jet CS is not developed anymore and this example may be out-of-date. If you still think the jet-by-jet CS is useful for you, please, contact the authors.


  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsRapMax(4.0)(hard_event);
  full_event = SelectorAbsRapMax(4.0)(full_event);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with rapidity |y|<4" << endl;

 
  // do the clustering with ghosts and get the jets
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.7);
  double ghost_area=0.01;  // the density of ghosts can be changed through this variable
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(4.0,1,ghost_area));  // in this step, the ghosts are added among the constituents of the jets

  ClusterSequenceArea clust_seq_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea clust_seq_full(full_event, jet_def, area_def);

  Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);

  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

 // create what we need for the background estimation
  //----------------------------------------------------------
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  Selector rho_range =  SelectorAbsRapMax(3.0);
  ClusterSequenceArea clust_seq_rho(full_event, jet_def, area_def);  

  JetMedianBackgroundEstimator bge_rho(rho_range, jet_def_for_rho, area_def);
  BackgroundJetScalarPtDensity *scalarPtDensity=new BackgroundJetScalarPtDensity();
  bge_rho.set_jet_density_class(scalarPtDensity); // this changes the computation of pt of patches from vector sum to scalar sum. The scalar sum seems more reasonable.
  bge_rho.set_particles(full_event);

  // subtractor:
  //----------------------------------------------------------
  contrib::ConstituentSubtractor subtractor(&bge_rho);  

 // by default, the masses of all particles are set to zero.
  // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
  //  subtractor.set_common_bge_for_rho_and_rhom();
  cout << subtractor.description() << endl;



  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout << setprecision(4);
  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    PseudoJet subtracted_jet = subtractor(jet);
    cout << "pt = " << subtracted_jet.pt() << ", rap = " << subtracted_jet.rap() << ", mass = " << subtracted_jet.m() << ", width = " << width(subtracted_jet) << endl;
    /*    for (unsigned int iconst=0; iconst<subtracted_jet.constituents().size(); iconst++){
      fastjet::PseudoJet cons=subtracted_jet.constituents().at(iconst);
      cout << "pt = " << cons.pt() << ", phi = " << cons.phi() << ", rap = " << cons.rap() << ", |mass| = " << cons.m() << endl;
      }*/
  }



  cout << endl;

  return 0;
}



