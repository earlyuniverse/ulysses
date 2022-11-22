// $Id: example_whole_event_using_charged_info.cc 1240 2020-02-23 13:51:05Z peter.berta $
//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_whole_event < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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


#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequence.hh"
#include "ConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/ConstituentSubtractor.hh                                     

#include "functions.hh"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  // set up before event loop:
  contrib::ConstituentSubtractor subtractor; // no need to provide background estimator in this case
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There  are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta  
  subtractor.set_max_distance(0.3); // free parameter for the maximal allowed distance between particle i and ghost k
  subtractor.set_alpha(1);  // free parameter for the distance measure (the exponent of particle pt). The larger the parameter alpha, the more are favoured the lower pt particles in the subtraction process
  subtractor.set_ghost_area(0.01); // free parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_do_mass_subtraction();   // use this line if also the mass term sqrt(pT^2+m^2)-pT should be corrected or not. It is necessary to specify it like this because the function set_common_bge_for_rho_and_rhom cannot be used in this case.
  double CBS=1.0;  // choose the scale for scaling the background charged particles
  double CSS=1.0;  // choose the scale for scaling the signal charged particles
  subtractor.set_remove_particles_with_zero_pt_and_mass(false); // set to false if you want to have also the zero pt and mtMinuspt particles in the output. Set to true, if not. The choice has no effect on the performance. By deafult, these particles are removed - this is the recommended way since then the output contains much less particles, and therefore the next step (e.g. clustering) is faster. In this example, it is set to false to make sure that this test is successful on all systems (mac, linux).
  subtractor.set_grid_size_background_estimator(0.6); // set the grid size (not area) for the background estimation with GridMedianBackgroundEstimation which is used within CS correction using charged info 
  cout << subtractor.description() << endl;



  // EVENT LOOP

  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  vector<PseudoJet> *hard_event_charged=new vector<PseudoJet>;
  vector<PseudoJet> *background_event_charged=new vector<PseudoJet>;

  read_event(hard_event, full_event, hard_event_charged, background_event_charged);

  double max_eta=4;  // specify the maximal pseudorapidity for the particles used in the subtraction
  double max_eta_jet=3; // the maximal pseudorapidity for selected jets

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);
  *hard_event_charged = SelectorAbsEtaMax(max_eta)(*hard_event_charged);
  *background_event_charged = SelectorAbsEtaMax(max_eta)(*background_event_charged);

  cout << "# read an event with " << hard_event.size() << " signal particles, " << full_event.size() - hard_event.size() << " background particles, " << hard_event_charged->size() << " signal charged particles, and " << background_event_charged->size() << " background charged particles" << " with pseudo-rapidity |eta|<4" << endl;

 
  // do the clustering with ghosts and get the jets
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.7);
  Selector sel_jets = SelectorNHardest(3) * SelectorAbsEtaMax(max_eta_jet);

  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // correction of the whole event using charged info
  vector<PseudoJet> corrected_event=subtractor.subtract_event_using_charged_info(full_event,CBS,*background_event_charged,CSS,*hard_event_charged,max_eta);

  ios::fmtflags f( cout.flags() );
  cout << setprecision(7) << fixed;
  cout << endl << "Corrected particles in the whole event:" << endl;
  for (unsigned int i=0; i<corrected_event.size(); i++){
    const PseudoJet &particle = corrected_event[i];
    cout << "pt = " << particle.pt()
	 << ", phi = " << particle.phi()
	 << ", rap = " << particle.rap()
	 << ", |mass| = " << fabs(particle.m()) << endl;
  }
  cout << endl;


  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
  cout << setprecision(7);
  cout << "Original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "Unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
  }
  cout << endl;

  cout << "Subtracted full jets" << endl;
  for (unsigned int i=0; i<corrected_jets.size(); i++){
    const PseudoJet &jet = corrected_jets[i];

    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", mass = " << jet.m()
	 << ", width = " << width(jet) << endl;
	 }
  cout << endl;

  return 0;
}



