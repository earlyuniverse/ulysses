//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_whole_event_iterative < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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


#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "IterativeConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/IterativeConstituentSubtractor.hh

#include "functions.hh"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  // This should be set up before event loop:
  double max_eta=4;   // specify the maximal pseudorapidity for the input particles. It is used for the subtraction. Particles with eta>|max_eta| are removed and not  used during the subtraction (they are not returned). The same parameter should be used for the GridMedianBackgroundEstimator as it is demonstrated in this example. If JetMedianBackgroundEstimator is used, then lower parameter should be used  (to avoid including particles outside this range).
  double max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not important for the subtraction.

  // background estimator
  GridMedianBackgroundEstimator bge_rho(max_eta,0.5); // maximal pseudo-rapidity cut is used inside ConstituentSubtraction, but in GridMedianBackgroundEstimator, the range is specified by maximal rapidity cut. Therefore, it is important to apply the same pseudo-rapidity cut also for particles used for background estimation (using function "set_particles") and also derive the rho dependence on rapidity using this max pseudo-rapidity cut to get the correct rescaling function!  

  contrib::IterativeConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  vector<double> max_distances;
  max_distances.push_back(0.15);
  max_distances.push_back(0.2);
  vector<double> alphas;
  alphas.push_back(1);
  alphas.push_back(1);
  subtractor.set_parameters(max_distances,alphas); // in this example, 2 CS corrections will be performed: 1. correction with max_distance of 0.1 and alpha of 0, 2. correction with max_distance of 0.15 and alpha of 0. After the first correction, the scalar sum of pt from remaining ghosts is evaluated and redistributed uniformly accross the event.
  subtractor.set_ghost_removal(true);  // set to true if the ghosts (proxies) which were not used in the previous CS procedure should be removed for the next CS procedure
  subtractor.set_ghost_area(0.004); // parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_max_eta(max_eta); // parameter for maximal |eta| cut. It is specified above.
  subtractor.set_background_estimator(&bge_rho); // specify the background estimator to estimate rho. You can use also specify background estimator for the mass term as an additional parameter


  /// For the treatment of massive particles, choose one of the following options (currently the default option is 3 which seems to be one of the most optimal in our s  tudies):   
/// 1. Keep the original mass and rapidity. Not recommended since jet mass is wrongly reconstructed. One can use it by adding before initialization: 
///    subtractor.set_keep_original_masses(); 
/// 2. Keep the original mass and pseudo-rapidity. Not recommended since jet mass is wrongly reconstructed. Use these two functions for this option:
///    subtractor.set_keep_original_masses();
///    subtractor.set_fix_pseudorapidity(); 
/// 3. Set all masses to zero. Keep rapidity unchanged. This is recommended and is the default option, since observed the best performance for high pt large-R jets. 
///    Nothing needs to be specified. 
/// 4. Set all masses to zero. Keep pseudo-rapidity unchanged. Also recommended, almost the same performance as for option 3. Use function:
///    subtractor.set_fix_pseudorapidity();
/// 5. Do correction of m_delta=sqrt(p_T^2+m^2)-p_t. Not recommended. One can use it by adding before initialization:
///    subtractor.set_do_mass_subtraction();
///    One must additionally provide option for the rho_m background estimation. There are several possibilities:
///    a) Same background estimator as for rho. Use this function:
///     subtractor.set_common_bge_for_rho_and_rhom();  // this must be used after set_background_estimator function.
///    b) Use a separate background estimator - you need to specify it as an additional parameter in function set_background_estimator
///    c) Use scalar background density using function set_scalar_background_density(double rho, double rhom).
/// 6. Same as 5 just keep pseudo-rapidity unchanged. Not recommended. Additionally to what is written for option 5, use:
///    subtractor.set_fix_pseudorapidity();
/// 7. Keep rapidity and pseudo-rapidity fixed (scale fourmomentum). Recommended - observed better performance than the mass correction. Use:
///    subtractor.set_scale_fourmomentum();                                        


  // Example selector for ConstituentSubtractor:                                                                                                                    
  //Selector sel_CS_correction = SelectorPhiRange(0,3) * SelectorEtaMin(-1.5) * SelectorEtaMax(0);        
  // subtractor.set_ghost_selector(&sel_CS_correction);             // uncomment in case you want to use ghost selector. Only ghosts fulfilling this selector will be constructed. No selection on input particles is done. The selection on input particles is the responsibility of the user, or the user can use the "set_particle_selector" function. However, the maximal eta cut (specified with "set_max_eta" function) is done always for both, ghosts and input particles.
  // subtractor.set_particle_selector(&sel_CS_correction);             // uncomment in case you want to use particle selector. Only particles fulfilling this selector will be corrected. The other particles will be unchanged. However, the maximal eta cut (specified with "set_max_eta" function) is done always for both, ghosts and input particles.                                                                                                                                                     
  Selector sel_max_pt = SelectorPtMax(15);
  subtractor.set_particle_selector(&sel_max_pt); // only particles with pt<15 will be corrected - the other particles will be copied without any changes.


  // If you want to use the information about hard proxies (typically charged tracks from primary vertex and/or high pt calorimeter clusters), then use the following lines:    
  /*vector<double> nearby_hard_radii;
  nearby_hard_radii.push_back(0.2);
  nearby_hard_radii.push_back(0.2);
  vector<double> nearby_hard_factors;
  nearby_hard_factors.push_back(2);
  nearby_hard_factors.push_back(5);
  subtractor.set_nearby_hard_parameters(nearby_hard_radii,nearby_hard_factors);  // In this example, during the first iteration, if there is a hard proxy within deltaR distance of 0.2, then the CS distance is multiplied by factor of 2, i.e. such particle is less favoured in the subtraction procedure. Similarly in the second iteration. If you uncomment these lines, then also uncomment line 120.
  */

  subtractor.initialize();  // this is new compared to previous usages of ConstituentSubtractor! It should be used after specifying all the parameters and before event loop.
  // print info (optional)
  cout << subtractor.description() << endl;




  // in event loop
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  vector<PseudoJet> *hard_event_charged=new vector<PseudoJet>;
  vector<PseudoJet> *background_event_charged=new vector<PseudoJet>;
  read_event(hard_event, full_event, hard_event_charged, background_event_charged);

  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with pseudo-rapidity |eta|<4" << endl;

 
  // jet definition and selector for jets
  JetDefinition jet_def(antikt_algorithm, 0.7);
  Selector sel_jets = SelectorNHardest(3) * SelectorAbsEtaMax(max_eta_jet);

  // clustering of the hard-scatter event and uncorrected event
  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());


  bge_rho.set_particles(full_event);

  // the correction of the whole event with ConstituentSubtractor
    vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);  // Do not use max_eta parameter here!!! It makes nothing and this parameter will be removed from here.
  // if you want to use the information about hard proxies, use this version: 
  //  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event,hard_event_charged);  // here all charged hard particles are used for hard proxies. In real experiment, this corresponds to charged tracks from primary vertex. Optionally, one can add among the hard proxies also high pt calorimeter clusters after some basic pileup correction.   


  // clustering of the corrected event
  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());


  ios::fmtflags f( cout.flags() );
  cout << setprecision(4) << fixed;
  cout << endl << "Corrected particles in the whole event:" << endl;
  for (unsigned int i=0; i<corrected_event.size(); i++){
    const PseudoJet &particle = corrected_event[i];
    cout << "pt = " << particle.pt()
         << ", phi = " << particle.phi()
         << ", rap = " << particle.rap()
         << ", |mass| = " << fabs(particle.m()) << endl;
  }
  cout << endl;

  // shape variables:
  //----------------------------------------------------------
  JetWidth width;

  // subtract and print the result
  //----------------------------------------------------------
  cout.flags( f );
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



