//
//----------------------------------------------------------------------
// Example on how to do pileup correction on the whole event
//
// run it with
//  ./example_background_rescaling < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
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
#include "IterativeConstituentSubtractor.hh" // In external code, this should be fastjet/contrib/IterativeConstituentSubtractor.hh
#include "RescalingClasses.hh" // In external code, this should be fastjet/contrib/RescalingClasses.hh
 

#include "functions.hh"
//#include "TH1D.h"
//#include "TH2D.h"
//#include "TF1.h"


using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
int main(){
  /// first several examples for rescaling classes are listed. Then the non-commented example is used with the Iterative CS.


  ///**** rescaling in heavy ion events using rapidity dependence from TH1 histogram: ****
  /// Set the five parameters for rescaling using function
  /// BackgroundRescalingYFromRootPhi(double v2, double v3, double v4, double psi, T* hist);
  /// which is parametrized as
  ///  f(y,phi) = phi_term(phi) * rap_term(y)                                                                                                           
  /// where                                                                                                                                             
  ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))                                     
  ///  with four parameters v2, v3, v4, and psi.                                                                                                        
  /// rap_term(y) = bin content of the histogram at position y         
  ///
  /// You need to set the parameters event-by-event. Example where the TH1 histogram is filled with a random distribution. Do not use this histogram!!! Fill your own histogram with all particles (topoclusters) weighted by their pt using minimum bias events.
  /*TH1D* hist=new TH1D("hist","",100,-5,5);
  TF1 polynom("polynom","1+0.03*x*x-0.003*x*x*x*x",-5,5);
  for (int i=1;i<hist->GetNbinsX()+1;i++){
    hist->SetBinContent(i,polynom(hist->GetBinCenter(i)));
  }
  contrib::BackgroundRescalingYFromRootPhi<TH1D> rescaling(0.1,0.1,0.001,0,hist);
  rescaling.use_rap_term(true);    // this is useful to make sure that the histogram with rapidity dependence is defined. Set to false, if you do not want to use the rapidity rescaling
  */


  ///**** rescaling in heavy ion events using rapidity dependence stored in vector<double>: ****
  /// Set the six parameters for rescaling using function
  /// BackgroundRescalingYPhiUsingVectorForY(double v2, double v3, double v4, double psi, std::vector<double> values, std::vector<double> rap_binning);
  /// which is parametrized as
  ///  f(y,phi) = phi_term(phi) * rap_term(y)                                                                                                           
  /// where                                                                                                                                             
  ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))                                     
  ///  with four parameters v2, v3, v4, and psi.                                                                                                        
  /// rap_term(y) = bin content of the vector at position y       
  ///
  /// The user must set the parameters event-by-event. Example:
  vector<double> values;
  vector<double> binning;
  values.push_back(5);
  values.push_back(6);
  values.push_back(5.5);
  values.push_back(5);
  binning.push_back(-4);
  binning.push_back(-2);
  binning.push_back(0);
  binning.push_back(2);
  binning.push_back(4);
  contrib::BackgroundRescalingYPhiUsingVectorForY rescaling(0.1,0.1,0.001,0,values,binning);
  rescaling.use_rap_term(true);    // this is useful to check if the vectors with rapidity dependence have good sizes, if one wants to use also the rapidity rescaling.
  


  ///**** rescaling in heavy ion events using parametrized rapidity dependence: ****
  // Set the seven parameters for rescaling using function
  //  BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2)
  /// which is parametrized as
  ///  f(y,phi) = phi_term(phi) * rap_term(y)                                                                                                           
  /// where                                                                                                                                             
  ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))                                     
  ///  with four parameters v2, v3, v4, and psi.                                                                                                        
  /// rap_term(y) = a1*exp(-pow(y,2)/(2*sigma1^2)) + a2*exp(-pow(y,2)/(2*sigma2^2))         
  ///  with four parameters sigma1, sigma2, 1a, and a2.                                     
  ///
  /// You need to set the parameters event-by-event. Example:
  //contrib::BackgroundRescalingYPhi rescaling(0.1,0.1,0.001,0,1,5,0,10);
  //rescaling.use_rap_term(false);    // set to true, if you have derived also the rapidity terms for the rescaling



  ///**** rescaling using rapidity dependence stored in root TH1 object: ****
  // find the rapidity distribution of pileup particles from minimum bias events in a separate run. Fill a root TH1 histogram with this distribution.
  // Here as an example where the TH1 histogram is filled with a random distribution. Do not use this!!! Fill it with all particles (topoclusters) weighted by their pt using minimum bias events.
  /* TH1D* hist=new TH1D("hist","",100,-5,5);
  TF1 polynom("polynom","1+0.03*x*x-0.003*x*x*x*x",-5,5);
  for (int i=1;i<hist->GetNbinsX()+1;i++){
    hist->SetBinContent(i,polynom(hist->GetBinCenter(i)));
  }

  contrib::BackgroundRescalingYFromRoot<TH1D> rescaling(hist);
*/


  ///**** rescaling using rapidity and azimuth dependence stored in root TH2 object: ****
  // find the rapidity and azimuth distribution of pileup particles from minimum bias events in a separate run. Fill a root TH2 histogram with this distribution.
  // Here as an example where the TH2 histogram is filled with a random distribution. Do not use this!!! Fill it with all particles (topoclusters) weighted by their pt using minimum bias events.
  /*  TH2D* hist=new TH2D("hist","",100,-5,5,100,0,6.28319);
  TF1 polynom("polynom","1+0.03*x*x-0.003*x*x*x*x",-5,5);
  for (int i=1;i<hist->GetNbinsX()+1;i++){
    for (int j=1;j<hist->GetNbinsY()+1;j++){
      hist->SetBinContent(i,j,polynom(hist->GetXaxis()->GetBinCenter(i)));
    }
  }

  contrib::BackgroundRescalingYFromRoot<TH2D> rescaling(hist);
  */



  // This should be set up before event loop:                                                                                                                        
  double max_eta=4;   // specify the maximal pseudorapidity for the input particles. It is used for the subtraction. Particles with eta>|max_eta| are removed and not used during the subtraction (they are not returned). The same parameter should be used for the GridMedianBackgroundEstimator as it is demonstrated in this example. If JetMedianBackgroundEstimator is used, then lower parameter should be used  (to avoid including particles outside this range).                                     
  double max_eta_jet=3; // the maximal pseudorapidity for selected jets. Not important for the subtraction.
  // background estimation

  GridMedianBackgroundEstimator bge_rho(max_eta,0.5);  // maximal pseudo-rapidity cut is used inside ConstituentSubtraction, but in GridMedianBackgroundEstimator, the range is specified by maximal rapidity cut. Therefore, it is important to apply the same pseudo-rapidity cut also for particles used for background estimation (using function "set_particles") and also derive the rho dependence on rapidity using this max pseudo-rapidity cut to get the correct rescaling function!

  contrib::IterativeConstituentSubtractor subtractor;
  subtractor.set_distance_type(contrib::ConstituentSubtractor::deltaR); // free parameter for the type of distance between particle i and ghost k. There are two options: "deltaR" or "angle" which are defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean angle between the momenta
  vector<double> max_distances;
  max_distances.push_back(0.1);
  max_distances.push_back(0.2);
  vector<double> alphas;
  alphas.push_back(0);
  alphas.push_back(0);
  subtractor.set_parameters(max_distances,alphas); // in this example, 2 CS corrections will be performed: 1. correction with max_distance of 0.1 and alpha of 0, 2. correction with max_distance of 0.15 and alpha of 0. After the first correction, the scalar sum of pt from remaining ghosts is evaluated and redistributed uniformly accross the event.
  subtractor.set_ghost_removal(true);  // set to true if the ghosts (proxies) which were not used in the previous CS procedure should be removed for the next CS procedure
  subtractor.set_ghost_area(0.004); // parameter for the density of ghosts. The smaller, the better - but also the computation is slower.
  subtractor.set_max_eta(max_eta); // parameter for maximal |eta| cut. It is pspecified above.

  subtractor.set_background_estimator(&bge_rho);


  // For examples how to treat massive particles or how to use Selector, see example_event_wide.cc and example_iterative.cc


  subtractor.initialize();  // this is new compared to previous usages of ConstituentSubtractor! It should be used after specifying all the parameters and before event loop.





  // in event loop
  // read in input particles
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsEtaMax(max_eta)(hard_event);
  full_event = SelectorAbsEtaMax(max_eta)(full_event);

  cout << "# read an event with " << hard_event.size() << " signal particles and " << full_event.size() - hard_event.size() << " background particles with pseudo-rapidity |eta|<4" << endl;

 
  // jet definition and jet selector
  JetDefinition jet_def(antikt_algorithm, 0.7);
  Selector sel_jets = SelectorNHardest(3) * SelectorAbsEtaMax(max_eta_jet);

  // jet clustering
  ClusterSequence clust_seq_hard(hard_event, jet_def);
  ClusterSequence clust_seq_full(full_event, jet_def);
  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());


  // setting the rescaling:
  bge_rho.set_rescaling_class(&rescaling);
  bge_rho.set_particles(full_event);

  // print info (optional)
  cout << subtractor.description() << endl;


  // correction of the whole event with ConstituentSubtractor
  vector<PseudoJet> corrected_event=subtractor.subtract_event(full_event);

  ClusterSequence clust_seq_corr(corrected_event, jet_def);
  vector<PseudoJet> corrected_jets = sel_jets(clust_seq_corr.inclusive_jets());

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



