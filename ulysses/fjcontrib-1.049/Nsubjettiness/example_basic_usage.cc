//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-13
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  Run this example with:
//     ./example_basic_usage < ../data/single-event.dat
//
//  $Id: example_basic_usage.cc 841 2015-08-19 17:33:59Z jthaler $
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


#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <sstream>
#include "Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "Njettiness.hh"
#include "NjettinessPlugin.hh"
#include "XConePlugin.hh"

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how Nsubjettiness contrib works
  analyze(event);

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

// Helper Function for Printing out Jet Information
void PrintJets(const vector <PseudoJet>& jets, const TauComponents & components, bool showTotal = true);
void PrintXConeJets(const vector <PseudoJet>& jets, bool commentOut = false);
void PrintXConeAxes(const vector <PseudoJet>& jets, bool commentOut = false);

////////
//
//  Main Routine for Analysis 
//
///////

void analyze(const vector<PseudoJet> & input_particles) {
   
   ////////
   //
   //  Start of analysis.  First find anti-kT jets, then find N-subjettiness values of those jets
   //
   ///////
   
   // Initial clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm;
   double jet_rad = 1.00; // jet radius for anti-kt algorithm
   JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
   ClusterSequence clust_seq(input_particles,jetDef);
   vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());
   
   for (int j = 0; j < 2; j++) { // Two hardest jets per event

      // get the jet for analysis
      PseudoJet this_jet = antikt_jets[j];
      
      // only look at if harder than 200 GeV
      if (this_jet.perp() < 200.0) continue;
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "Analyzing Jet " << j + 1 << ":" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      
      ////////
      //
      //  Basic checks of tau values first
      //
      //  If you don't want to know the directions of the subjets,
      //  then you can use the simple function Nsubjettiness.
      //
      //  Recommended usage for Nsubjettiness:
      //  AxesMode:  KT_Axes(), WTA_KT_Axes(), OnePass_KT_Axes(), or OnePass_WTA_KT_Axes()
      //  MeasureMode:  Unnormalized_Measure(beta)
      //  beta with KT_Axes: 2.0
      //  beta with WTA_KT_Axes: anything greater than 0.0 (particularly good for 1.0)
      //  beta with OnePass_KT_Axes or OnePass_WTA_KT_Axes:  between 1.0 and 3.0
      //
      ///////
      
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << "N-subjettiness with Unnormalized Measure (in GeV)" << endl;
      cout << "beta = 1.0:  One-pass Winner-Take-All kT Axes" << endl;
      cout << "beta = 2.0:  One-pass E-Scheme kT Axes" << endl;
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      
      // Now loop through all options
      cout << setprecision(6) << right << fixed;
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      cout << setw(15) << "beta"
         << setw(14) << "tau1"
         << setw(14) << "tau2"
         << setw(14) << "tau3"
         << setw(14) << "tau2/tau1"
         << setw(14) << "tau3/tau2"
         << endl;
      
      
      // Define Nsubjettiness functions for beta = 1.0 using one-pass WTA KT axes
      double beta = 1.0;
      Nsubjettiness         nSub1_beta1(1,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      Nsubjettiness         nSub2_beta1(2,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      Nsubjettiness         nSub3_beta1(3,   OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      NsubjettinessRatio   nSub21_beta1(2,1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      NsubjettinessRatio   nSub32_beta1(3,2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
      
      // calculate Nsubjettiness values (beta = 1.0)
      double  tau1_beta1 =  nSub1_beta1(this_jet);
      double  tau2_beta1 =  nSub2_beta1(this_jet);
      double  tau3_beta1 =  nSub3_beta1(this_jet);
      double tau21_beta1 = nSub21_beta1(this_jet);
      double tau32_beta1 = nSub32_beta1(this_jet);
      
      // Output results (beta = 1.0)
      cout << setw(15) << 1.0
         << setw(14) << tau1_beta1
         << setw(14) << tau2_beta1
         << setw(14) << tau3_beta1
         << setw(14) << tau21_beta1
         << setw(14) << tau32_beta1
         << endl;
      
      // Repeat the above for beta = 2.0
      
      // Define Nsubjettiness functions for beta = 2.0, using one-pass ordinary KT axes
      beta = 2.0;
      Nsubjettiness         nSub1_beta2(1,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
      Nsubjettiness         nSub2_beta2(2,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
      Nsubjettiness         nSub3_beta2(3,   OnePass_KT_Axes(), UnnormalizedMeasure(beta));
      NsubjettinessRatio   nSub21_beta2(2,1, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
      NsubjettinessRatio   nSub32_beta2(3,2, OnePass_KT_Axes(), UnnormalizedMeasure(beta));
      
      // calculate Nsubjettiness values (beta = 2.0)
      double  tau1_beta2 =  nSub1_beta2(this_jet);
      double  tau2_beta2 =  nSub2_beta2(this_jet);
      double  tau3_beta2 =  nSub3_beta2(this_jet);
      double tau21_beta2 = nSub21_beta2(this_jet);
      double tau32_beta2 = nSub32_beta2(this_jet);

      // Output results (beta = 2.0)
      cout << setw(15) << 2.0
         << setw(14) << tau1_beta2
         << setw(14) << tau2_beta2
         << setw(14) << tau3_beta2
         << setw(14) << tau21_beta2
         << setw(14) << tau32_beta2
         << endl;

      
      // Additional information

      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Subjets found using beta = 1.0 tau values" << endl;
      PrintJets(nSub1_beta1.currentSubjets(),nSub1_beta1.currentTauComponents()); // these subjets have valid constituents
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub2_beta1.currentSubjets(),nSub2_beta1.currentTauComponents()); // these subjets have valid constituents
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub3_beta1.currentSubjets(),nSub3_beta1.currentTauComponents()); // these subjets have valid constituents
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Axes used for above beta = 1.0 tau values" << endl;
      
      PrintJets(nSub1_beta1.currentAxes(),nSub1_beta1.currentTauComponents(),false);
      //PrintJets(nSub1_beta1.seedAxes());  // For one-pass minimization, this would show starting seeds
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub2_beta1.currentAxes(),nSub2_beta1.currentTauComponents(),false);
      //PrintJets(nSub2_beta1.seedAxes());  // For one-pass minimization, this would show starting seeds
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub3_beta1.currentAxes(),nSub3_beta1.currentTauComponents(),false);
      //PrintJets(nSub3_beta1.seedAxes());  // For one-pass minimization, this would show starting seeds
      
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      
      
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Subjets found using beta = 2.0 tau values" << endl;
      PrintJets(nSub1_beta2.currentSubjets(),nSub1_beta2.currentTauComponents()); // these subjets have valid constituents
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub2_beta2.currentSubjets(),nSub2_beta2.currentTauComponents()); // these subjets have valid constituents
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub3_beta2.currentSubjets(),nSub3_beta2.currentTauComponents()); // these subjets have valid constituents
      cout << "-------------------------------------------------------------------------------------" << endl;
      
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Axes used for above beta = 2.0 tau values" << endl;
      
      PrintJets(nSub1_beta2.currentAxes(),nSub1_beta2.currentTauComponents(),false);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub2_beta2.currentAxes(),nSub2_beta2.currentTauComponents(),false);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(nSub3_beta2.currentAxes(),nSub3_beta2.currentTauComponents(),false);
      cout << "-------------------------------------------------------------------------------------" << endl;
      
   }


   ////////// The XCone Jet Algorithm ///////////////////////////

   ////////
   //
   //  We define a specific implementation of N-jettiness as a jet algorithm, which we call "XCone".
   //  This is the recommended version for all users.  
   //
   //  Recommended usage of XConePlugin is with beta = 2.0
   //  Beta = 1.0 is also useful as a recoil-free variant in the face of pile-up.
   //  
   ///////
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Using the XCone Jet Algorithm" << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;

   // Jet radius to use throughout
   double R = 0.5;

   // Define the jet finding plugins for beta = 1.0
   double beta = 1.0;
   
   // define the plugins
   XConePlugin xcone_plugin2_beta1(2, R, beta);
   XConePlugin xcone_plugin3_beta1(3, R, beta);
   XConePlugin xcone_plugin4_beta1(4, R, beta);
   
   // and the jet definitions
   JetDefinition xcone_jetDef2_beta1(&xcone_plugin2_beta1);
   JetDefinition xcone_jetDef3_beta1(&xcone_plugin3_beta1);
   JetDefinition xcone_jetDef4_beta1(&xcone_plugin4_beta1);
   
   // and the cluster sequences
   ClusterSequence xcone_seq2_beta1(input_particles, xcone_jetDef2_beta1);
   ClusterSequence xcone_seq3_beta1(input_particles, xcone_jetDef3_beta1);
   ClusterSequence xcone_seq4_beta1(input_particles, xcone_jetDef4_beta1);

   // and find the jets
   vector<PseudoJet> xcone_jets2_beta1 = xcone_seq2_beta1.inclusive_jets();
   vector<PseudoJet> xcone_jets3_beta1 = xcone_seq3_beta1.inclusive_jets();
   vector<PseudoJet> xcone_jets4_beta1 = xcone_seq4_beta1.inclusive_jets();

   // Do exactly the same for beta = 2.0 (which is the default, so no argument needed)
   
   // define the plugins
   XConePlugin xcone_plugin2_beta2(2, R);
   XConePlugin xcone_plugin3_beta2(3, R);
   XConePlugin xcone_plugin4_beta2(4, R);
   
   // and the jet definitions
   JetDefinition xcone_jetDef2_beta2(&xcone_plugin2_beta2);
   JetDefinition xcone_jetDef3_beta2(&xcone_plugin3_beta2);
   JetDefinition xcone_jetDef4_beta2(&xcone_plugin4_beta2);
   
   // and the cluster sequences
   ClusterSequence xcone_seq2_beta2(input_particles, xcone_jetDef2_beta2);
   ClusterSequence xcone_seq3_beta2(input_particles, xcone_jetDef3_beta2);
   ClusterSequence xcone_seq4_beta2(input_particles, xcone_jetDef4_beta2);

   // and find the jets
   vector<PseudoJet> xcone_jets2_beta2 = xcone_seq2_beta2.inclusive_jets();
   vector<PseudoJet> xcone_jets3_beta2 = xcone_seq3_beta2.inclusive_jets();
   vector<PseudoJet> xcone_jets4_beta2 = xcone_seq4_beta2.inclusive_jets();
   
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Using beta = " << setprecision(2) << 1.0 << ", R = " << setprecision(2) << R << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;

   PrintXConeJets(xcone_jets2_beta1);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeJets(xcone_jets3_beta1);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeJets(xcone_jets4_beta1);
   
   // The axes might point in a different direction than the jets
   // Using the njettiness_extras pointer (ClusterSequence::Extras) to access that information
   vector<PseudoJet> xcone_axes2_beta1 = njettiness_extras(xcone_seq2_beta1)->axes();
   vector<PseudoJet> xcone_axes3_beta1 = njettiness_extras(xcone_seq3_beta1)->axes();
   vector<PseudoJet> xcone_axes4_beta1 = njettiness_extras(xcone_seq4_beta1)->axes();
   
   cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
   cout << "Axes Used for Above Jets" << endl;
   
   PrintXConeAxes(xcone_axes2_beta1);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeAxes(xcone_axes3_beta1);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeAxes(xcone_axes4_beta1);

   
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Using beta = " << setprecision(2) << 2.0 << ", R = " << setprecision(2) << R << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;
   
   PrintXConeJets(xcone_jets2_beta2);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeJets(xcone_jets3_beta2);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeJets(xcone_jets4_beta2);
   
   // The axes might point in a different direction than the jets
   // Using the njettiness_extras pointer (ClusterSequence::Extras) to access that information
   vector<PseudoJet> xcone_axes2_beta2 = njettiness_extras(xcone_seq2_beta2)->axes();
   vector<PseudoJet> xcone_axes3_beta2 = njettiness_extras(xcone_seq3_beta2)->axes();
   vector<PseudoJet> xcone_axes4_beta2 = njettiness_extras(xcone_seq4_beta2)->axes();
   
   cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
   cout << "Axes Used for Above Jets" << endl;
   
   PrintXConeAxes(xcone_axes2_beta2);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeAxes(xcone_axes3_beta2);
   cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
   PrintXConeAxes(xcone_axes4_beta2);
   
   
   
   bool calculateArea = false;
   if (calculateArea) {
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Adding Area Information for beta = 2.0 (quite slow)" << endl;
      
      double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
      AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));
      
      // Defining cluster sequences with area
      ClusterSequenceArea xcone_seq_area2(input_particles, xcone_jetDef2_beta2, area_def);
      ClusterSequenceArea xcone_seq_area3(input_particles, xcone_jetDef3_beta2, area_def);
      ClusterSequenceArea xcone_seq_area4(input_particles, xcone_jetDef4_beta2, area_def);
      
      vector<PseudoJet> xcone_jets_area2 = xcone_seq_area2.inclusive_jets();
      vector<PseudoJet> xcone_jets_area3 = xcone_seq_area3.inclusive_jets();
      vector<PseudoJet> xcone_jets_area4 = xcone_seq_area4.inclusive_jets();
      
      PrintXConeJets(xcone_jets_area2);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintXConeJets(xcone_jets_area3);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintXConeJets(xcone_jets_area4);
   }
   
   cout << "-------------------------------------------------------------------------------------" << endl;
   cout << "Done Using the XCone Jet Algorithm" << endl;
   cout << "-------------------------------------------------------------------------------------" << endl;

}

void PrintJets(const vector <PseudoJet>& jets, const TauComponents & components, bool showTotal) {
   
   string commentStr = "";
   
   // gets extras information
   if (jets.size() == 0) return;
   
   // For printing out component tau information
   vector<double> subTaus = components.jet_pieces();
   double totalTau = components.tau();
   
   bool useArea = jets[0].has_area();
   
   // define nice tauN header
   int N = jets.size();
   stringstream ss(""); ss << "tau" << N; string tauName = ss.str();
   
   cout << fixed << right;
   
   cout << commentStr << setw(5) << "jet #" << "   "
   <<  setw(10) << "rap"
   <<  setw(10) << "phi"
   <<  setw(11) << "pt"
   <<  setw(11) << "m"
   <<  setw(11) << "e";
   if (jets[0].has_constituents()) cout <<  setw(11) << "constit";
   cout << setw(13) << tauName;
   if (useArea) cout << setw(10) << "area";
   cout << endl;
   
   
   // print out individual jet information
   for (unsigned i = 0; i < jets.size(); i++) {
      cout << commentStr << setw(5) << i+1  << "   "
      << setprecision(4) <<  setw(10) << jets[i].rap()
      << setprecision(4) <<  setw(10) << jets[i].phi()
      << setprecision(4) <<  setw(11) << jets[i].perp()
      << setprecision(4) <<  setw(11) << max(jets[i].m(),0.0) // needed to fix -0.0 issue on some compilers.
      << setprecision(4) <<  setw(11) << jets[i].e();
      if (jets[i].has_constituents()) cout << setprecision(4) <<  setw(11) << jets[i].constituents().size();
      cout << setprecision(6) <<  setw(13) << max(subTaus[i],0.0);
      if (useArea) cout << setprecision(4) << setw(10) << (jets[i].has_area() ? jets[i].area() : 0.0 );
      cout << endl;
   }
   
   // print out total jet
   if (showTotal) {
      fastjet::PseudoJet total = join(jets);
      
      cout << commentStr << setw(5) << "total" << "   "
      <<  setprecision(4) << setw(10) << total.rap()
      <<  setprecision(4) << setw(10) << total.phi()
      <<  setprecision(4) << setw(11) << total.perp()
      <<  setprecision(4) << setw(11) << max(total.m(),0.0) // needed to fix -0.0 issue on some compilers.
      <<  setprecision(4) <<  setw(11) << total.e();
      if (jets[0].has_constituents()) cout << setprecision(4)  <<  setw(11) << total.constituents().size();
      cout <<  setprecision(6) << setw(13) << totalTau;
      if (useArea) cout << setprecision(4) << setw(10) << (total.has_area() ? total.area() : 0.0);
      cout << endl;
   }
}


void PrintXConeJets(const vector <PseudoJet>& jets, bool commentOut) {
   
   string commentStr = "";
   if (commentOut) commentStr = "#";
   
   // gets extras information
   if (jets.size() == 0) return;
   const NjettinessExtras * extras = njettiness_extras(jets[0]);
   
   // bool useExtras = true;
   bool useExtras = (extras != NULL);
   bool useArea = jets[0].has_area();
   bool useConstit = jets[0].has_constituents();

   // define nice tauN header
   int N = jets.size();
   stringstream ss(""); ss << "tau" << N; string tauName = ss.str();
   
   cout << fixed << right;
   
   cout << commentStr << setw(5) << "jet #" << "   "
   <<  setw(10) << "rap"
   <<  setw(10) << "phi"
   <<  setw(11) << "pt"
   <<  setw(11) << "m"
   <<  setw(11) << "e";
   if (useConstit) cout <<  setw(11) << "constit";
   if (useExtras) cout << setw(14) << tauName;
   if (useArea) cout << setw(10) << "area";
   cout << endl;
   
   fastjet::PseudoJet total(0,0,0,0);
   int total_constit = 0;
   
   // print out individual jet information
   for (unsigned i = 0; i < jets.size(); i++) {
      cout << commentStr << setw(5) << i+1  << "   "
      << setprecision(4) <<  setw(10) << jets[i].rap()
      << setprecision(4) <<  setw(10) << jets[i].phi()
      << setprecision(4) <<  setw(11) << jets[i].perp()
      << setprecision(4) <<  setw(11) << max(jets[i].m(),0.0) // needed to fix -0.0 issue on some compilers.
      << setprecision(4) <<  setw(11) << jets[i].e();
      if (useConstit) cout << setprecision(4) <<  setw(11) << jets[i].constituents().size();
      if (useExtras) cout << setprecision(6) <<  setw(14) << max(extras->subTau(jets[i]),0.0);
      if (useArea) cout << setprecision(4) << setw(10) << (jets[i].has_area() ? jets[i].area() : 0.0 );
      cout << endl;
      total += jets[i];
      if (useConstit) total_constit += jets[i].constituents().size();
   }
   
   // print out total jet
   if (useExtras) {
      double beamTau = extras->beamTau();
      
      if (beamTau > 0.0) {
         cout << commentStr << setw(5) << " beam" << "   "
         <<  setw(10) << ""
         <<  setw(10) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(14) << setprecision(6) << beamTau
         << endl;
      }
      
      cout << commentStr << setw(5) << "total" << "   "
      <<  setprecision(4) << setw(10) << total.rap()
      <<  setprecision(4) << setw(10) << total.phi()
      <<  setprecision(4) << setw(11) << total.perp()
      <<  setprecision(4) << setw(11) << max(total.m(),0.0) // needed to fix -0.0 issue on some compilers.
      <<  setprecision(4) <<  setw(11) << total.e();
      if (useConstit) cout << setprecision(4) <<  setw(11) << total_constit;
      if (useExtras) cout <<  setprecision(6) << setw(14) << extras->totalTau();
      if (useArea) cout << setprecision(4) << setw(10) << (total.has_area() ? total.area() : 0.0);
      cout << endl;
   }
   
}


void PrintXConeAxes(const vector <PseudoJet>& jets, bool commentOut) {
   
   string commentStr = "";
   if (commentOut) commentStr = "#";
   
   // gets extras information
   if (jets.size() == 0) return;
   const NjettinessExtras * extras = njettiness_extras(jets[0]);
   
   // bool useExtras = true;
   bool useExtras = (extras != NULL);
   bool useArea = jets[0].has_area();

   // define nice tauN header
   int N = jets.size();
   stringstream ss(""); ss << "tau" << N; string tauName = ss.str();
   
   cout << fixed << right;
   
   cout << commentStr << setw(5) << "jet #" << "   "
   <<  setw(10) << "rap"
   <<  setw(10) << "phi"
   <<  setw(11) << "pt"
   <<  setw(11) << "m"
   <<  setw(11) << "e";
   if (useExtras) cout << setw(14) << tauName;
   if (useArea) cout << setw(10) << "area";
   cout << endl;
   
   fastjet::PseudoJet total(0,0,0,0);
   
   // print out individual jet information
   for (unsigned i = 0; i < jets.size(); i++) {
      cout << commentStr << setw(5) << i+1  << "   "
      << setprecision(4) <<  setw(10) << jets[i].rap()
      << setprecision(4) <<  setw(10) << jets[i].phi()
      << setprecision(4) <<  setw(11) << jets[i].perp()
      << setprecision(4) <<  setw(11) << max(jets[i].m(),0.0) // needed to fix -0.0 issue on some compilers.
      << setprecision(4) <<  setw(11) << jets[i].e();
      if (useExtras) cout << setprecision(6) <<  setw(14) << max(extras->subTau(jets[i]),0.0);
      if (useArea) cout << setprecision(4) << setw(10) << (jets[i].has_area() ? jets[i].area() : 0.0 );
      cout << endl;
      total += jets[i];
   }
   
   // print out total jet
   if (useExtras) {
      double beamTau = extras->beamTau();
      
      if (beamTau > 0.0) {
         cout << commentStr << setw(5) << " beam" << "   "
         <<  setw(10) << ""
         <<  setw(10) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(11) << ""
         <<  setw(14) << setprecision(6) << beamTau
         << endl;
      }
      
      cout << commentStr << setw(5) << "total" << "   "
      <<  setprecision(4) << setw(10) << total.rap()
      <<  setprecision(4) << setw(10) << total.phi()
      <<  setprecision(4) << setw(11) << total.perp()
      <<  setprecision(4) << setw(11) << max(total.m(),0.0) // needed to fix -0.0 issue on some compilers.
      <<  setprecision(4) <<  setw(11) << total.e()
      <<  setprecision(6) << setw(14) << extras->totalTau();
      if (useArea) cout << setprecision(4) << setw(10) << (total.has_area() ? total.area() : 0.0);
      cout << endl;
   }
   
}