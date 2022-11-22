//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-13
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
// Run this example with
//     ./example_advanced_usage < ../data/single-event.dat
//
//  $Id: example_advanced_usage.cc 933 2016-04-04 22:23:32Z jthaler $
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

// Simple class to store Axes along with a name for display
class AxesStruct {
   
private:
   // Shared Ptr so it handles memory management
   SharedPtr<AxesDefinition> _axes_def;
   
public:
   AxesStruct(const AxesDefinition & axes_def)
   : _axes_def(axes_def.create()) {}

   // Need special copy constructor to make it possible to put in a std::vector
   AxesStruct(const AxesStruct& myStruct)
   : _axes_def(myStruct._axes_def->create()) {}
   
   const AxesDefinition & def() const {return *_axes_def;}
   string description() const {return _axes_def->description();}
   string short_description() const {return _axes_def->short_description();}

};


// Simple class to store Measures to make it easier to put in std::vector
class MeasureStruct {
   
private:
   // Shared Ptr so it handles memory management
   SharedPtr<MeasureDefinition> _measure_def;
   
public:
   MeasureStruct(const MeasureDefinition& measure_def)
   : _measure_def(measure_def.create()) {}
   
   // Need special copy constructor to make it possible to put in a std::vector
   MeasureStruct(const MeasureStruct& myStruct)
   : _measure_def(myStruct._measure_def->create()) {}
   
   const MeasureDefinition & def() const {return *_measure_def;}
   string description() const {return _measure_def->description();}
   
};


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
void PrintJets(const vector <PseudoJet>& jets, bool commentOut = false);
void PrintAxes(const vector <PseudoJet>& jets, bool commentOut = false);
void PrintJetsWithComponents(const vector <PseudoJet>& jets, bool commentOut = false);

////////
//
//  Main Routine for Analysis 
//
///////

void analyze(const vector<PseudoJet> & input_particles) {
   
   ////////
   //
   //  This code will check multiple axes/measure modes
   //  First thing we do is establish the various modes we will check
   //
   ///////

   //Define characteristic test parameters to use here
   double p = 0.5;
   double delta = 10.0; // close to winner-take-all.  TODO:  Think about right value here.
   double R0 = 0.2;
   double Rcutoff = 0.5;
   double infinity = std::numeric_limits<int>::max();
   int nExtra = 2;
   int NPass = 10;

   // A list of all of the available axes modes
   vector<AxesStruct> _testAxes;
   _testAxes.push_back(KT_Axes());
   _testAxes.push_back(CA_Axes());
   _testAxes.push_back(AntiKT_Axes(R0));
   _testAxes.push_back(WTA_KT_Axes());
   _testAxes.push_back(WTA_CA_Axes());
   _testAxes.push_back(GenKT_Axes(p, R0));
   _testAxes.push_back(WTA_GenKT_Axes(p, R0));
   _testAxes.push_back(GenET_GenKT_Axes(delta, p, R0));

   _testAxes.push_back(OnePass_KT_Axes());
   _testAxes.push_back(OnePass_AntiKT_Axes(R0));
   _testAxes.push_back(OnePass_WTA_KT_Axes());
   _testAxes.push_back(OnePass_GenKT_Axes(p, R0));
   _testAxes.push_back(OnePass_WTA_GenKT_Axes(p, R0));
   _testAxes.push_back(OnePass_GenET_GenKT_Axes(delta, p, R0));

   _testAxes.push_back(Comb_GenKT_Axes(nExtra, p, R0));
   _testAxes.push_back(Comb_WTA_GenKT_Axes(nExtra, p, R0));
   _testAxes.push_back(Comb_GenET_GenKT_Axes(nExtra, delta, p, R0));

   // manual axes (should be identical to kt axes)
   _testAxes.push_back(Manual_Axes());
   _testAxes.push_back(OnePass_Manual_Axes());

   // these axes are not checked during make check since they do not give reliable results
   _testAxes.push_back(OnePass_CA_Axes()); // not recommended 
   _testAxes.push_back(OnePass_WTA_CA_Axes()); // not recommended
   _testAxes.push_back(MultiPass_Axes(NPass));
   _testAxes.push_back(MultiPass_Manual_Axes(NPass));
   int num_unchecked = 4; // number of unchecked axes

   //
   // Note:  Njettiness::min_axes is not guarenteed to give a global
   // minimum, only a local minimum, and different choices of the random
   // number seed can give different results.  For that reason,
   // the one-pass minimization are recommended over min_axes.
   //
   
   // Getting a smaller list of recommended axes modes
   // These are the ones that are more likely to give sensible results (and are all IRC safe)
   vector<AxesStruct> _testRecommendedAxes;
   _testRecommendedAxes.push_back(KT_Axes());
   _testRecommendedAxes.push_back(WTA_KT_Axes());
   _testRecommendedAxes.push_back(OnePass_KT_Axes());
   _testRecommendedAxes.push_back(OnePass_WTA_KT_Axes());

   // new axes options added in most recent version of Nsubjettiness
   // these are separate from above since they should only be defined with a cutoff value for sensible results
   vector<AxesStruct> _testAlgorithmRecommendedAxes;
   _testAlgorithmRecommendedAxes.push_back(GenET_GenKT_Axes(1.0, 1.0, Rcutoff));
   _testAlgorithmRecommendedAxes.push_back(GenET_GenKT_Axes(infinity, 1.0, Rcutoff));
   _testAlgorithmRecommendedAxes.push_back(GenET_GenKT_Axes(1.0, 0.5, Rcutoff));
   _testAlgorithmRecommendedAxes.push_back(OnePass_GenET_GenKT_Axes(1.0, 1.0, Rcutoff));
   _testAlgorithmRecommendedAxes.push_back(OnePass_GenET_GenKT_Axes(infinity, 1.0, Rcutoff));
   _testAlgorithmRecommendedAxes.push_back(OnePass_GenET_GenKT_Axes(1.0, 0.5, Rcutoff));

   // Getting some of the measure modes to test
   // (When applied to a single jet we won't test the cutoff measures,
   // since cutoffs aren't typically helpful when applied to single jets)
   // Note that we are calling measures by their MeasureDefinition
   vector<MeasureStruct> _testMeasures;
   _testMeasures.push_back(  NormalizedMeasure(1.0, 1.0, pt_R));
   _testMeasures.push_back(UnnormalizedMeasure(1.0     , pt_R));
   _testMeasures.push_back(  NormalizedMeasure(2.0, 1.0, pt_R));
   _testMeasures.push_back(UnnormalizedMeasure(2.0     , pt_R));

   // When doing Njettiness as a jet algorithm, want to test the cutoff measures.
   // (Since they are not senisible without a cutoff)
   vector<MeasureStruct> _testCutoffMeasures;
   _testCutoffMeasures.push_back(UnnormalizedCutoffMeasure(1.0, Rcutoff, pt_R));
   _testCutoffMeasures.push_back(UnnormalizedCutoffMeasure(2.0, Rcutoff, pt_R));
   // new measures added in the most recent version of NSubjettiness
   _testCutoffMeasures.push_back(ConicalMeasure(1.0, Rcutoff));
   _testCutoffMeasures.push_back(ConicalMeasure(2.0, Rcutoff));
   _testCutoffMeasures.push_back(OriginalGeometricMeasure(Rcutoff));
   _testCutoffMeasures.push_back(ModifiedGeometricMeasure(Rcutoff));
   _testCutoffMeasures.push_back(ConicalGeometricMeasure(1.0, 1.0, Rcutoff));
   _testCutoffMeasures.push_back(ConicalGeometricMeasure(2.0, 1.0, Rcutoff));
   _testCutoffMeasures.push_back(XConeMeasure(1.0, Rcutoff)); // Should be identical to ConicalGeometric
   _testCutoffMeasures.push_back(XConeMeasure(2.0, Rcutoff));

   /////// N-subjettiness /////////////////////////////

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
   
   // clust_seq.delete_self_when_unused();
   // small number to show equivalence of doubles
   double epsilon = 0.0001;
   
   for (int j = 0; j < 2; j++) { // Two hardest jets per event
      if (antikt_jets[j].perp() < 200) continue;
            
      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Analyzing Jet " << j + 1 << ":" << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;

      
      ////////
      //
      //  Basic checks of tau values first
      //
      //  If you don't want to know the directions of the subjets,
      //  then you can use the simple function Nsubjettiness.
      //
      //  Recommended usage for Nsubjettiness:
      //  AxesMode:  kt_axes, wta_kt_axes, onepass_kt_axes, or onepass_wta_kt_axes
      //  MeasureMode:  unnormalized_measure
      //  beta with kt_axes: 2.0
      //  beta with wta_kt_axes: anything greater than 0.0 (particularly good for 1.0)
      //  beta with onepass_kt_axes or onepass_wta_kt_axes:  between 1.0 and 3.0
      //
      ///////
      
      
      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Outputting N-subjettiness Values" << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;

      
      // Now loop through all options
      cout << setprecision(6) << right << fixed;
      for (unsigned iM = 0; iM < _testMeasures.size(); iM++) {
         
         cout << "-----------------------------------------------------------------------------------------------" << endl;
         cout << _testMeasures[iM].description() << ":" << endl;
         cout << setw(25) << "AxisMode"
            << setw(14) << "tau1"
            << setw(14) << "tau2"
            << setw(14) << "tau3"
            << setw(14) << "tau2/tau1"
            << setw(14) << "tau3/tau2"
            << endl;
         
         for (unsigned iA = 0; iA < _testAxes.size(); iA++) {

            // Current axes/measure modes and particles
            const PseudoJet         & my_jet      = antikt_jets[j];
            const vector<PseudoJet>   particles   = my_jet.constituents();
            const AxesDefinition    & axes_def    = _testAxes[iA].def();
            const MeasureDefinition & measure_def = _testMeasures[iM].def();
            
            // This case doesn't work, so skip it.
            // if (axes_def.givesRandomizedResults()) continue;

            // define Nsubjettiness functions
            Nsubjettiness        nSub1(1,    axes_def, measure_def);
            Nsubjettiness        nSub2(2,    axes_def, measure_def);
            Nsubjettiness        nSub3(3,    axes_def, measure_def);

            // define manual axes when they are necessary (should be identical to KT_Axes)
            if (axes_def.needsManualAxes()) {
               JetDefinition manual_jetDef(fastjet::kt_algorithm,
                                             fastjet::JetDefinition::max_allowable_R, 
                                             // new WinnerTakeAllRecombiner(), 
                                             fastjet::E_scheme, 
                                             fastjet::Best);

               fastjet::ClusterSequence manual_clustSeq(particles, manual_jetDef);

               nSub1.setAxes(manual_clustSeq.exclusive_jets(1));
               nSub2.setAxes(manual_clustSeq.exclusive_jets(2));
               nSub3.setAxes(manual_clustSeq.exclusive_jets(3));
            }

            // calculate Nsubjettiness values
            double tau1 = nSub1(my_jet);
            double tau2 = nSub2(my_jet);
            double tau3 = nSub3(my_jet);
                        
            //These should only happen if the axes are not manual and are not multipass
            double tau21, tau32;
            if (!axes_def.needsManualAxes() && !axes_def.givesRandomizedResults()) {
               // An entirely equivalent, but painful way to calculate is:
               double tau1alt = measure_def(particles,axes_def(1,particles,&measure_def));
               double tau2alt = measure_def(particles,axes_def(2,particles,&measure_def));
               double tau3alt = measure_def(particles,axes_def(3,particles,&measure_def));
               assert(tau1alt == tau1);
               assert(tau2alt == tau2);
               assert(tau3alt == tau3);

               NsubjettinessRatio   nSub21(2,1, axes_def, measure_def);
               NsubjettinessRatio   nSub32(3,2, axes_def, measure_def);
               tau21 = nSub21(my_jet);
               tau32 = nSub32(my_jet);

               // Make sure calculations are consistent
               if (!_testAxes[iA].def().givesRandomizedResults()) {
                  assert(abs(tau21 - tau2/tau1) < epsilon);
                  assert(abs(tau32 - tau3/tau2) < epsilon);
               }
            }
            else {
               tau21 = tau2/tau1;
               tau32 = tau3/tau2;
            }

            string axesName = _testAxes[iA].short_description();
            string left_hashtag;

            // comment out with # because MultiPass uses random number seed, or because axes do not give reliable results (those at the end of axes vector)
            if (_testAxes[iA].def().givesRandomizedResults() || iA >= (_testAxes.size() - num_unchecked)) left_hashtag = "#";
            else left_hashtag = " ";

            // Output results:
            cout << std::right
               << left_hashtag
               << setw(23)
               << axesName
               << ":"
               << setw(14) << tau1
               << setw(14) << tau2
               << setw(14) << tau3
               << setw(14) << tau21
               << setw(14) << tau32
               << endl;
         }
      }

      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Done Outputting N-subjettiness Values" << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;

      
      ////////
      //
      //  Finding axes/jets found by N-subjettiness partitioning
      //
      //  This uses the component_results function to get the subjet information
      //
      ///////

      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Outputting N-subjettiness Subjets" << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;

      
      // Loop through all options, this time setting up jet finding
      cout << setprecision(6) << left << fixed;
      for (unsigned iM = 0; iM < _testMeasures.size(); iM++) {
         
         for (unsigned iA = 0; iA < _testRecommendedAxes.size(); iA++) {

            const PseudoJet         & my_jet      = antikt_jets[j];
            const AxesDefinition    & axes_def    = _testRecommendedAxes[iA].def();
            const MeasureDefinition & measure_def = _testMeasures[iM].def();
            
            // This case doesn't work, so skip it.
            if (axes_def.givesRandomizedResults()) continue;
            
            // define Nsubjettiness functions
            Nsubjettiness        nSub1(1,    axes_def, measure_def);
            Nsubjettiness        nSub2(2,    axes_def, measure_def);
            Nsubjettiness        nSub3(3,    axes_def, measure_def);
            
            // get component results
            TauComponents tau1comp = nSub1.component_result(my_jet);
            TauComponents tau2comp = nSub2.component_result(my_jet);
            TauComponents tau3comp = nSub3.component_result(my_jet);
            
            vector<PseudoJet> jets1 = tau1comp.jets();
            vector<PseudoJet> jets2 = tau2comp.jets();
            vector<PseudoJet> jets3 = tau3comp.jets();

            vector<PseudoJet> axes1 = tau1comp.axes();
            vector<PseudoJet> axes2 = tau2comp.axes();
            vector<PseudoJet> axes3 = tau3comp.axes();

            cout << "-----------------------------------------------------------------------------------------------" << endl;
            cout << measure_def.description() << ":" << endl;
            cout << axes_def.description() << ":" << endl;
            
            bool commentOut = false;
            if (axes_def.givesRandomizedResults()) commentOut = true;  // have to comment out min_axes, because it has random values
            
            // This helper function tries to find out if the jets have tau information for printing
            PrintJetsWithComponents(jets1,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJetsWithComponents(jets2,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJetsWithComponents(jets3,commentOut);
            
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
            cout << "Axes Used for Above Subjets" << endl;

            PrintAxes(axes1,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintAxes(axes2,commentOut);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintAxes(axes3,commentOut);
            
         }
      }
      
      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Done Outputting N-subjettiness Subjets" << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;

   }
   

   ////////// the XCone Jet Algorithm ///////////////////////////

   ////////
   //
   //  We define a specific implementation of N-jettiness as a jet algorithm, which we call "XCone".
   //  This is the recommended version for all users.  
   //
   //  Recommended usage of XConePlugin is with beta = 2.0
   //  Beta = 1.0 is also useful as a recoil-free variant in the face of pile-up.
   //
   ///////
   
   cout << "-----------------------------------------------------------------------------------------------" << endl;
   cout << "Using the XCone Jet Algorithm" << endl;
   cout << "-----------------------------------------------------------------------------------------------" << endl;

   //create list of various values of beta
   vector<double> betalist;
   betalist.push_back(1.0);
   betalist.push_back(2.0);
   unsigned int n_betas = betalist.size();   

   for (unsigned iB = 0; iB < n_betas; iB++) {
   
      double beta = betalist[iB];

      // define the plugins
      XConePlugin xcone_plugin2(2, Rcutoff, beta);
      XConePlugin xcone_plugin3(3, Rcutoff, beta);
      XConePlugin xcone_plugin4(4, Rcutoff, beta);

      // and the jet definitions
      JetDefinition xcone_jetDef2(&xcone_plugin2);
      JetDefinition xcone_jetDef3(&xcone_plugin3);
      JetDefinition xcone_jetDef4(&xcone_plugin4);

      // and the cluster sequences
      ClusterSequence xcone_seq2(input_particles, xcone_jetDef2);
      ClusterSequence xcone_seq3(input_particles, xcone_jetDef3);
      ClusterSequence xcone_seq4(input_particles, xcone_jetDef4);

      // and associated extras for more information
      const NjettinessExtras * extras2 = njettiness_extras(xcone_seq2);
      const NjettinessExtras * extras3 = njettiness_extras(xcone_seq3);
      const NjettinessExtras * extras4 = njettiness_extras(xcone_seq4);

      // and find the jets
      vector<PseudoJet> xcone_jets2 = xcone_seq2.inclusive_jets();
      vector<PseudoJet> xcone_jets3 = xcone_seq3.inclusive_jets();
      vector<PseudoJet> xcone_jets4 = xcone_seq4.inclusive_jets();

      // (alternative way to find the jets)
      //vector<PseudoJet> xcone_jets2 = extras2->jets();
      //vector<PseudoJet> xcone_jets3 = extras3->jets();
      //vector<PseudoJet> xcone_jets4 = extras4->jets();

      cout << "-----------------------------------------------------------------------------------------------" << endl;
      cout << "Using beta = " << setprecision(2) << beta << ", Rcut = " << setprecision(2) << Rcutoff << endl;
      cout << "-----------------------------------------------------------------------------------------------" << endl;
      
      PrintJets(xcone_jets2);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(xcone_jets3);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintJets(xcone_jets4);

      // The axes might point in a different direction than the jets
      // Using the NjettinessExtras pointer (ClusterSequence::Extras) to access that information         
      vector<PseudoJet> xcone_axes2 = extras2->axes();
      vector<PseudoJet> xcone_axes3 = extras3->axes();
      vector<PseudoJet> xcone_axes4 = extras4->axes();
      
      cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
      cout << "Axes Used for Above Jets" << endl;
      
      PrintAxes(xcone_axes2);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintAxes(xcone_axes3);
      cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      PrintAxes(xcone_axes4);
      
      bool calculateArea = false;
      if (calculateArea) {
         cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
         cout << "Adding Area Information (quite slow)" << endl;
         
         double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
         AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));
         
         // Defining cluster sequences with area
         ClusterSequenceArea xcone_seq_area2(input_particles, xcone_jetDef2, area_def);
         ClusterSequenceArea xcone_seq_area3(input_particles, xcone_jetDef3, area_def);
         ClusterSequenceArea xcone_seq_area4(input_particles, xcone_jetDef4, area_def);
         
         vector<PseudoJet> xcone_jets_area2 = xcone_seq_area2.inclusive_jets();
         vector<PseudoJet> xcone_jets_area3 = xcone_seq_area3.inclusive_jets();
         vector<PseudoJet> xcone_jets_area4 = xcone_seq_area4.inclusive_jets();
         
         PrintJets(xcone_jets_area2);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(xcone_jets_area3);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(xcone_jets_area4);
      }
   }
   
   cout << "-----------------------------------------------------------------------------------------------" << endl;
   cout << "Done Using the XCone Jet Algorithm" << endl;
   cout << "-----------------------------------------------------------------------------------------------" << endl;



   ////////// N-jettiness as a jet algorithm ///////////////////////////

   ////////
   //
   //  The user can also defined N-jettiness as a jet algorithm more generally, using different choice
   //  for measures and for axis finding. 
   //
   //  Recommended usage of NjettinessPlugin (event-wide)
   //  AxesMode:  wta_kt_axes or onepass_wta_kt_axes
   //  MeasureMode:  unnormalized_measure
   //  beta with wta_kt_axes: anything greater than 0.0 (particularly good for 1.0)
   //  beta with onepass_wta_kt_axes:  between 1.0 and 3.0
   //
   //  Note that the user should find that the usage of Conical Geometric Measure beta = 1.0 with 
   //  GenET_GenKT_Axes(std::numeric_limits<int>::max(), 1.0, Rcutoff) should be identical to XCone beta = 1.0,
   //  and Conical Geometric Measure beta = 2.0 with GenET_GenKT_Axes(1.0, 0.5, Rcutoff) should be identical to
   //  XCone beta = 2.0.
   //
   ///////
   
   cout << "-----------------------------------------------------------------------------------------------" << endl;
   cout << "Using N-jettiness as a Jet Algorithm" << endl;
   cout << "-----------------------------------------------------------------------------------------------" << endl;

   
   for (unsigned iM = 0; iM < _testCutoffMeasures.size(); iM++) {
      
      for (unsigned iA = 0; iA < _testAlgorithmRecommendedAxes.size(); iA++) {
         
         const AxesDefinition    & axes_def    = _testAlgorithmRecommendedAxes[iA].def();
         const MeasureDefinition & measure_def = _testCutoffMeasures[iM].def();
         
         // define the plugins
         NjettinessPlugin njet_plugin2(2, axes_def,measure_def);
         NjettinessPlugin njet_plugin3(3, axes_def,measure_def);
         NjettinessPlugin njet_plugin4(4, axes_def,measure_def);
   
         // and the jet definitions
         JetDefinition njet_jetDef2(&njet_plugin2);
         JetDefinition njet_jetDef3(&njet_plugin3);
         JetDefinition njet_jetDef4(&njet_plugin4);

         // and the cluster sequences
         ClusterSequence njet_seq2(input_particles, njet_jetDef2);
         ClusterSequence njet_seq3(input_particles, njet_jetDef3);
         ClusterSequence njet_seq4(input_particles, njet_jetDef4);

         // and associated extras for more information
         const NjettinessExtras * extras2 = njettiness_extras(njet_seq2);
         const NjettinessExtras * extras3 = njettiness_extras(njet_seq3);
         const NjettinessExtras * extras4 = njettiness_extras(njet_seq4);

         // and find the jets
         vector<PseudoJet> njet_jets2 = njet_seq2.inclusive_jets();
         vector<PseudoJet> njet_jets3 = njet_seq3.inclusive_jets();
         vector<PseudoJet> njet_jets4 = njet_seq4.inclusive_jets();

         // (alternative way to find the jets)
         //vector<PseudoJet> njet_jets2 = extras2->jets();
         //vector<PseudoJet> njet_jets3 = extras3->jets();
         //vector<PseudoJet> njet_jets4 = extras4->jets();

         cout << "-----------------------------------------------------------------------------------------------" << endl;
         cout << measure_def.description() << ":" << endl;
         cout << axes_def.description() << ":" << endl;
         
         PrintJets(njet_jets2);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_jets3);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintJets(njet_jets4);

         // The axes might point in a different direction than the jets
         // Using the NjettinessExtras pointer (ClusterSequence::Extras) to access that information         
         vector<PseudoJet> njet_axes2 = extras2->axes();
         vector<PseudoJet> njet_axes3 = extras3->axes();
         vector<PseudoJet> njet_axes4 = extras4->axes();
         
         cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
         cout << "Axes Used for Above Jets" << endl;
         
         PrintAxes(njet_axes2);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintAxes(njet_axes3);
         cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
         PrintAxes(njet_axes4);
         
         bool calculateArea = false;
         if (calculateArea) {
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
            cout << "Adding Area Information (quite slow)" << endl;
            
            double ghost_maxrap = 5.0; // e.g. if particles go up to y=5
            AreaDefinition area_def(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxrap));
            
            // Defining cluster sequences with area
            ClusterSequenceArea njet_seq_area2(input_particles, njet_jetDef2, area_def);
            ClusterSequenceArea njet_seq_area3(input_particles, njet_jetDef3, area_def);
            ClusterSequenceArea njet_seq_area4(input_particles, njet_jetDef4, area_def);

            vector<PseudoJet> njet_jets_area2 = njet_seq_area2.inclusive_jets();
            vector<PseudoJet> njet_jets_area3 = njet_seq_area3.inclusive_jets();
            vector<PseudoJet> njet_jets_area4 = njet_seq_area4.inclusive_jets();

            PrintJets(njet_jets_area2);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(njet_jets_area3);
            cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
            PrintJets(njet_jets_area4);
         }
         
      }
   }
   
   cout << "-----------------------------------------------------------------------------------------------" << endl;
   cout << "Done Using N-jettiness as a Jet Algorithm" << endl;
   cout << "-----------------------------------------------------------------------------------------------" << endl;
   
   
   
   // Below are timing tests for the developers
   double do_timing_test = false;
   if (do_timing_test) {
      
      clock_t clock_begin, clock_end;
      double num_iter;
      
      cout << setprecision(6);
      
      num_iter = 1000;
      
      double R0 = 0.5;
      double beta = 2.0;
      double N = 6;
      
      
      // AKT
      JetDefinition aktDef = JetDefinition(antikt_algorithm,R0,E_scheme,Best);
      
      // XC
      XConePlugin xconePlugin(N, R0, beta);
      JetDefinition xconeDef = JetDefinition(&xconePlugin);
      
      // pXC
      PseudoXConePlugin pseudoxconePlugin(N, R0, beta);
      JetDefinition pseudoxconeDef = JetDefinition(&pseudoxconePlugin);
      
      //AKT
      cout << "Timing for " << aktDef.description() << endl;
      clock_begin = clock();
      for (int t = 0; t < num_iter; t++) {
         ClusterSequence clust_seq(input_particles,aktDef);
         clust_seq.inclusive_jets();
      }
      clock_end = clock();
      cout << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per AKT"<< endl;
      
      // XC
      cout << "Timing for " << xconeDef.description() << endl;
      clock_begin = clock();
      for (int t = 0; t < num_iter; t++) {
         ClusterSequence clust_seq(input_particles,xconeDef);
         clust_seq.inclusive_jets();
      }
      clock_end = clock();
      cout << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per XCone"<< endl;
      
      // pXC
      cout << "Timing for " << pseudoxconePlugin.description() << endl;
      clock_begin = clock();
      for (int t = 0; t < num_iter; t++) {
         ClusterSequence clust_seq(input_particles,pseudoxconeDef);
         clust_seq.inclusive_jets();
      }
      clock_end = clock();
      cout << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per PseudoXCone"<< endl;
      
      
   }
}


void PrintJets(const vector <PseudoJet>& jets, bool commentOut) {
   
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


void PrintAxes(const vector <PseudoJet>& jets, bool commentOut) {
   
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

void PrintJetsWithComponents(const vector <PseudoJet>& jets, bool commentOut) {
   
   string commentStr = "";
   if (commentOut) commentStr = "#";
   
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
   cout << setw(14) << tauName;
   if (useArea) cout << setw(10) << "area";
   cout << endl;
   
   fastjet::PseudoJet total(0,0,0,0);
   double total_tau = 0;
   int total_constit = 0;
   
   
   // print out individual jet information
   for (unsigned i = 0; i < jets.size(); i++) {
      double thisTau = jets[i].structure_of<TauComponents>().tau();
      
      cout << commentStr << setw(5) << i+1  << "   "
      << setprecision(4) <<  setw(10) << jets[i].rap()
      << setprecision(4) <<  setw(10) << jets[i].phi()
      << setprecision(4) <<  setw(11) << jets[i].perp()
      << setprecision(4) <<  setw(11) << max(jets[i].m(),0.0) // needed to fix -0.0 issue on some compilers.
      << setprecision(4) <<  setw(11) << jets[i].e();
      if (jets[i].has_constituents()) cout << setprecision(4) <<  setw(11) << jets[i].constituents().size();
      cout << setprecision(6) <<  setw(14) << max(thisTau,0.0);
      if (useArea) cout << setprecision(4) << setw(10) << (jets[i].has_area() ? jets[i].area() : 0.0 );
      cout << endl;
      total += jets[i];
      total_tau += thisTau;
      if (jets[i].has_constituents()) total_constit += jets[i].constituents().size();
   }
   
   cout << commentStr << setw(5) << "total" << "   "
   <<  setprecision(4) << setw(10) << total.rap()
   <<  setprecision(4) << setw(10) << total.phi()
   <<  setprecision(4) << setw(11) << total.perp()
   <<  setprecision(4) << setw(11) << max(total.m(),0.0) // needed to fix -0.0 issue on some compilers.
   <<  setprecision(4) <<  setw(11) << total.e();
   if (jets[0].has_constituents()) cout << setprecision(4)  <<  setw(11) << total_constit;
   cout <<  setprecision(6) << setw(14) << total_tau;
   if (useArea) cout << setprecision(4) << setw(10) << (total.has_area() ? total.area() : 0.0);
   cout << endl;
   
}



