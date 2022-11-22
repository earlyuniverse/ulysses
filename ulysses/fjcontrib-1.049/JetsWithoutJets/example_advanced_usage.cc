// Example showing advanced usage of JetsWithoutJets classes.
// Please see example_basic_usage for usage of basic classes first.
//
// Compile it with "make example_advanced_usage" and run it with
//
//   ./example_advanced_usage < ../data/single-event.dat
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: example_advanced_usage.cc 554 2014-02-21 19:02:08Z danbert $
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
#include <iostream>
#include <sstream>
#include <ctime>


#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/tools/Filter.hh"

#include "JetsWithoutJets.hh" // In external code, this should be fastjet/contrib/JetsWithoutJets.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::jwj;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){
   
   //----------------------------------------------------------
   // read in input particles
   vector<PseudoJet> event;
   read_event(event);
   cout << "#########"  << endl;
   cout << "## Read an event with " << event.size() << " particles" << endl;
   
   //----------------------------------------------------------
   // illustrate how this JetsWithoutJets contrib works
   
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

// Main code
void analyze(const vector<PseudoJet> & input_particles) {
   
   // Jet parameters. 
   double Rjet = 1.0;
   double pTcut = 200.0;
   
   // Subjet parameters. Need for trimming.
   double Rsub = 0.2;
   double fcut = 0.05;
   double ptsubcut = 50.0;  // additional ptsubcut for subjet multiplicity
   
   //////////
   // Multiple pT analysis.
   //////////
   
   // These use classes derived from JetLikeEventShape_MultiplePtCutValues, where 
   // the full pT information is kept in calculating the event shape.
   // One can calculate the value of the event shape at a given pT cut, or
   // one can calculate the inverse function to find the value of the pT cut
   // needed to return a given value of the event shape.
   
   ShapeJetMultiplicity_MultiplePtCutValues Nj_allPt(Rjet);
   ShapeJetMultiplicity_MultiplePtCutValues Nj_allPt_trim(Rjet,Rsub,fcut); // with trimming
   Nj_allPt.set_input(input_particles);
   Nj_allPt_trim.set_input(input_particles);
   
   cout << "#########" << endl;
   cout << "## Example of multiple pT analysis" << endl;
   cout << "#########" << endl;
   cout << "Jet parameters: R_jet=" << Rjet << ", pTcut=" << pTcut << " (unless otherwise stated)" <<endl;
   cout << "Trimming parameters: R_sub=" << Rsub << ", fcut=" << fcut <<endl;
   cout << "#####" << endl;
   cout << "N_jet from full pT_cut function=" <<  Nj_allPt.eventShapeFor(pTcut) << endl;
   cout << "N_jet_trim from full pT_cut function=" << Nj_allPt_trim.eventShapeFor(pTcut) << endl;
   cout << "#####" << endl;
   cout << "Multiple pT_cut" << endl;
   double pTcut_values[] = {10,20,50,100,150}; 
   for (unsigned int i=0; i<5; i++){
   	cout << "pT_cut=" << pTcut_values[i] << ", N_jet=" << Nj_allPt.eventShapeFor(pTcut_values[i]) << endl;
   }
   cout << "#####" << endl;
   cout << "pT of the hardest jet=" <<  Nj_allPt.ptCutFor(1) << endl;
   cout << "pT of the hardest trimmed jet=" << Nj_allPt_trim.ptCutFor(1) << endl;
  
   
   // Can get the full function array, though this is not output on the example.ref file
   // This is a vector of vectors of double. Each vector contains ptcut/Nj pair
   std::vector< std::vector<double> > functionArray_pt =  Nj_allPt.functionArray();
   
   //////////
   // Multiple R analysis.
   //////////
   
   // These use the class ShapeJetMultiplicity_MultipleRValues (at present there is no
   // general class for multiple R observables because of the high computational time). 
   // Here, the full R information is kept.
   
   ShapeJetMultiplicity_MultipleRValues Nj_allR(pTcut);
   ShapeJetMultiplicity_MultipleRValues Nj_allR_trim(pTcut,Rsub,fcut); // with trimming
   Nj_allR.set_input(input_particles);
   Nj_allR_trim.set_input(input_particles);
   
   cout << "#########" << endl;
   cout << "## Example of multiple R analysis" << endl;
   cout << "#########" << endl;
   cout << "Jet parameters: R_jet=" << Rjet << " (unless otherwise stated), pTcut=" << pTcut <<endl;
   cout << "Trimming parameters: R_sub=" << Rsub << ", fcut=" << fcut <<endl;
   cout << "#####" << endl;
   cout << "N_jet from full R function=" <<  Nj_allR.eventShapeFor(Rjet) << endl;
   cout << "N_jet_trim from full R function=" <<  Nj_allR_trim.eventShapeFor(Rjet) << endl;
   cout << "#####" << endl;
   cout << "Multiple Rjet" << endl;
   double Rjet_values[] = {0.2,0.4,0.6,0.8,1}; 
   for(unsigned int i=0; i<5; i++){
   	cout << "R=" << Rjet_values[i] << ", N_jet=" << Nj_allR.eventShapeFor(Rjet_values[i]) << endl;
   }
   
   // Can get the full function array, though this is not output on the example.ref file
   // This is a vector of vectors of double. Each vector contains Rjet/Nj pair
   std::vector< std::vector<double> > functionArray_R =  Nj_allR.functionArray();
   
   
   //////////
   // Finding Jet Axes
   //////////
   
   EventShapeDensity_JetAxes axes_eventShape(Rjet,pTcut); //Default local recombination metric is antikt with winner-take-all scheme
   axes_eventShape.set_input(input_particles);  //Loads particles and finds axes and weights
   
   //Get jet axes (by default sorted by pt). Pt of each axis corresponds to its pt weight. 	
   vector<PseudoJet> axes = axes_eventShape.axes();
   //Get N_jet weights
   vector<double> Njw = axes_eventShape.Njet_weights();
   
   //If you want turn on global consistency and recalculate axes and weights
   axes_eventShape.setGlobalConsistencyCheck(true);
   axes_eventShape.find_axes_and_weights(); // uses stored information to save time; don't need to set_input again
   vector<PseudoJet> axes_gc = axes_eventShape.axes();
   vector<double> Njw_gc = axes_eventShape.Njet_weights();
   
   
   //If you want to use a different metric for local recombination, Cambridge-Aachen for example
   EventShapeDensity_JetAxes axes_eventShape_CA(Rjet,pTcut,cambridge_algorithm);
   axes_eventShape_CA.set_input(input_particles);
   
   //Get jet axes and Njet weights	
   vector<PseudoJet> axes_CA = axes_eventShape_CA.axes();
   vector<double> Njw_CA = axes_eventShape_CA.Njet_weights();
   
   
   //Output hardest jet eta/phi/pT/N_jet weight
   cout << "#########" << endl;
   cout << "## Hardest jet axis and weights" << endl;
   cout << "#########" << endl;
   cout << "Without global consistency and anti-kt metric:" << endl;
   cout << "Rap_hardest_jet=" <<  axes[0].rap() <<endl;
   cout << "Phi_hardest_jet=" <<  axes[0].phi() <<endl;
   cout << "pt_weight_hardest_jet=" <<  axes[0].pt() <<endl;
   cout << "N_jet_weight_hardest_jet=" <<  Njw[0] <<endl;
   cout << "#" << endl;
   cout << "With global consistency and anti-kt metric:" << endl;
   cout << "Rap_hardest_jet=" <<  axes_gc[0].rap() <<endl;
   cout << "Phi_hardest_jet=" <<  axes_gc[0].phi() <<endl;
   cout << "pt_weight_hardest_jet=" <<  axes_gc[0].pt() <<endl;
   cout << "N_jet_weight_hardest_jet=" <<  Njw_gc[0] <<endl;
   cout << "#" << endl;
   cout << "Without global consistency and cambridge metric:" << endl;
   cout << "Rap_hardest_jet=" <<  axes_CA[0].rap() <<endl;
   cout << "Phi_hardest_jet=" <<  axes_CA[0].phi() <<endl;
   cout << "pt_weight_hardest_jet=" <<  axes_CA[0].pt() <<endl;
   cout << "N_jet_weight_hardest_jet=" <<  Njw_CA[0] <<endl;
   
   //////////
   // Showing a wide range of event shapes
   //////////
   
   // This emphasizes that the jet-like event shapes are all derived from
   // a common base class
   vector<JetLikeEventShape*> myShapes;
   vector<string > myNames;
   
   myShapes.push_back( new ShapeJetMultiplicity(Rjet,pTcut) );
   myNames.push_back( "Njet");
   myShapes.push_back( new ShapeJetMultiplicity(Rjet,pTcut,Rsub,fcut) ); //with trimming
   myNames.push_back( "Njet_trim");
   
   myShapes.push_back( new ShapeScalarPt(Rjet,pTcut) );
   myNames.push_back( "HT");
   myShapes.push_back( new ShapeScalarPt(Rjet,pTcut,Rsub,fcut) ); //with trimming
   myNames.push_back( "HT_trim");
   
   myShapes.push_back( new ShapeMissingPt(Rjet,pTcut) );
   myNames.push_back( "MissHT");
   myShapes.push_back( new ShapeMissingPt(Rjet,pTcut,Rsub,fcut) ); //with trimming
   myNames.push_back( "MissHT_trim");
   
   myShapes.push_back( new ShapeSummedMass(Rjet,pTcut) );
   myNames.push_back( "SigmaM");
   myShapes.push_back( new ShapeSummedMass(Rjet,pTcut,Rsub,fcut) ); //with trimming
   myNames.push_back( "SigmaM_trim");
   
   myShapes.push_back( new ShapeSummedMassSquared(Rjet,pTcut) );
   myNames.push_back( "SigmaM2");
   myShapes.push_back( new ShapeSummedMassSquared(Rjet,pTcut,Rsub,fcut) ); //with trimming
   myNames.push_back( "SigmaM2_trim");
   
   myShapes.push_back( new ShapeTrimmedSubjetMultiplicity(Rjet,pTcut,Rsub,fcut,ptsubcut) ); //with trimming
   myNames.push_back( "SigmaNsub_trim");
   
   for (int n = -2; n <= 2; n++) {
      myShapes.push_back( new ShapeScalarPtToN(n,Rjet,pTcut) );
      myShapes.push_back( new ShapeScalarPtToN(n,Rjet,pTcut,Rsub,fcut) ); //with trimming
      
      stringstream myStream;
      myStream << "GenHT_" << n;
      myNames.push_back( myStream.str());
      myNames.push_back( myStream.str() + "_trim");
   }
   
   // outputing the variables
   cout << "#########" << endl;
   cout << "## Examples showing more general shapes" << endl;
   cout << "#########" << endl;
   cout << "Jet parameters: R_jet=" << Rjet << ", pTcut=" << pTcut << endl;
   cout << "Trimming parameters: R_sub=" << Rsub << ", fcut=" << fcut << endl;
   cout << "#####" << endl;
   for (unsigned int i = 0; i < myShapes.size(); i++) {
      cout << myShapes.at(i)->description() << endl; // ideally each even shape should have a description
      
      myShapes.at(i)->setUseLocalStorage(false);
      double valueNoStorage = myShapes.at(i)->result(input_particles);
      myShapes.at(i)->setUseLocalStorage(true);
      double valueYesStorage = myShapes.at(i)->result(input_particles);
      
      cout << myNames.at(i) << "="
      << valueNoStorage << endl;
      cout << "(Local Storage works? " << (valueNoStorage == valueYesStorage ? "Yes" : "No!!") << ")" << endl;
      cout << "#" << endl;
      
      delete myShapes.at(i);
   }
   
   
   //////////
   // User defined JetLikeEventShape and using EventStorage
   //////////
   
   // Need to write a measurement function. For example, use FunctionScalarPtSum (already pre-defined) to get the ShapeScalarPt
   class ShapeScalarPt_new: public JetLikeEventShape {
	public:
      ShapeScalarPt_new(double Rjet, double ptcut): JetLikeEventShape(new FunctionScalarPtSum(),Rjet,ptcut) {}
   	ShapeScalarPt_new(double Rjet, double ptcut, double Rsub, double fcut) : JetLikeEventShape(new FunctionScalarPtSum(),Rjet,ptcut,Rsub,fcut) {}
   	~ShapeScalarPt_new(){}
   };
	
   ShapeScalarPt_new HT_new(Rjet,pTcut);

   // Built-in HT
   ShapeScalarPt HT(Rjet,pTcut);
   
   cout << "#########" << endl;
   cout << "## Example of JetLikeEventShape defined with an external measurement" << endl;
   cout << "#########" << endl;
   cout << "HT_new = " << HT_new(input_particles) << endl;
   cout << "compare to built-in HT = "<< HT(input_particles)<<endl;
   
   // Also, if you want to do multiple measurements with the same parameters ( Rjet,pTcut and Rsub,fcut if trimming is on) can build 
   // one event storage for all and reduce computation time
   
   // Input of EventStorage is Rjet,pTcut,bool for use of LocalStorage (true by default), bool for storing neighbours(true by default)
   // or Rjet,pTcut,Rsub,fcut,bool for use of LocalStorage (true by default), bool for storing neighbours(true by default).
   // Use the former if don't need to do trimming.
   // In this example we don't do trimming.
   EventStorage myStorage(Rjet,pTcut);
   myStorage.establishStorage(input_particles);
   
   vector<MyFunctionOfVectorOfPseudoJets<double>*> myMeasurements;
   vector<string> myNames1;
   myMeasurements.push_back( new FunctionUnity() );
   myNames1.push_back("Njet");
   myMeasurements.push_back( new FunctionScalarPtSum() );
   myNames1.push_back("HT");
   myMeasurements.push_back( new FunctionScalarPtSumToN(2) );
   myNames1.push_back("HT^2");
   myMeasurements.push_back( new FunctionInvariantMass() );
   myNames1.push_back("SigmaM");
   myMeasurements.push_back( new FunctionInvariantMassSquared() );
   myNames1.push_back("SigmaM2");
   
   cout << "#########" << endl;
   cout << "## Example of multiple measurements using a single storage" << endl;
   cout << "#########" << endl;
   
   for(unsigned int i=0; i<myMeasurements.size(); i++){
      JetLikeEventShape myShape(myMeasurements[i],Rjet,pTcut);
      cout<<"#"<<endl;
      cout<<myShape.description()<<endl;
      cout<<myNames1[i]<<"="<<myShape.result(myStorage)<<endl;
   }
   cout << "#########" << endl;	
   
   // timing tests for the developers
   double do_timing_test = false;
   if (do_timing_test) {
      clock_t clock_begin, clock_end;
      double num_iter;
      
      cout << setprecision(6);
      
      num_iter = 1000;
      clock_begin = clock();
      ShapeScalarPt NjetA(Rjet,pTcut);
      cout << "timing " << NjetA.description() << endl;
      for (int t = 0; t < num_iter; t++) {
         NjetA(input_particles);
      }
      clock_end = clock();
      cout << "Method A: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per Njet"<< endl;
      
      num_iter = 1000;
      clock_begin = clock();
      JetLikeEventShape NjetB(new FunctionScalarPtSum(),Rjet,pTcut);
      NjetB.setUseLocalStorage(false);
      cout << "timing " << NjetB.description() << endl;
      for (int t = 0; t < num_iter; t++) {
         NjetB(input_particles);
      }
      clock_end = clock();
      cout << "Method B: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per Njet"<< endl;
      
      
   }
   
}
