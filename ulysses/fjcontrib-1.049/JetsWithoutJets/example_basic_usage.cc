// Example showing basic usage of JetsWithoutJets classes.
//
// Compile it with "make example_basic_usage" and run it with
//
//   ./example_basic_usage < ../data/single-event.dat
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: example_basic_usage.cc 554 2014-02-21 19:02:08Z danbert $
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
   
   //////////
   // Basic event shape analysis
   //////////
   
   // Classes derived from JetLikeEventShape take as input a vector<PseudoJet>
   // representing the event, and return a number.
   // The simplest constructor requires you to specify Rjet and pTcut.

   
   ShapeJetMultiplicity Nj(Rjet,pTcut);// ShapeJetMultiplicity = jet counting
  
   ShapeScalarPt HT(Rjet,pTcut);// ShapeScalarPt = HT (summed jet pT)
  
   ShapeMissingPt HTmiss(Rjet,pTcut);// ShapeMissingPt = missing pT
  
   cout << setprecision(6);
   cout << "#########" << endl;
   cout << "## Example of basic event shape analysis" << endl;
   cout << "#########" << endl;
   cout << "Jet parameters: R_jet=" << Rjet << ", pTcut=" << pTcut <<endl;
   cout << "#####" << endl;
   cout << "N_jet=" << Nj(input_particles) << endl;
   cout << "H_T=" << HT(input_particles) << endl;
   cout << "Missing H_T=" << HTmiss(input_particles) << endl;
   
   //////////
   // Testing different trimming methods
   //////////
   
   // Different trimming mechanisms
   Selector eventShapeTrimmer=SelectorShapeTrimming(Rjet,pTcut,Rsub,fcut); // Event shape trimming
   JetShapeTrimmer jetShapeTrimmer(Rsub,fcut);// Jet shape trimming
   Filter treeTrimmer(Rsub, SelectorPtFractionMin(fcut));// Standard trimming
   
   // Clustering with anti-kt algorithm
   JetAlgorithm algorithm = antikt_algorithm; 
   JetDefinition jetDef = JetDefinition(algorithm,Rjet);
   ClusterSequence clust_seq(input_particles,jetDef);
   PseudoJet hardestJet  = sorted_by_pt(clust_seq.inclusive_jets(pTcut))[0];
   
   // Clustering after applying event shape trimming
   ClusterSequence clust_seq_estrim(eventShapeTrimmer(input_particles),jetDef);
   PseudoJet hardestJet_estrim  = sorted_by_pt(clust_seq_estrim.inclusive_jets(pTcut))[0];
   
   double untrimmedMass = hardestJet.m();
   double treeTrimmedMass = treeTrimmer(hardestJet).m();
   double jetShapeTrimmedMass = jetShapeTrimmer(hardestJet).m();
   double eventShapeTrimmedMass = hardestJet_estrim.m();
   
   cout << setprecision(6);
   cout << "#########" << endl;
   cout << "## Example of different trimming methods" << endl;
   cout << "#########" << endl;
   cout << "Anti-kT jet parameters: R_jet=" << Rjet << ", pTcut=" << pTcut <<endl;
   cout << "Trimming parameters: R_sub=" << Rsub << ", fcut=" << fcut <<endl;
   cout << "#########" << endl;
   cout << "Untrimmed Mass = " <<  untrimmedMass << endl;
   cout << "Tree Trimmed Mass = " <<  treeTrimmedMass << endl;
   cout << "Jet Shape Trimmed Mass = " <<  jetShapeTrimmedMass << endl;
   cout << "Event Shape Trimmed Mass = " <<  eventShapeTrimmedMass << endl;


   //////////
   // Event shapes with built-in trimming
   //////////

   // One can calculate the trimmed version of an event shape.
   // In this case the constructor requires you to specify Rjet, pTcut, Rsub and fcut.
   
   // N.B.: trimming an event and then calculating the event shape is NOT
   // the same as calculating the trimmed event shape.  See the README for
   // more information.

   ShapeJetMultiplicity Nj_trim(Rjet,pTcut,Rsub,fcut);  // jet counting with trimming
 
   ShapeScalarPt HT_trim(Rjet,pTcut,Rsub,fcut); // HT with trimming
 
   ShapeMissingPt HTmiss_trim(Rjet,pTcut,Rsub,fcut); // missing pT with trimming

   cout << setprecision(6);
   cout << "#########" << endl;
   cout << "## Example of trimmed event shapes" << endl;
   cout << "#########" << endl;
   cout << "Jet parameters: R_jet=" << Rjet << ", pTcut=" << pTcut <<endl;
   cout << "Trimming parameters: R_sub=" << Rsub << ", fcut=" << fcut <<endl;
   cout << "#########" << endl;
   cout << "N_jet_trim (using built-in)=" << Nj_trim(input_particles) << endl;
   cout << "N_jet_trim (using trimmer)=" << Nj(eventShapeTrimmer(input_particles)) << endl; // N.B.: Not exactly equivalent
   cout << "H_T_trim (using built-in)=" << HT_trim(input_particles) << endl;
   cout << "H_T_trim (using trimmer)=" << HT(eventShapeTrimmer(input_particles)) << endl; // N.B.: Not exactly equivalent
   cout << "Missing H_T_trim (using built-in)=" << HTmiss_trim(input_particles) << endl;
   cout << "Missing H_T_trim (using trimmer)=" << HTmiss(eventShapeTrimmer(input_particles)) << endl; // N.B.: Not exactly equivalent
   cout << "#########" << endl;


}
