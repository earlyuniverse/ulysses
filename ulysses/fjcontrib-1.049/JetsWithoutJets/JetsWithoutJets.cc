// JetsWithoutJets Package
// Questions/Comments? danbert@mit.edu jthaler@mit.edu
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: JetsWithoutJets.cc 554 2014-02-21 19:02:08Z danbert $
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

#include "JetsWithoutJets.hh"

using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace jwj {

//////
//
// Shape Trimming
//
//////


//----------------------------------------------------------------------
// Helper for SelectorShapeTrimming. 
// Class derived from SelectorWorker: pass, terminator, applies_jet_by_jet, and description are overloaded. 
// It can't be applied to an individual jet since it requires information about the neighborhood of the jet.

class SW_ShapeTrimming: public SelectorWorker {
   
private:
   
   double _Rjet, _ptcut, _Rsub, _fcut;
   bool _useLocalStorage;
   
public:
   
   SW_ShapeTrimming(double Rjet, double ptcut, double Rsub, double fcut, bool useLocalStorage = true) : _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _useLocalStorage(useLocalStorage) {};
   
   virtual bool pass(const PseudoJet &) const {
      if (!applies_jet_by_jet())	
         throw Error("Cannot apply this selector worker to an individual jet");
      return false;
   }
   
   virtual void terminator(vector<const PseudoJet *> & jets) const {
      
      // Copy pointers content to a vector of PseudoJets. Check if the pointer has been already nullified. 
      vector<PseudoJet> my_jets;
      vector<unsigned int> indices; 
      for (unsigned int i=0; i < jets.size(); i++){
         if (jets[i]){
            indices.push_back(i);
            my_jets.push_back(*jets[i]);
         } 	
      }
      
      EventStorage myEventStorage(_Rjet,_ptcut,_Rsub,_fcut,_useLocalStorage,false);
      myEventStorage.establishStorage(my_jets);
      
      
      for (unsigned int i=0; i < myEventStorage.size(); i++){
         
         if (!myEventStorage[i].includeParticle()) {
            jets[indices[i]] = NULL;
            continue;
         }
      }
   }
   
   virtual bool applies_jet_by_jet() const {return false;}
   
   std::string jetParameterString() const {
      std::stringstream stream;
      stream << "R_jet=" << _Rjet << ", pT_cut=" << _ptcut << ", R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   
   virtual string description() const {
      return "Shape trimmer, " + jetParameterString();
   }
}; 

//----------------------------------------------------------------------
// Helper for SelectorJetShapeTrimming. 
// Class derived from SelectorWorker: pass, terminator, applies_jet_by_jet, and description are overloaded. 
// This cannot be applied to an individual jet.
class SW_JetShapeTrimming: public SelectorWorker {
   
private:
   
   double _Rsub, _fcut;
   
public:
   
   SW_JetShapeTrimming(double Rsub, double fcut) : _Rsub(Rsub), _fcut(fcut) {};
   
   virtual bool pass(const PseudoJet &) const {
      if (!applies_jet_by_jet())	
         throw Error("Cannot apply this selector worker to an individual jet");
      return false;
   }
   
   virtual void terminator(vector<const PseudoJet *> & jets) const {
      
      // Copy pointers content to a vector of PseudoJets. Check if the pointer has been already nullified. 
      vector<PseudoJet> my_jets;
      vector<unsigned int> indices; 
      for (unsigned int i=0; i < jets.size(); i++){
         if (jets[i]){
            indices.push_back(i);
            my_jets.push_back(*jets[i]);
         } 	
      }
      
      FunctionScalarPtSum sumPt = FunctionScalarPtSum();
      
      // take sum pt of given jet
      double pt_Rjet = sumPt(my_jets);
      EventStorage myJetStorage(_Rsub,pt_Rjet*_fcut,false,false);
      myJetStorage.establishStorage(my_jets);
      
      for (unsigned int i=0; i < myJetStorage.size(); i++){
         
         if (!myJetStorage[i].includeParticle()) {
            jets[indices[i]] = NULL;
            continue;
         }
         
      }
   }
   
   virtual bool applies_jet_by_jet() const {return false;}
   
   std::string jetParameterString() const {
      std::stringstream stream;
      stream << "R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   
   virtual string description() const {
      return "Jet shape trimmer, " + jetParameterString();
   }
}; 


Selector SelectorShapeTrimming(double Rjet, double ptcut, double Rsub, double fcut){
   return Selector(new SW_ShapeTrimming(Rjet,ptcut,Rsub,fcut));
} 

Selector SelectorJetShapeTrimming(double Rsub, double fcut){
   return Selector(new SW_JetShapeTrimming(Rsub,fcut));
} 


//////
//
// Multiple pT_cut values
//
//////

// helper to sort a vector of vector<double> by comparing just the first entry.
bool _mySortFunction (std::vector<double> v_0, std::vector<double> v_1) { return (v_0[0] > v_1[0]); }


// As described in the appendix of the physics paper, one can construct the inverse 
// of an event shape with respect to pT by calculating all of the possible pTi,R values and sorting them.
// This helper function accomplishes that task.

void JetLikeEventShape_MultiplePtCutValues::_storeLocalInfo(const std::vector<PseudoJet> particles) {
   
   EventStorage myEventStorage(_Rjet,0.0,_Rsub,_fcut,_useLocalStorage);
   myEventStorage.establishStorage(particles);
   
   _functionArray.resize(0);
   
   
   for (unsigned int i = 0; i < myEventStorage.size(); i++){ 
      std::vector<double> point(2);
      point[0] = myEventStorage[i].pt_in_Rjet();
      point[1] = myEventStorage[i].weight()*_measurement->result(myEventStorage.particles_near_to(i));
      _functionArray.push_back(point);
   }
}

void JetLikeEventShape_MultiplePtCutValues::_buildStepFunction() {
   
   std::sort (_functionArray.begin(), _functionArray.end(), _mySortFunction);
   
   // Calculating the cumulative sum of the measurements up to a given pt value  
   if (!_functionArray.empty()){
      for (unsigned int i = 1; i < _functionArray.size(); i++) _functionArray[i][1] += _functionArray[i-1][1];
      
   }
}


// helper to compare in a vector of vector<double> the first entry with an external val
bool _myCompFunction_0 (std::vector<double> v, double val) { return (v[0] < val); }	

// get the event shape for a given ptcut
double JetLikeEventShape_MultiplePtCutValues::eventShapeFor(const double ptcut_0) const {
   
   double eventShape = 0.0;  
   
   if(ptcut_0 <= _functionArray.front()[0]) eventShape = (*lower_bound(_functionArray.rbegin(),_functionArray.rend(),ptcut_0,_myCompFunction_0))[1];
   
   return (eventShape);
}

// helper to compare in a vector of vector<double> the second entry with an external val
bool _myCompFunction_1 (std::vector<double> v, double val) { return (v[1] < val); }

// get the ptcut for a given event shape
double JetLikeEventShape_MultiplePtCutValues::ptCutFor(const double eventShape_0) const {
   
   double ptcut = 0.0;  
   
   double new_eventShape_0 = eventShape_0 - _offset;
   
   if ( new_eventShape_0 <= 0 || new_eventShape_0 > _functionArray.back()[1] ) {
      throw Error("Event shape value not valid");
   } 
   
   else ptcut =  (*lower_bound(_functionArray.begin(),_functionArray.end(),new_eventShape_0,_myCompFunction_1))[0];
   
   return(ptcut);
}

//////
//
// Njet with multiple R_jet
//
//////

// Similar to the same named function in JetLikeEventShape_MultiplePtCutValues, except now we have
// to store a lot more information.

void ShapeJetMultiplicity_MultipleRValues::_buildStepFunction(const std::vector<PseudoJet> particles) {
   
   EventStorage myEventStorage(_Rsub,0.0,false,false);
   myEventStorage.establishStorage(particles);
   
   
   unsigned int N = myEventStorage.size();
   
   std::vector< std::vector<double> > myValues;
   myValues.resize(0);
   
   for(unsigned int i = 0; i < N; i++) {
      unsigned int j = i+1;
      while(j < N) {
         vector<double> myPair(3);
         myPair[0] = sqrt(myEventStorage[i].deltaRsq(myEventStorage[j]));
         myPair[1] = i;
         myPair[2] = j;
         myValues.push_back(myPair);
         j++;
      }
   }
   
   std::sort (myValues.begin(), myValues.end(), _mySortFunction);
   
   _functionArray.clear();
   _functionArray.resize(myValues.size());
   
   // Initialize
   double _tot_pt=0.0;
   for(unsigned int i=0; i<N; i++) _tot_pt += myEventStorage[i].pt();
   
   std::vector<double> pTR(N);
   std::fill(pTR.begin(), pTR.end(), _tot_pt);
   
   if (_trim) {
      std::vector<double> pTRsub;
      for (unsigned int i = 0; i < N; i++) pTRsub.push_back(myEventStorage[i].pt_in_Rsub());
      
      for (unsigned int i = 0; i < myValues.size(); i++) {
         
         _functionArray[i].push_back(myValues[i][0]);
         
         int id1 = (int)myValues[i][1];
         int id2 = (int)myValues[i][2];
         
         pTR[id1] -= myEventStorage[id2].pt();
         pTR[id2] -= myEventStorage[id1].pt();
         
         double myEventShape = 0;
         
         for(unsigned int j = 0; j < N; j++) {
            if(pTR[j] >= _ptcut && pTRsub[j] / pTR[j] >= _fcut) myEventShape += myEventStorage[j].pt() / pTR[j];
         }
         _functionArray[i].push_back(myEventShape);
      }
   }
   
   else {
      
      for (unsigned int i = 0; i < myValues.size(); i++) {
         
         _functionArray[i].push_back(myValues[i][0]);
         
         unsigned int id1 = (int)myValues[i][1];
         unsigned int id2 = (int)myValues[i][2];
         
         pTR[id1] -= myEventStorage[id2].pt();
         pTR[id2] -= myEventStorage[id1].pt();
         
         double myEventShape = 0;
         
         for(unsigned int j = 0; j < N; j++) {
            if(pTR[j] >= _ptcut) myEventShape += myEventStorage[j].pt() / pTR[j];
         }
         _functionArray[i].push_back(myEventShape);
      }
   }
   
   
   myValues.clear(); 
}		

// returns event shape for a given Rjet value
double ShapeJetMultiplicity_MultipleRValues::eventShapeFor(const double Rjet_0) const {
   
   double eventShape = 0.0;  
   
   if( Rjet_0 < _Rsub ) throw Error("Rjet < Rsub");
   else if (Rjet_0 < 0) throw Error("Negative Rjet");
   
   else if( Rjet_0 > _functionArray.front()[0] ) eventShape = _functionArray.front()[1];
   else eventShape = (*lower_bound(_functionArray.rbegin(),_functionArray.rend(),Rjet_0,_myCompFunction_0))[1];
   
   return(eventShape);
}



//////
//
// Finding Jet Axes with the Event Shape Density
//
//////

// Find the local axes in a region of size R around each particle
void EventShapeDensity_JetAxes::_find_local_axes(const std::vector<PseudoJet> & particles) {
   
   _myParticles = particles;
   _N = particles.size();
   
   _axes.resize(0);
   _Njet_weights.resize(0);
   _pt_weights.resize(0);
   
   // Indexing particles (note that this is happening internally, so shouldn't modify user input)
   for (unsigned int i=0; i<_N; i++) _myParticles[i].set_user_index(i);
   
   EventStorage myEventStorage(_Rjet,_ptcut,_Rjet,1.0,_useLocalStorage);
   myEventStorage.establishStorage(_myParticles);	
   
   
   // Find axis associated with each particle.	
   for (unsigned int i=0; i<_N; i++) {
      //Recombiner has to carry user_index information.
      FunctionJetAxis axis(_jetDef); 
      
      vector<PseudoJet> near_particles=myEventStorage.particles_near_to(i);
      
      int myAxis = -1;
      double pt_w = 0;	
      double Njet_w = 0;	
      
      if(myEventStorage[i].includeParticle()) {
         myAxis = axis(near_particles).user_index();
         pt_w = myEventStorage[i].pt(); 	
         Njet_w = myEventStorage[i].weight();
      } 
      
      _axes.push_back(myAxis);
      _pt_weights.push_back(pt_w);
      _Njet_weights.push_back(Njet_w);
   }
}

// Take the axes calculated by _find_local_axes and figure out what the final
// global axes will be.  With the global consistency condition, this decreases
// the number of axes under consideration. 
void EventShapeDensity_JetAxes::find_axes_and_weights() {
   
   
   //If requested establish global consistency.
   //Loop over axes and reassign unstable axes until no unstable axis is left.
   if(_applyGlobalConsistency) {
      int n_unstable_axes;
      do {
         n_unstable_axes = 0;	
         for (unsigned int i=0; i<_N; i++) {
            if(_axes[i]!= -1 && !_isStable(_axes[i])){ 
               n_unstable_axes++;
               _axes[i] = _axes[_axes[i]];
            }	
         }
      } while(n_unstable_axes>0);
   }
   
   //Find distinct axes. 
   //Output a vector of Pseudojets (sorted by pT), with pT of each axis corresponding to the pT weight.
   vector<double> tot_Njet_weights(_N,0), tot_pt_weights(_N,0);
   
   for (unsigned int i=0; i<_N; i++) {
      if (_axes[i] != -1) {
         tot_pt_weights[_axes[i]] += _pt_weights[i];
         tot_Njet_weights[_axes[i]] += _Njet_weights[i];
      }
   }
   
   _distinctAxes.resize(0);
   _tot_Njet_weights.resize(0);
   
   for (unsigned int i=0; i<_N; i++){
      if(tot_pt_weights[i] > 0) {
         PseudoJet myAxis = _lightLikeVersion(_myParticles[i]) * tot_pt_weights[i];
         _distinctAxes.push_back(myAxis);
      }
   }
   
   _distinctAxes = sorted_by_pt(_distinctAxes);
   
   //Find the corresponding N_jet weights.
   for(unsigned int i=0; i<_distinctAxes.size(); i++) _tot_Njet_weights.push_back(tot_Njet_weights[_distinctAxes[i].user_index()]);		
}

// Check stability: a particle defines a stable wta axis if it is its own wta or if it is pointing to a particle that does not pass pTcut.
bool EventShapeDensity_JetAxes::_isStable(const int thisAxis) const {
   
   bool answer = false;
   if(_axes[thisAxis] == thisAxis || _axes[thisAxis] == -1) answer = true;
   return(answer);
}



} // namespace jwj

FASTJET_END_NAMESPACE
