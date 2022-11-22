// JetsWithoutJets Package
// Questions/Comments? danbert@mit.edu jthaler@mit.edu
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: JetsWithoutJets.hh 554 2014-02-21 19:02:08Z danbert $
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

#ifndef __FASTJET_CONTRIB_JETSWITHOUTJETS_HH__
#define __FASTJET_CONTRIB_JETSWITHOUTJETS_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/tools/Transformer.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SharedPtr.hh"
#include <string>
#include <algorithm> 

#include "EventStorage.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace jwj {

//////
//
// Defining a function of a vector of PseudoJets
//
//////

//----------------------------------------------------------------------
// MyFunctionOfVectorOfPseudoJets
// In current version of FastJet there is no standard interface for a function
// that takes a vector of PseudoJets as argument. This class serves as a base class
// and it is the analog of FunctionOfPseudoJet, it only differs in taking a vector
// of PseudoJets as argument. Should a standard interface become available, this
// class could be removed and the classes below (like JetLikeEventShape,
// FunctionUnity, etc.) could derive from the standard base class.

template<typename TOut>
class MyFunctionOfVectorOfPseudoJets {
   
public:
   
   MyFunctionOfVectorOfPseudoJets(){}
   
   MyFunctionOfVectorOfPseudoJets(const TOut &constant_value);
   
   virtual ~MyFunctionOfVectorOfPseudoJets(){}
   
   virtual std::string description() const{ return "";}
   
   virtual TOut result(const std::vector<PseudoJet> &pjs) const = 0;
   
   TOut operator()(const std::vector<PseudoJet> &pjs) const { return result(pjs);}
};


//////
//
// Extendable Jet-like Event Shapes
// (works on any function of vector<PseudoJet> that returns double)
//
//////


//----------------------------------------------------------------------
// JetLikeEventShape
// Event shape that corresponds to a jet-based observable (summed over all jets).
// It takes a Function (derived from MyFunctionOfVectorOfPseudoJets) that would 
// correspond to a measurement done on each jet as an argument and it applies it to cones or radius R around each particle, 
// then sums over all such "jets" which are above pTcut threshold. 
// Optional trimming built in.  Note that in general trimming first,
// then applying the event shape gives a different answer.
// Only works with quantities that are doubles for each jet
// (could templatize, but easier to deal with case by case)
// If one wants to derive one's own JetLikeEventShape that is not based on a measurement,
// then one has to overload the virtual result(EventStorage & myStorage) function.

class JetLikeEventShape: public MyFunctionOfVectorOfPseudoJets<double>{
   
public:
   
   JetLikeEventShape(): _measurement(NULL),_useLocalStorage(true) {}
   
   // Constructor without trimming. By default it uses LocalStorage, it stores neighbors and and it does not store mass of neighbors. 
   // Measurement is automatically deleted when destructor is called.
   JetLikeEventShape(MyFunctionOfVectorOfPseudoJets<double> * measurement, double Rjet, double ptcut)
   : _measurement(measurement), _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rjet), _fcut(1.0), _trim(false),_useLocalStorage(true),_storeNeighbors(true),_storeMass(false) {}
   
   // Constructor with trimming. By default it uses LocalStorage, it stores neighbors and it does not store mass of neighbors. 
   // Measurement is automatically deleted when destructor is called.
   JetLikeEventShape(MyFunctionOfVectorOfPseudoJets<double> * measurement, double Rjet, double ptcut, double Rsub, double fcut)
   : _measurement(measurement), _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _trim(true),_useLocalStorage(true),_storeNeighbors(true),_storeMass(false) {}
   
   virtual ~JetLikeEventShape(){ if(_measurement) delete _measurement;}
   
   // This result function takes an EventStorage as input and returns the event shape value.
   // Call directly this function if you want to use a pre-built EventStorage 
   // This function is overloaded for built-in JetLikeEventShapes
   virtual double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // check if storage can be used for this shape
      if (!_check_storage_parameters(myStorage) || !myStorage.storeNeighbors() ) {
         throw Error("Storage cannot be used for this shape"); 
      } else {
         // loop through all particles
         for (unsigned int i = 0; i<myStorage.size(); i++) {
            // for particles above the cut add the measurement weighted by
            // the fraction of pT carried by particle i in a neightborhood of radius R
            if (myStorage[i].includeParticle()) myMeasurement += myStorage[i].weight()*_measurement->result(myStorage.particles_near_to(i));
         }
      }
      return myMeasurement;
   }
   
   // Standard result function. It takes a vector<PseudoJet> as input and returns the event shape.
   // It creates the appropriate EventStorage and calls result(EventStorage & myStorage)
   double result(const std::vector<PseudoJet> & particles) const {
      EventStorage myStorage(_Rjet,_ptcut,_Rsub,_fcut,_useLocalStorage,_storeNeighbors,_storeMass);
      myStorage.establishStorage(particles);		
      
      return (result(myStorage));
   }
   
   // Choose whether to use LocalStorage or not (on by default).
   // In general, no reason not turn it off except for debugging.
   void setUseLocalStorage(bool useLocalStorage) {_useLocalStorage = useLocalStorage;}      
   
   // Description
   std::string jetParameterString() const {
      std::stringstream stream;
      stream << "R_jet=" << _Rjet << ", pT_cut=" << _ptcut;
      if (_trim) stream << ", trimming with R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   
   virtual std::string description() const {
      return "Summed " + _measurement->description() + " as event shape, " + jetParameterString();
   }
   
   
protected:
   
   MyFunctionOfVectorOfPseudoJets<double> * _measurement;
   
   double _Rjet, _ptcut, _Rsub, _fcut;
   bool _trim;
   
   // optional storage array to cache the 2R x 2R partitions.
   bool _useLocalStorage;
   
   // Choose whether to store neighbors or not, and wheter to store mass of neighbors or not.
   // By default it stores neighbors, but not needed for the built-in JetLikeEventShapes
   // Also, by default it does not store mass of neighbors, but needed for built-in jet mass or mass squared shapes.
   bool _storeNeighbors,_storeMass;
   void _setStoreNeighbors(bool storeNeighbors) {_storeNeighbors=storeNeighbors;}
   void _setStoreMass(bool storeMass) {_storeMass=storeMass;}
   
   // helper to check if an external EventStorage can be used for this JetLikeEventShape (i.e. whether it has been built using the same parameters)
   bool _check_storage_parameters(EventStorage & myStorage) const {
      bool answer = false;
      if(!_trim) answer = myStorage.Rjet()==_Rjet && myStorage.ptcut()==_ptcut;
      else answer = myStorage.Rjet()==_Rjet && myStorage.ptcut()==_ptcut && myStorage.Rsub()==_Rsub && myStorage.fcut()==_fcut;
      
      return answer;
   }
   
};



//////
//
// Example Functions to Use
//
//////

//----------------------------------------------------------------------
// Below are examples of function to use as measurement in JetLikeEventShape.
// Most of these are irrevant since there are built-in JetLikeEventShapes,
// but we have maintained them for debugging purposes.

//----------------------------------------------------------------------
// FunctionUnity
// Returns 1 for each jet

class FunctionUnity : public MyFunctionOfVectorOfPseudoJets<double> {
   
   double result(const std::vector<PseudoJet> & particles) const { return 1.0; }
   
   std::string description() const { return "Jet unit weight";}
};

//----------------------------------------------------------------------
// FunctionScalarPtSum
// Finds the summed pt of the jet constituents (NOT the pT of the jet)

class FunctionScalarPtSum : public MyFunctionOfVectorOfPseudoJets<double> {
   
   double result(const std::vector<PseudoJet> & particles) const {
      double myPt = 0.0;
      for (unsigned int i = 0; i<particles.size(); i++) {
         myPt += particles[i].pt();
      }
      return myPt;
   }
   
   std::string description() const { return "Jet Scalar Pt";}
};

//----------------------------------------------------------------------
// FunctionScalarPtSumToN
// Raises FunctionScalarPtSum to an arbitrary Power

class FunctionScalarPtSumToN : public MyFunctionOfVectorOfPseudoJets<double> {
   
private:
   int _n;
   
public:
   
   FunctionScalarPtSumToN(int n) : _n(n) {}
   double result(const std::vector<PseudoJet> & particles) const {
      FunctionScalarPtSum sumPt = FunctionScalarPtSum();
      return pow(sumPt(particles),_n);
   }
   
   std::string description() const {
      std::stringstream myStream;
      myStream << "Jet Scalar Pt^" << _n;
      return myStream.str();
   }
};


//----------------------------------------------------------------------
// FunctionInvariantMass
// Finds the mass of the jet

class FunctionInvariantMass : public MyFunctionOfVectorOfPseudoJets<double> {
   
   double result(const std::vector<PseudoJet> & particles) const {
      PseudoJet myJet(0,0,0,0);
      for (unsigned int i = 0; i<particles.size(); i++) {
         myJet += particles[i];
      }
      return myJet.m();
   }
   
   std::string description() const { return "Jet Mass";}
};

//----------------------------------------------------------------------
// FunctionInvariantMassSquared
// Finds the squared mass of the jet

class FunctionInvariantMassSquared : public MyFunctionOfVectorOfPseudoJets<double> {
   
   double result(const std::vector<PseudoJet> & particles) const {
      PseudoJet myJet(0,0,0,0);
      for (unsigned int i = 0; i<particles.size(); i++) {
         myJet += particles[i];
      }
      return myJet.m2();
   }
   
   std::string description() const { return "Jet Mass^2";}
};



//////
//
// Built-In Jet-like Event Shapes
//
//////

//----------------------------------------------------------------------
// The following JetLikeEventShapes are hard coded.
// They don't use an external measurement function,
// but use information stored in EventStorage to reduce computation time


//----------------------------------------------------------------------
// ShapeJetMultiplicity
// Defines event shape jet counting

class ShapeJetMultiplicity: public JetLikeEventShape {
   
public:
   
   ShapeJetMultiplicity(double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut) {_setStoreNeighbors(false);} //Don't need to store neighbors
   ShapeJetMultiplicity(double Rjet, double ptcut, double Rsub, double fcut) : 
   JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut) {_setStoreNeighbors(false);}
   ~ShapeJetMultiplicity(){}
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else {
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()) myMeasurement += myStorage[i].weight(); //counting is just summing the weights.
         }
      }
      return myMeasurement;
   }
   
   std::string description() const {
      return "Jet multiplicity as event shape, " + jetParameterString();
   }
   
};

//----------------------------------------------------------------------
// ShapeScalarPt
// Defines event shape scalar pT sum

class ShapeScalarPt: public JetLikeEventShape{
   
public:
   
   ShapeScalarPt(double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut) {_setStoreNeighbors(false);} //Don't need to store neighbors
   ShapeScalarPt(double Rjet, double ptcut, double Rsub, double fcut) : JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut) {_setStoreNeighbors(false);}
   ~ShapeScalarPt(){}
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()) myMeasurement += myStorage[i].pt();
         }
      }
      return myMeasurement;
   }
   
   std::string description() const {
      return "Summed scalar pt as event shape, " + jetParameterString();
   }
};


//----------------------------------------------------------------------
// ShapeScalarPtToN
// Sum of pT^n of jets 

class ShapeScalarPtToN: public JetLikeEventShape {
   
public:
   
   ShapeScalarPtToN(double n, double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut), _n(n) {_setStoreNeighbors(false);} //Don't need to store neighbors
   ShapeScalarPtToN(double n, double Rjet, double ptcut, double Rsub, double fcut) : 
   JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut), _n(n) {_setStoreNeighbors(false);}
   ~ShapeScalarPtToN(){}
   
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()) myMeasurement += myStorage[i].pt()*pow(myStorage[i].pt_in_Rjet(),_n-1);
         }
      }
      return myMeasurement;
   }		
   
   
   std::string description() const {
      std::stringstream myStream;
      myStream << "Jet Scalar Pt^" << _n <<"as event shape, ";
      return myStream.str()+jetParameterString();
   }
   
protected:
   
   int _n;
};

//----------------------------------------------------------------------
// ShapeSummedMass
// Sum over jet masses 

class ShapeSummedMass: public JetLikeEventShape {
   
public:
   
   ShapeSummedMass(double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut) {
      _setStoreNeighbors(false); 
      _setStoreMass(true);} //Don't need to store neighbors, but need to store mass of neighbors
   
   ShapeSummedMass(double Rjet, double ptcut, double Rsub, double fcut) : JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut) {
      _setStoreNeighbors(false);
      _setStoreMass(true);}
   ~ShapeSummedMass(){}
   
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else if(!myStorage.storeMass()) throw Error("Mass is not stored, can't use this EventStorage"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()) myMeasurement += myStorage[i].weight()*myStorage[i].m_in_Rjet();
         }
      }
      return myMeasurement;
   }		
   
   std::string description() const {
      return "Summed jet mass as event shape, " + jetParameterString();
   }
};

//----------------------------------------------------------------------
// ShapeSummedMassSquared
// Sum over jet mass^2

class ShapeSummedMassSquared: public JetLikeEventShape {
   
public:
   
   ShapeSummedMassSquared(double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut) {
      _setStoreNeighbors(false);
      _setStoreMass(true);} //Don't need to store neighbors, but need to store mass of neighbors.
   
   ShapeSummedMassSquared(double Rjet, double ptcut, double Rsub, double fcut) : 
   JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut) {
      _setStoreNeighbors(false);
      _setStoreMass(true);}
   virtual ~ShapeSummedMassSquared(){}
   
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else if(!myStorage.storeMass()) throw Error("Mass is not stored, can't use this EventStorage"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()) myMeasurement += myStorage[i].weight()*myStorage[i].m_in_Rjet()*myStorage[i].m_in_Rjet();
         }
      }
      return myMeasurement;
   }			
   
   std::string description() const {
      return "Summed squared jet mass as event shape, " + jetParameterString();
   }
};


//----------------------------------------------------------------------
// ShapeMissingPt is pT of the vector sum of jets momenta.

class ShapeMissingPt : public JetLikeEventShape {
   
public:
   
   ShapeMissingPt(double Rjet, double ptcut): JetLikeEventShape(NULL,Rjet,ptcut) {_setStoreNeighbors(false);} //Don't need to store neighbors
   ShapeMissingPt(double Rjet, double ptcut, double Rsub, double fcut) : JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut) {_setStoreNeighbors(false);}
   ~ShapeMissingPt(){}
   
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double tot_px = 0.0;
      double tot_py = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle()){
               tot_px += myStorage[i].px();
               tot_py += myStorage[i].py();
            }
         }
      }
      return sqrt(tot_px*tot_px+tot_py*tot_py);
   }		
   
   
   std::string description() const {
      return "Missing pt as event shape, " + jetParameterString();
   }
};

//----------------------------------------------------------------------
// ShapeTrimmedSubjetMultiplicity
// A counter for subjets within trimmed fat jets.  

class ShapeTrimmedSubjetMultiplicity: public JetLikeEventShape {
   
public:
   
   ShapeTrimmedSubjetMultiplicity(double Rjet, double ptcut, double Rsub, double fcut,double ptsubcut = 0.0): 
   JetLikeEventShape(NULL,Rjet,ptcut,Rsub,fcut), _ptsubcut(ptsubcut) {_setStoreNeighbors(false);} //Don't need to store neighbors
   ~ShapeTrimmedSubjetMultiplicity(){}
   
   // Call directly this function if you want to use a pre-built EventStorage 
   double result(EventStorage & myStorage) const {
      
      double myMeasurement = 0.0;
      
      // if this result function is called with a pre-built storage need to check if storage can be used for this shape
      if(!_check_storage_parameters(myStorage)) throw Error("Storage parameters are not consistent with shape parameters"); 
      else{
         for (unsigned int i = 0; i<myStorage.size(); i++){
            if(myStorage[i].includeParticle() && myStorage[i].pt_in_Rsub() > _ptsubcut ) myMeasurement += myStorage[i].pt()/myStorage[i].pt_in_Rsub();
         }
      }
      return myMeasurement;
   }		
   
   
   std::string description() const {
      std::stringstream myStream;
      myStream << "Trimmed subjet multiplicity as event shape, ptsub_cut=" << _ptsubcut <<"and, ";
      return myStream.str()+jetParameterString();
   }
   
protected:
   
   double _ptsubcut; // additional subjet pT cut, if desired
};


//////
//
// Shape Trimming
//
//////


//----------------------------------------------------------------------
// SelectorShapeTrimming
// Event-wide trimming is implemented as a selector

Selector SelectorShapeTrimming(double Rjet, double ptcut, double Rsub, double fcut);

//----------------------------------------------------------------------
// SelectorJetShapeTrimming
// Jet-based trimming, implimented as a selector.
// It takes the jet pT as the sum of the given four-vectors
Selector SelectorJetShapeTrimming(double Rsub, double fcut);

//----------------------------------------------------------------------
// JetShapeTrimmer
// For convenience, a Transformer version of jet-shape trimming which 
// acts on individual jets.  Simply a wrapper for SelectorJetShapeTrimming

class JetShapeTrimmer : public Transformer {
   
private:
   
   double _Rsub, _fcut;
   Selector _selector;
   
public:
   
   JetShapeTrimmer(double Rsub, double fcut) : _Rsub(Rsub), _fcut(fcut){
      _selector = SelectorJetShapeTrimming(_Rsub,_fcut);
   }
   
   virtual PseudoJet result(const PseudoJet & original) const {
      return join(_selector(original.constituents()));
   }
   
   std::string jetParameterString() const {
      std::stringstream stream;
      stream << "R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   
   virtual std::string description() const {
      return "Jet shape trimmer, " + jetParameterString();
   }
};


//////
//
// Multiple pT_cut values
//
//////


//----------------------------------------------------------------------
// JetLikeEventShape_MultiplePtCutValues
// Generic class to get the whole range of event shape values as ptcut 
// is changed.

class JetLikeEventShape_MultiplePtCutValues {
   
public:
   
   JetLikeEventShape_MultiplePtCutValues(): _measurement(NULL), _useLocalStorage(true) {}
   
   // Constructor without trimming. By default it uses LocalStorage. Measurement is automatically deleted when destructor is called.
   JetLikeEventShape_MultiplePtCutValues(MyFunctionOfVectorOfPseudoJets<double> * measurement, double Rjet, double offset = 0.0): _measurement(measurement), _Rjet(Rjet), _Rsub(Rjet), _fcut(1.0), _offset(offset), _trim(false),_useLocalStorage(true) {}
   
   // Constructor with trimming. By default it uses LocalStorage. Measurement is automatically deleted when destructor is called.
   JetLikeEventShape_MultiplePtCutValues(MyFunctionOfVectorOfPseudoJets<double> * measurement, double Rjet, double Rsub, double fcut, double offset = 0.0): _measurement(measurement), _Rjet(Rjet), _Rsub(Rsub), _fcut(fcut), _offset(offset), _trim(true), _useLocalStorage(true) {}
   
   virtual ~JetLikeEventShape_MultiplePtCutValues() {if(_measurement) delete _measurement;}
   
   // Initialization: input particles and build step function.
   virtual void set_input(const std::vector<PseudoJet> & particles) {
      _storeLocalInfo(particles);	
      _buildStepFunction();
   }
   
   // Get the event shape for any value of ptcut.  
   virtual double eventShapeFor(const double ptcut_0) const;
   
   
   // Pseudo-inverse: get ptcut for a given event shape value. If initialized it uses an offset, default offset is zero.
   virtual double ptCutFor(const double eventShape_0) const;
   
   // Get the full function array. This is a vector of vector<double>. Each vector contains a ptcut/eventShape pair
   virtual const std::vector< std::vector<double> >& functionArray() const {return _functionArray;}
   
   void setUseLocalStorage(bool useLocalStorage) {_useLocalStorage = useLocalStorage;}
   
   std::string ParameterString() const {
      std::stringstream stream;
      stream << "R_jet=" << _Rjet;
      if (_trim) stream << ", trimming with R_sub=" << _Rsub <<", fcut=" << _fcut;
      stream << ", offset for inverse function=" << _offset;
      return stream.str();
   }
   
   virtual std::string description() const {
      return _measurement->description() + "as function of pT_cut, " + ParameterString();
   }
   
   
protected:
   
   // The jet measurement of interest
   MyFunctionOfVectorOfPseudoJets<double> * _measurement;
   
   double _Rjet, _Rsub, _fcut, _offset;
   bool _trim, _useLocalStorage;   
   
   // Build up the event shape as a function of ptcut.
   void _storeLocalInfo(const std::vector<PseudoJet> particles);
   void _buildStepFunction();
   std::vector< std::vector<double> > _functionArray;
};


//----------------------------------------------------------------------
// ShapeJetMultiplicity_MultiplePtCutValues
// Example of an event shape with muiltiple ptcut values.  Here, we are just doing
// jet multiplicity, but the same could be done for any of the other event shapes.
// Again ShapeJetMultiplicity_MultiplePtCutValues does not use an external measurement but it is hard-coded
// and it uses information cached in EventStorage to improve speed.

class ShapeJetMultiplicity_MultiplePtCutValues: public JetLikeEventShape_MultiplePtCutValues {
   
public:
   
   // As per the physics paper, the default offset is 0.5  
   ShapeJetMultiplicity_MultiplePtCutValues(double Rjet, double offset = 0.5): JetLikeEventShape_MultiplePtCutValues(NULL, Rjet, offset) {};
   
   ShapeJetMultiplicity_MultiplePtCutValues(double Rjet, double Rsub, double fcut, double offset = 0.5): JetLikeEventShape_MultiplePtCutValues(NULL, Rjet, Rsub, fcut, offset) {};
   
   ~ShapeJetMultiplicity_MultiplePtCutValues(){}
   
   
   void set_input(const std::vector<PseudoJet> & particles) {
      
      EventStorage myEventStorage(_Rjet,0.0,_Rsub,_fcut,_useLocalStorage,false);
      myEventStorage.establishStorage(particles);
      
      _functionArray.resize(myEventStorage.size());
      for (unsigned int i = 0; i < myEventStorage.size(); i++){ 
         _functionArray[i].push_back(myEventStorage[i].pt_in_Rjet());
         _functionArray[i].push_back(myEventStorage[i].weight());
      }	
      _buildStepFunction();
   }
   
   std::string description() const {
      return "Shape jet multiplicity as function of pT_cut, " + ParameterString();
   }
   
};



//////
//
// Njet with multiple R_jet values
//
//////


//----------------------------------------------------------------------
// ShapeJetMultiplicity_MultipleRValues
// It calculates jet multiplicity for all values of R_jet.
// It needs to store all particle mutual distances, 
// so memory required scales as N^2 for N particles in the event.
// This is the only event shape we have coded up with the ability to do
// multiple R.

class ShapeJetMultiplicity_MultipleRValues {
   
   
public:
   
   ShapeJetMultiplicity_MultipleRValues(){}
   
   // constructor without trimming
   ShapeJetMultiplicity_MultipleRValues(double ptcut):_ptcut(ptcut), _Rsub(0.0), _fcut(1.0), _trim(false) {}
   
   // constructor with trimming
   ShapeJetMultiplicity_MultipleRValues(double ptcut, double Rsub, double fcut): _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _trim(true) {}
   
   ~ShapeJetMultiplicity_MultipleRValues() {}
   
   // Initialization: input particles and build step function.
   void set_input(const std::vector<PseudoJet> & particles) {_buildStepFunction(particles);}
   
   // Get the event shape for any value of R.  
   virtual double eventShapeFor(const double Rjet_0) const;
   
   // Get the full function array. This is a vector of vectors. Each vector contains a Rjet/eventShape pair 
   const std::vector< std::vector<double> >& functionArray() const {return _functionArray;}
   
   std::string ParameterString() const {
      std::stringstream stream;
      stream << "pT_cut=" << _ptcut;
      if (_trim) stream << ", trimming with R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   
   std::string description() const {
      return "Shape jet multiplicity as function of Rjet, " + ParameterString();
   }
   
protected:
   
   double _ptcut, _Rsub, _fcut;
   bool _trim;  
   
   // Build up the event shape as a function of Rjet
   void _buildStepFunction(const std::vector<PseudoJet> particles); 
   
   std::vector< std::vector<double> > _functionArray; 
};


//////
//
// Finding Jet Axes with the Event Shape Density
//
//////


//--------------------------------------------------------
// FunctionJetAxis
// Function to find jet axis according to a given jet definition.
// Needed for the event shape densities

class FunctionJetAxis : public MyFunctionOfVectorOfPseudoJets< PseudoJet > {
   
public:
   
   FunctionJetAxis(){}
   FunctionJetAxis(fastjet::JetDefinition &jetDef) : _jetDef(jetDef) {} 
   ~FunctionJetAxis(){}
   
   PseudoJet result(const std::vector<PseudoJet> & particles) const {
      fastjet::ClusterSequence clustSeq(particles,_jetDef);
      fastjet::PseudoJet myAxis = clustSeq.inclusive_jets(0.0)[0];  // should cluster everything
      return(myAxis);
   }
   
   std::string description() const { return "Jet axis with " + _jetDef.description();}
   
private:
   
   fastjet::JetDefinition _jetDef;
};

//--------------------------------------------------------
// LightLikeAxis
// Takes a pseudojet and returns a light-like version with no mass,
// pT is set to 1, and it maintains the input rapidity and azimuth

class LightLikeAxis : public FunctionOfPseudoJet< PseudoJet > {
   
public:
   
   LightLikeAxis() {}
   
   // Convert to light-like object
   PseudoJet result(const PseudoJet & jet) const {
      PseudoJet p = jet;
      p.reset_momentum_PtYPhiM(1.0, jet.rap(), jet.phi(),0.0);
      return p;
   }
   
   std::string description() const { return "Light-like version with pT=1";}
};

//--------------------------------------------------------
// Winner-take-all recombiner.  
// Returns a recombined pseudojet whose pt is the sum of the input pts,
// but direction is the hardest jet direction.

class WinnerTakeAllRecombiner : public  fastjet::JetDefinition::Recombiner {
   
public:
   
   WinnerTakeAllRecombiner() {}
   
   virtual std::string description() const {return "WinnerTakeAll Recombination Scheme";}
   
   // recombine pa and pb and put result into pab
   virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, PseudoJet & pab) const {
      PseudoJet harder = pa;
      PseudoJet softer = pb;
      double pTa = pa.pt();
      double pTb = pb.pt();
      
      if (pTa < pTb) std::swap(harder,softer);
      
      PseudoJet direction = lightLikeVersion(harder);
      double newpT = pTa+pTb;
      pab = newpT * direction;
   }
   
protected:
   
   LightLikeAxis lightLikeVersion;
};


//--------------------------------------------------------
// EventShapeDensity_JetAxes
// This is the hybrid probability density shape for jet axes position
// It calculates axis directions using the winner-take-all recombination,
// and also returns weights

class EventShapeDensity_JetAxes {
   
public:
   
   EventShapeDensity_JetAxes(){}
   
   // using default cluster measure (anti-kT with winner-take-all recombination). By default it doesn't apply global consistency
   // and it uses LocalStorage 
   EventShapeDensity_JetAxes(double Rjet, double ptcut, bool applyGlobalConsistency = false): _Rjet(Rjet), _ptcut(ptcut), _jetDef(fastjet::JetDefinition(fastjet::antikt_algorithm, 2.0*Rjet, new WinnerTakeAllRecombiner(), fastjet::Best)), _applyGlobalConsistency(applyGlobalConsistency), _useLocalStorage(true) {_jetDef.delete_recombiner_when_unused();}
   
   // Using arbitrary jet clustering procedure
   EventShapeDensity_JetAxes(double Rjet, double ptcut, const fastjet::JetAlgorithm &jetAlgo, bool applyGlobalConsistency = false): _Rjet(Rjet), _ptcut(ptcut), _jetDef(fastjet::JetDefinition(jetAlgo, 2.0*Rjet, new WinnerTakeAllRecombiner(), fastjet::Best)), _applyGlobalConsistency(applyGlobalConsistency), _useLocalStorage(true) { _jetDef.delete_recombiner_when_unused();}
   
   ~EventShapeDensity_JetAxes(){}
   
   // Set input and find local axes	
   void set_input(const std::vector<PseudoJet> & particles) {
      _find_local_axes(particles);
      find_axes_and_weights();
   }   
   
   //Turn on/off global consistency requirement
   void setGlobalConsistencyCheck(bool applyGlobalConsistency) {_applyGlobalConsistency = applyGlobalConsistency;}
   
   //(re)find axes and weights
   void find_axes_and_weights();
   
   //Turn on/off global use of LocalStorage
   void setUseLocalStorage(bool useLocalStorage) {_useLocalStorage = useLocalStorage;}
   
   //Return a vector of PseudoJets sorted by pT. pT of each axis is its corresponding pT weight.	
   std::vector<PseudoJet> axes() const { return _distinctAxes; } 
   
   //Return the corresponding vector of Njet weights
   std::vector<double> Njet_weights() const { return _tot_Njet_weights; }
   
   //Description
   std::string ParameterString() const {
      std::stringstream stream;
      stream << "R_jet=" << _Rjet << ", pT_cut=" << _ptcut;
      stream << ". Local clustering with " << _jetDef.description()<<".";
      if (_applyGlobalConsistency) stream << " Global consistency condition on.";
      return stream.str();
   }
   
   
   std::string description() const { return "Hybrid event shape density for finding jet axes." + ParameterString();}
   
   
private:
   
   double _Rjet, _ptcut;
   fastjet::JetDefinition _jetDef;
   bool _applyGlobalConsistency;
   
   void _find_local_axes(const std::vector<PseudoJet> & particles);
   bool _isStable(const int thisAxis) const;
   
   unsigned int _N;
   std::vector<double> _tot_Njet_weights,_Njet_weights, _pt_weights;
   std::vector<int> _axes;
   
   std::vector<PseudoJet> _myParticles, _distinctAxes;
   
   LightLikeAxis _lightLikeVersion;
   
   bool _useLocalStorage;
   LocalStorage _myLocalStorage;
};


} // namespace jwj

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_JETSWITHOUTJETS_HH__
