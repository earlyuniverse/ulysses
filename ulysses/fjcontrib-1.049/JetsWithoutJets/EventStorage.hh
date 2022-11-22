// JetsWithoutJets Package
// Questions/Comments? danbert@mit.edu jthaler@mit.edu
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: EventStorage.hh 554 2014-02-21 19:02:08Z danbert $
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

#ifndef __FASTJET_CONTRIB_EVENTSTORAGE_HH__
#define __FASTJET_CONTRIB_EVENTSTORAGE_HH__

#include "fastjet/PseudoJet.hh"
#include <sstream>

FASTJET_BEGIN_NAMESPACE

namespace jwj {


//////
//
// ParticleStorage 
//
//////


//----------------------------------------------------------------------
// ParticleStorage
// ParticleStorage stores information about single particles.
// By default it stores the full pseudojet, rapidity, phi, pt, mass, px and py.
// This information is used multiple times, so it is convenient to cache it. 
// By design, we are storing redundant information (i.e. the full pseudo jet as well as rap/phi, px/py, etc.)
// because we don't want to be ever rely on the internal way the PseudoJet stores/computes this information
// In addition, information about the neighborhood of the particle (i.e. within Rjet/Rsub circles),
// such as pT_in_Rjet, pT_in_Rsub, mass_in_Rjet, indices identifying neighbors, can be stored.
// Storing of this additional information is done externally by the EventStorage class.
// Note that the EventStorage class uses data already cached in ParticleStorage to calculate information about neighborhoods.


class ParticleStorage {

private:
   
   PseudoJet _pj;
   double _rap,_phi,_pt,_m,_px,_py,_pt_in_Rjet,_pt_in_Rsub,_m_in_Rjet,_weight;
   bool _includeParticle;
   std::vector<unsigned int> _neighbors;
   
public:
   
   ParticleStorage(){}
   
   // by default it stores the full pseudojet, rapidity, phi, pt, mass, px and py (need px and py to calculate vector pt) of the particle
   // This information is typically set by EventStorage::_establishBasicStorage().
   // Also, _pt_in_Rjet,_pt_in_Rsub,_m_in_Rjet,_weight,_includeParticle are initialized to default values. They can be reset through
   // corresponding set functions (typically done by EventStorage::_establishDerivedStorage()).
   ParticleStorage(const PseudoJet & myParticle): _pj(myParticle),_rap(myParticle.rap()),_phi(myParticle.phi()),_pt(myParticle.pt()),_m(myParticle.m()),_px(myParticle.px()), _py(myParticle.py()),_pt_in_Rjet(0.0),_pt_in_Rsub(0.0),_m_in_Rjet(0.0),_weight(0.0),_includeParticle(false) {}  
   ~ParticleStorage(){}
   
   // return stored information
   const PseudoJet & pseudoJet() const {return _pj;}	
   const double & rap() const {return _rap;}
   const double & phi() const {return _phi;}
   const double & pt() const {return _pt;}
   const double & m() const {return _m;}
   const double & px() const {return _px;}
   const double & py() const {return _py;}
   
   // deltaRsq returns the squared distance in rapidity-azimuth plane between this particle and other
   double deltaRsq(const ParticleStorage & other) const;
   
   // The following functions set and call information about the neighborhood of the particle.	
   // This information will typically be set by EventStorage::_establishDerivedStorage().
   
   // scalar sum pt in a cone of radius Rjet around this particle
   void set_pt_in_Rjet(const double & pt_in_Rjet) {_pt_in_Rjet=pt_in_Rjet;}
   const double & pt_in_Rjet() const {return _pt_in_Rjet;}
   
   // scalar sum pt in a cone of radius Rsub around this particle
   void set_pt_in_Rsub(const double & pt_in_Rsub) {_pt_in_Rsub=pt_in_Rsub;}
   const double & pt_in_Rsub() const {return _pt_in_Rsub;}
   
   // mass in a cone of radius Rjet around this particle
   void set_m_in_Rjet(const double & m_in_Rjet) {_m_in_Rjet=m_in_Rjet;}
   const double & m_in_Rjet() const {return _m_in_Rjet;}
   
   // weight=pt/pt_in_Rjet for this particle
   void set_weight(const double & weight) {_weight=weight;}
   const double & weight() const {return _weight;}
   
   // true if pt_in_Rjet >= ptcut
   void set_includeParticle(const bool & includeParticle) {_includeParticle=includeParticle;}
   const bool & includeParticle() const {return _includeParticle;}
   
   // vector of indices of particles within a distance Rjet from this particle 
   void set_neighbors(const std::vector<unsigned int> & neighbors) {_neighbors=neighbors;}
   const std::vector<unsigned int> & neighbors() const {return _neighbors;}

};


//////
//
// Local Storage (to reduce computation time)
//
//////

//----------------------------------------------------------------------
// LocalStorage
// This class divides an event into (overlapping) regions of size 2R x 2R
// and caches the partition for later use.  Partitions that are below
// the ptcut can be ignored entirely.  This saves some computational time by
// avoiding having to calculate irrelevant distances.  For speed purposes, this
// class uses ParticleStorage information (instead of PseudoJets)

class LocalStorage {

private:
   
   double _Rjet;
   double _ptcut;
   
public:
   
   LocalStorage() : _Rjet(0.0) {
      // at the moment, this value is hard coded
      _rapmax = 10.0;
   } 
   
   // make the local storage
   void establishStorage(const std::vector<ParticleStorage> & myParticles, double Rjet, double ptcut);
   
   // give the 2R x 2R region relevant for a given particle
   const std::vector<unsigned int> & getStorageFor(const ParticleStorage & myParticle) const;
   
   // give whether the 2R x 2R region has enough pT to merit further examination 
   bool aboveCutFor(const ParticleStorage & myParticle);
   
private:
      
   // storage for partitions of particles
   std::vector<std::vector<std::vector<unsigned int> > > _regionStorage; 
   
   // bool for whether region needs to be considered
   std::vector<std::vector<bool> > _aboveCutBool;
   
   double _rapmax;
   int _maxRapIndex;
   double _rapSpread;
   int _maxPhiIndex;
   double _phiSpread;
   
   // helper functions to determine which region a pseudo jet is in   
   int getRapIndex(const ParticleStorage & myParticle) const;
   int getPhiIndex(const ParticleStorage & myParticle) const;
   
   // helper function to determine sum scalar pt of a subset of particles identifed by myIds 
   double getSumPt(const std::vector<ParticleStorage> & Particles, const std::vector<unsigned int> myIds) const;
   
};


//////
//
// EventStorage 
//
//////

//----------------------------------------------------------------------
// EventStorage
// This class is a wrapper for a vectors of ParticleStorages.
// It uses default information stored in ParticleStorage to calculate additional infomation
// such pt_in_Rjet,pt_in_Rsub etc and it stores this information back to ParticleStorage.
// Can choose if want to use LocalStorage to reduce computation time (on by default)

class EventStorage {
   
private:
   
   double _Rjet,_ptcut,_Rsub,_fcut;
   std::vector<unsigned int> _ids;	
   std::vector<ParticleStorage> _storage;
   bool _useLocalStorage;
   bool _storeNeighbors, _storeMass;

public:
   
   EventStorage(){}
   // standard constructor, by default it uses LocalStorage and for each particle it stores particles within a distance Rjet (neighbors) and it does
   // not store mass of the jet made of neighbors
   EventStorage(double Rjet, double ptcut, double Rsub, double fcut, bool useLocalStorage=true, bool storeNeighbors=true, bool storeMass=false):
   _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rsub), _fcut(fcut), _useLocalStorage(useLocalStorage), _storeNeighbors(storeNeighbors), _storeMass(storeMass) {}

   // constructor with Rsub and fcut automatically initialized to Rjet and 1.0 respectively. 
   // By default it uses LocalStorage and for each particle it stores particles within a distance Rjet (neighbors) and it does
   // not store mass of the jet made of neighbors.
   // Use this constructor if don't need to do trimming.	
   EventStorage(double Rjet, double ptcut, bool useLocalStorage=true, bool storeNeighbors=true, bool storeMass=false):
   _Rjet(Rjet), _ptcut(ptcut), _Rsub(Rjet), _fcut(1.0), _useLocalStorage(useLocalStorage), _storeNeighbors(storeNeighbors), _storeMass(storeMass) {}
   ~EventStorage() {}
   
   // establish the storage
   void establishStorage(const std::vector<PseudoJet> & particles){
      _establishBasicStorage(particles);
      _establishDerivedStorage();
   }
   
   // get general info about the storage
   unsigned int size() const {return _storage.size();}
   const double & Rjet() const {return _Rjet;}
   const double & ptcut() const {return _ptcut;}
   const double & Rsub() const {return _Rsub;}
   const double & fcut() const {return _fcut;}
   const bool & storeNeighbors() const {return _storeNeighbors;}
   const bool & storeMass() const {return _storeMass;}
   
   // access individual ParticleStorage through [] operator
   ParticleStorage operator[](const unsigned int i) const {return _storage[i];}
   
   // convert stored indices of neighbors (particles within a distance Rjet) into a vector<PseudoJet>
   std::vector<PseudoJet> particles_near_to(const unsigned int id) const {
      std::vector<unsigned int> neighbors=_storage[id].neighbors();
      std::vector<PseudoJet> answer;
      for(unsigned int i=0; i<neighbors.size(); i++) answer.push_back(_storage[neighbors[i]].pseudoJet());
      return answer;
   }
   
   // Description
   std::string parameterString() const {
      std::stringstream stream;
      stream << "R_jet=" << _Rjet << ", pT_cut=" << _ptcut << ", R_sub=" << _Rsub <<", fcut=" << _fcut;
      return stream.str();
   }
   std::string description() const {
      return "Event Storage with "+parameterString();
   }
   
private:
      
   void _establishBasicStorage(const std::vector<PseudoJet> & particles);
   void _establishDerivedStorage();
   void _get_local_info(const unsigned int id, const std::vector<unsigned int>* myLocalRegion, double & pt_in_Rjet, double & pt_in_Rsub, double & m_in_Rjet, std::vector<unsigned int> & neighbors) const;
   
};

} // namespace jwj

FASTJET_END_NAMESPACE

#endif
