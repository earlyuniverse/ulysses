// JetsWithoutJets Package
// Questions/Comments? danbert@mit.edu jthaler@mit.edu
//
// Copyright (c) 2013
// Daniele Bertolini and Jesse Thaler
//
// $Id: EventStorage.cc 554 2014-02-21 19:02:08Z danbert $
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

#include "EventStorage.hh"


using namespace std;


FASTJET_BEGIN_NAMESPACE

namespace jwj {


//////
//
// ParticleStorage 
//
//////

// find distance in rapidity-azimuth plane between this particle and other 
double ParticleStorage::deltaRsq(const ParticleStorage & other) const {
   double deltaRap=rap()-other.rap();
   double deltaPhi=abs(phi()-other.phi());
   if(deltaPhi > pi) deltaPhi=twopi-deltaPhi;
   return(deltaRap*deltaRap+deltaPhi*deltaPhi);
}	


//////
//
// LocalStorage
//
//////


// Creates storage array
void LocalStorage::establishStorage(const std::vector<ParticleStorage> & myParticles, double Rjet, double ptcut) {
   
   // set radius and ptcut
   _Rjet = Rjet;
   _ptcut = ptcut;
   
   // Dividing up Phase Space to bins of (approximate) size 2R by 2R
   // Actually makes bins of size bigger than 2R by 2R such that /
   // we have equal spacing of regions.
   // (This index scheme need to be synced with getRapIndex/getPhiIndex)
   _maxRapIndex = (int) floor((_rapmax)/_Rjet);// 2 rapmax / 2 Rjet
   _rapSpread = 2.0*_rapmax/(floor((_rapmax)/_Rjet));
   _maxPhiIndex = (int) floor((M_PI)/_Rjet); // 2 pi / 2 Rjet
   _phiSpread = 2.0*M_PI/(floor((M_PI)/_Rjet));
   
   // Initializing storage
   _regionStorage.resize(_maxRapIndex);
   _aboveCutBool.resize(_maxRapIndex);
   for (int rapIndex = 0 ; rapIndex < _maxRapIndex ; rapIndex++) {
      _regionStorage[rapIndex].resize(_maxPhiIndex);
      _aboveCutBool[rapIndex].resize(_maxPhiIndex);
      for (int phiIndex = 0; phiIndex < _maxPhiIndex; phiIndex++) {
         _regionStorage[rapIndex][phiIndex].clear();
      }
   }
   
   // looping through particles and assigning to bins 
   for (unsigned int i = 0; i < myParticles.size(); i++) {
      double rap = myParticles[i].rap();
      double phi = myParticles[i].phi();
      
      // find index for particle
      int lowRap = (int) floor((rap + _rapmax)/_rapSpread);
      int highRap = (int) ceil((rap + _rapmax)/_rapSpread);
      int lowPhi = (int) floor(phi/_phiSpread);
      int highPhi = (int) ceil(phi/_phiSpread);
      if (highPhi >= _maxPhiIndex) highPhi = highPhi - _maxPhiIndex; // loop around
      
      // if out of bounds, set to be in bounds.
      if (lowRap < 0) lowRap = 0;
      if (lowRap >= _maxRapIndex) lowRap = _maxRapIndex-1;
      if (highRap < 0) highRap = 0;
      if (highRap >= _maxRapIndex) highRap = _maxRapIndex-1;
      
      // Fill in storage.  For particles at periphery in rapidity, only fill two bins instead of four.  If phi overlaps, do the same thing.
      _regionStorage[lowRap][lowPhi].push_back(i);
      if (lowPhi != highPhi) _regionStorage[lowRap][highPhi].push_back(i);
      if (lowRap != highRap) _regionStorage[highRap][lowPhi].push_back(i);
      if (lowRap != highRap && lowPhi != highPhi) _regionStorage[highRap][highPhi].push_back(i);
      
   }
   
   // store whether above ptCut values
   for (int rapIndex = 0; rapIndex < _maxRapIndex; rapIndex++) {
      for (int phiIndex = 0; phiIndex < _maxPhiIndex; phiIndex++) {
         _aboveCutBool[rapIndex][phiIndex] = (getSumPt(myParticles,_regionStorage[rapIndex][phiIndex]) >= ptcut);
      }
   }
}

// What is the rapidity index?  (Sync with LocalStorage::establishStorage)
int LocalStorage::getRapIndex(const ParticleStorage & myParticle) const {
   double rap = myParticle.rap();
   int rapIndex = round((rap + _rapmax)/_rapSpread);
   
   if (rapIndex < 0) rapIndex = 0;
   if (rapIndex >= _maxRapIndex) rapIndex = _maxRapIndex - 1;
   
   return rapIndex;
}

// what is the phi index? (Sync with LocalStorage::establishStorage)
int LocalStorage::getPhiIndex(const ParticleStorage & myParticle) const {
   double phi = myParticle.phi();
   int phiIndex = round(phi/_phiSpread);
   
   if (phiIndex >= _maxPhiIndex) phiIndex = phiIndex - _maxPhiIndex;
   
   return phiIndex;
}

const vector<unsigned int> & LocalStorage::getStorageFor(const ParticleStorage & myParticle) const {
   return _regionStorage[getRapIndex(myParticle)][getPhiIndex(myParticle)];
}

bool LocalStorage::aboveCutFor(const ParticleStorage & myParticle) {
   return _aboveCutBool[getRapIndex(myParticle)][getPhiIndex(myParticle)];
}

double LocalStorage::getSumPt(const std::vector<ParticleStorage> & Particles, const std::vector<unsigned int> myIds) const {
   double myPt = 0;
   for(unsigned int i=0; i<myIds.size(); i++) myPt += Particles[myIds[i]].pt();
   return myPt;
}


//////
//
// EventStorage 
//
//////


// build a ParticleStorage for each particle and store default info,
// also build an internal vector containing particle indices (need this in _establishDerivedStorage())
void EventStorage::_establishBasicStorage(const std::vector<PseudoJet> & particles){
   
   _storage.resize(0);
   _ids.resize(0);
   
   for(unsigned int i=0; i<particles.size(); i++){
      ParticleStorage myParticleStorage(particles[i]);
      _storage.push_back(myParticleStorage);
      _ids.push_back(i); // internal vector to store particle indices
   }
}


// calculate local info and store it back in ParticleStorage.
// Right now this is calculating and storing by default pt_in_Rjet, pt_in_Rsub,m_in_Rjet,and the weight=pt/pt_in_Rjet. 
// A single EventStorage can be shared by multiple JetLikeEventShapes to reduce computation time.
// If bool flag _storeNeighbors==true it will also store a vector of indices corresponing to particles within a distance Rjet. 
// If bool flag _storeMass==true it will also store mass of the jet made of neighbors 
void EventStorage::_establishDerivedStorage() {
   
   LocalStorage myLocalStorage;
   // if requested establish LocalStorage (use info from _storage to do so)
   if(_useLocalStorage) myLocalStorage.establishStorage(_storage,_Rjet,_ptcut);
   
   // myLocalRegion defines a region of interest where to look for 
   // the neighborhood (i.e. Rjet/Rsub circles) of the particle. 
   // By default it is the whole event, it will be reduced to a smaller 2Rjet by 2Rjet region by LocalStorage
   // Note that the default vector _ids is only created once per event
   const vector<unsigned int>* myLocalRegion = &_ids;

   for (unsigned int i=0; i<_storage.size(); i++){      
      
      double pt_in_Rjet,pt_in_Rsub,m_in_Rjet;
      vector<unsigned int> neighbors;
      
      _storage[i].set_includeParticle(false);
      
      if (_useLocalStorage) { // start from rapidity/phi blocks
         if (!myLocalStorage.aboveCutFor(_storage[i])) continue; //don't do any further analysis if not above jet pTcut
         
         // if using local storage (and the region is above ptcut) then region of interest
         // where to look for the neighborhood (i.e. Rjet/Rsub circles) of the particle is 
         // given by the local storage (2R_jet by 2R_jet region)
         else myLocalRegion = &myLocalStorage.getStorageFor(_storage[i]);
      }
      
      // helper function to get all of the relevant information about neighbors
      _get_local_info(i, myLocalRegion, pt_in_Rjet, pt_in_Rsub, m_in_Rjet,neighbors);
      
      // See if there is enough pt in the neighborhood to pass _ptcut
      if (pt_in_Rjet < _ptcut) continue;
      
      assert(_Rsub <= _Rjet);
      if (pt_in_Rsub/pt_in_Rjet < _fcut) continue;
      
      _storage[i].set_includeParticle(true);
      _storage[i].set_pt_in_Rjet(pt_in_Rjet);
      _storage[i].set_pt_in_Rsub(pt_in_Rsub);
      _storage[i].set_m_in_Rjet(m_in_Rjet); // empty information if _storeMass is off
      _storage[i].set_neighbors(neighbors); // empty information if _storeNeighbors is off
      _storage[i].set_weight(_storage[i].pt()/pt_in_Rjet);
      
   }
}

// helper function to calculate local info around each particle.
// myLocalregion tells where to look, either the whole set of particles or a 2Rx2R
// local region if LocalStorage is in use
// Note that this function is calculating more information than is strictly necessary
// but this is helpful if information has to be reused.  
void EventStorage::_get_local_info(const unsigned int id, const vector<unsigned int>* myLocalRegion, double & pt_in_Rjet, double & pt_in_Rsub, double & m_in_Rjet, std::vector<unsigned int> & neighbors) const {
   
   double Rjetsq = _Rjet*_Rjet;
   double Rsubsq = _Rsub*_Rsub;
   pt_in_Rjet = 0.0;
   pt_in_Rsub = 0.0;
   m_in_Rjet = 0.0;
   neighbors.resize(0);
   PseudoJet pj_in_Rjet(0.0,0.0,0.0,0.0);
   
   
   for (unsigned int i=0; i<myLocalRegion->size(); i++) {
      double deltaRsq=_storage[id].deltaRsq(_storage[myLocalRegion->at(i)]);
      if(deltaRsq <= Rjetsq) {
         pt_in_Rjet += _storage[myLocalRegion->at(i)].pt();
         if (_storeMass) pj_in_Rjet += _storage[myLocalRegion->at(i)].pseudoJet();
         if (_storeNeighbors) neighbors.push_back(myLocalRegion->at(i));
         if (deltaRsq <= Rsubsq) pt_in_Rsub += _storage[myLocalRegion->at(i)].pt();
      }
   }
   
   m_in_Rjet = pj_in_Rjet.m();
}

} // namespace jwj

FASTJET_END_NAMESPACE
