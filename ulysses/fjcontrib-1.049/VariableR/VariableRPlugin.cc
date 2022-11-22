//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2016
//  David Krohn, Gregory Soyez, Jesse Thaler, and Lian-Tao Wang
//
//  $Id: VariableR.cc 596 2014-04-16 22:15:15Z jthaler $
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

#include "VariableRPlugin.hh"

#include <cstdio>
#include "math.h"
#include <iomanip>
#include <cmath>
#include <map>
#include <sstream>
#include <queue>

#include <fastjet/NNH.hh>
#if FASTJET_VERSION_NUMBER >= 30200
#include <fastjet/NNFJN2Plain.hh>
#include <fastjet/NNFJN2Tiled.hh>
#endif

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

  //----------------------------------------------------------------------
  // Defining static constants
  const double VariableRPlugin::CALIKE  =  0.0;
  const double VariableRPlugin::KTLIKE  =  1.0;
  const double VariableRPlugin::AKTLIKE = -1.0;
   
  //----------------------------------------------------------------------
  // classes to help run a Variable R algorithm using NN-type classes

  // class carrying particle-independent information
  class VariableRNNInfo {
  public:
    VariableRNNInfo(double rho2_in, double min_r2_in, double max_r2_in,
                    VariableRPlugin::ClusterType clust_type_in)
      : _rho2(rho2_in), _min_r2(min_r2_in), _max_r2(max_r2_in),
        _clust_type(clust_type_in) {}
    
    double rho2()  const  {return _rho2; }
    double min_r2() const {return _min_r2;}
    double max_r2() const {return _max_r2;}
    double momentum_scale_of_pt2(double pt2) const {
      return pow(pt2,_clust_type);
    }
    
  private:
    double _rho2;   ///< constant that controls the overall R magnitude
    double _min_r2; ///< minimal allowed radius
    double _max_r2; ///< maximal allowed radius
    VariableRPlugin::ClusterType _clust_type;
  };

  // class carrying the minimal info required for the clustering
  class VariableRBriefJet {
  public:
    void init(const PseudoJet & jet, VariableRNNInfo *info) {
      _rap = jet.rap();
      _phi = jet.phi();

      double pt2 = jet.pt2();

      // get the appropriate "radius"
      _beam_R2 = info->rho2()/pt2;
      if      (_beam_R2 > info->max_r2()){ _beam_R2 = info->max_r2();}
      else if (_beam_R2 < info->min_r2()){ _beam_R2 = info->min_r2();}

      // get the appropriate momentum scale
      _mom_factor2 = info->momentum_scale_of_pt2(pt2);
    }

    double geometrical_distance(const VariableRBriefJet * jet) const {
      double dphi = std::abs(_phi - jet->_phi);
      double deta = (_rap - jet->_rap);
      if (dphi > pi) {dphi = twopi - dphi;}
      return dphi*dphi + deta*deta;
    }

    double geometrical_beam_distance() const { return _beam_R2; }

    double momentum_factor() const{ return _mom_factor2; }

    /// make this BJ class compatible with the use of NNH
    double distance(const VariableRBriefJet * other_bj_jet){
      double mom1 = momentum_factor();
      double mom2 = other_bj_jet->momentum_factor();
      return (mom1<mom2 ? mom1 : mom2) * geometrical_distance(other_bj_jet);
    }
    double beam_distance(){
      return momentum_factor() * geometrical_beam_distance();
    }

    // the following are required by N2Tiled
    inline double rap() const{ return _rap;}
    inline double phi() const{ return _phi;}

  private:
    double _rap, _phi, _mom_factor2, _beam_R2;
  };

  //----------------------------------------------------------------------
  // now the implementation of VariableR itself
  //----------------------------------------------------------------------
  LimitedWarning VariableRPlugin::_preclustering_deprecated_warning;
  
  // Constructor that sets VR algorithm parameters
  //  - rho         mass scale for effective radius (i.e. R ~ rho/pT)
  //  - min_r       minimum jet radius
  //  - max_r       maximum jet radius
  //  - clust_type  whether to use CA-like, kT-like, or anti-kT-like distance measure
  //  - strategy    one of Best (the default), N2Tiled , N2Plain or NNH
  VariableRPlugin::VariableRPlugin(double rho, double min_r, double max_r,
                                   ClusterType clust_type, bool precluster,
                                   Strategy requested_strategy)
    :_rho2(rho*rho), _min_r2(min_r*min_r), _max_r(max_r), _max_r2(max_r*max_r),
     _clust_type(clust_type), _requested_strategy(requested_strategy),
     _precluster(precluster), _pre_jet_def(kt_algorithm, min_r){
    // Added errors for user input.
    if (min_r < 0.0)              throw Error("VariableRPlugin: Minimum radius must be positive.");
    if (precluster && min_r==0.0) throw Error("VariableRPlugin: To apply preclustering, minimum radius must be non-zero.");
    if (max_r < 0.0)              throw Error("VariableRPlugin: Maximum radius must be positive.");
    if (min_r > max_r)            throw Error("VariableRPlugin: Minimum radius must be bigger than or equal to maximum radius.");

    // decide the strategy

#if FASTJET_VERSION_NUMBER < 30200
    // this is only supported for the Best and Native strategies
    if ((requested_strategy!=Best) && (requested_strategy!=Native) && (requested_strategy!=NNH))
      throw Error("VariableRPlugin: with FastJet<3.2.0, Native, Best, and NNH are the only supported strategies");
#endif

    // with pre-clustering we're forced into the native clustering
    // Q: do we issue a warning? throw?
    if (precluster){
      // this is only supported for the Best and Native strategies
      if ((requested_strategy!=Best) && (requested_strategy!=Native))
        throw Error("VariableRPlugin: pre-clustering is only supported for the Native and Best strategies");

      _preclustering_deprecated_warning.warn("VariableRPlugin: internal pre-clustering is deprecated; use the NestedDefs FastJet plugin instead.");
    }  
  }

  //----------------------------------------------------------------------
  // Description of algorithm, including parameters
  string VariableRPlugin::description() const{
    stringstream myStream("");
      
    myStream << "Variable R (0903.0392), ";
      
    if (_clust_type == AKTLIKE){
      myStream << "AKT";
    } else if (_clust_type == CALIKE){
      myStream << "CA";
    } else if (_clust_type == KTLIKE){
      myStream << "KT";
    } else {
      myStream << "GenKT(p=" << _clust_type << ")";      
    }
      
    myStream << fixed << setprecision(1) << ", rho=" << sqrt(_rho2);
    myStream << ", min_r=" << sqrt(_min_r2);
    myStream << ", max_r=" << sqrt(_max_r2);
    myStream << (_precluster ? ", with precluster" : "");
    switch (_requested_strategy){
    case Best:    myStream << ", strategy=Best"; break;
    case N2Tiled: myStream << ", strategy=N2Tiled"; break;
    case N2Plain: myStream << ", strategy=N2Plain"; break;
    case NNH:     myStream << ", strategy=NNH"; break;
    case Native:  myStream << ", strategy=Native"; break;
    };
    
    return myStream.str();
  }

  //----------------------------------------------------------------------
  // do the clustering
  void VariableRPlugin::run_clustering(ClusterSequence & cs) const {
    Strategy strategy = _requested_strategy;

    // decide the best option upon request
    if (_requested_strategy==Best){
      strategy = _best_strategy(cs.jets().size());
    }

    // handle the case of native clustering
    if (strategy==Native){
      _native_clustering(cs);
      return;
    }

    // handle the NN-type clustering
    VariableRNNInfo info(_rho2, _min_r2, _max_r2, _clust_type);

#if FASTJET_VERSION_NUMBER >= 30200
    if (strategy==N2Tiled){
      NNFJN2Tiled<VariableRBriefJet,VariableRNNInfo> nnt(cs.jets(), _max_r, &info);
      _NN_clustering(cs, nnt);
    } else if (strategy==N2Plain){
      NNFJN2Plain<VariableRBriefJet,VariableRNNInfo> nnp(cs.jets(), &info);
      _NN_clustering(cs, nnp);
    } else { // NNH is the only option left
#endif
      fastjet::NNH<VariableRBriefJet,VariableRNNInfo> nnh(cs.jets(), &info);
      _NN_clustering(cs, nnh);
#if FASTJET_VERSION_NUMBER >= 30200
    }
#endif
  }

  //---------------------------------------------------------------------
  // decide the optimal strategy
  VariableRPlugin::Strategy VariableRPlugin::_best_strategy(unsigned int N) const{
#if FASTJET_VERSION_NUMBER >= 30200
    // pre-clustering requires the native implementation
    if (_precluster) return Native;
    // use the FastJet (v>/3.1) transition between N2Plain and N2Tiled
    if (N <= 30 || N <= 39.0/(max(_max_r, 0.1) + 0.6)) return N2Plain;
    return N2Tiled;
#else
    return Native;
#endif
  }

  //---------------------------------------------------------------------
  // the "new" NN-style clustering
  template<typename NN>
  void VariableRPlugin::_NN_clustering(ClusterSequence &cs, NN &nn) const{
    int njets = cs.jets().size();
    while (njets > 0) {
      int i, j, k;
      double dij = nn.dij_min(i, j);
      
      if (j >= 0) {
        cs.plugin_record_ij_recombination(i, j, dij, k);
        nn.merge_jets(i, j, cs.jets()[k], k);
      } else {
        cs.plugin_record_iB_recombination(i, dij);
        nn.remove_jet(i);
      }
      njets--;
    }
  }

  
  
  //----------------------------------------------------------------------
  // the code below is for the "native" clustering only
  
  // precluster into kT subjets if desired.
  void VariableRPlugin::_preclustering(ClusterSequence & cs, set<int>& unmerged_jets) const {
    int cntr(0);
    for(vector<PseudoJet>::const_iterator it = cs.jets().begin(); it != cs.jets().end(); it++)
      unmerged_jets.insert(unmerged_jets.end(), cntr++);
      
    // Make preclusters
    ClusterSequence pre_cs(cs.jets(), _pre_jet_def);
    vector<PseudoJet> preclustered_jets = pre_cs.inclusive_jets();
    vector<int> particle_jet_indices = pre_cs.particle_jet_indices(preclustered_jets);
      
    // This code take preclustered objects and puts them into the ClusterSequence tree of VR
    // This algorithm is likely to be deprecated in a future version of the code
    for(int i = 0 ; i < (int)preclustered_jets.size(); i++){
      queue<int> constit_indices;
      for(int j = 0 ; j < (int)particle_jet_indices.size(); j++)
        if(particle_jet_indices[j] == i)
          constit_indices.push(j);
         
      int final_jet;
      while(constit_indices.size() > 1){
        int indx1 = constit_indices.front();
        unmerged_jets.erase(indx1);
        constit_indices.pop();
        int indx2 = constit_indices.front();
        unmerged_jets.erase(indx2);
        constit_indices.pop();
        cs.plugin_record_ij_recombination(indx1, indx2, 0., final_jet);
        constit_indices.push(final_jet);
        unmerged_jets.insert(unmerged_jets.end(), final_jet);
      }
    }
  }
   
  // Implements VR alorithm (virtual function from JetDefinition::Plugin)
  void VariableRPlugin::_native_clustering(ClusterSequence & cs) const {
    set<int> unmerged_jets;
      
    if(_precluster){ // do kT preclustering
      assert(_min_r2 > 0.);  // since this is (min_r)^2, this should never happen
      _preclustering(cs, unmerged_jets);
    } else // make a list of the unmerged jets
      for(int i = 0 ; i < (int)cs.jets().size(); i++)
        unmerged_jets.insert(unmerged_jets.end(), i);
      
    // priority_queue is part of the C++ standard library.
    // The objects being sorted are JetDistancePairs
    // The comparison function is just asking who has the smallest distance
    // The distances are set initially by _setup_distance_measures and then updated by the while loop below.
    vector<JetDistancePair> jet_vec;
    _setup_distance_measures(cs, jet_vec, unmerged_jets);
    priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > jet_queue(jet_vec.begin(),jet_vec.end());
      
    // go through the jet_queue until empty
    while(!jet_queue.empty()){
         
      // find the closest pair
      JetDistancePair jdpair = jet_queue.top();
      jet_queue.pop();
         
      // Rebuild the jet_queue
      // DK - the 1.5 below was just found empirically
      // It rebuilds the jet_queue instead of letting the below functions do the hard work.
      if(jet_queue.size() > 50 && jet_queue.size() > 1.5*unmerged_jets.size()*unmerged_jets.size()){
        jet_vec.clear();
        _setup_distance_measures(cs, jet_vec, unmerged_jets);
        priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > tmp_jet_queue(jet_vec.begin(),jet_vec.end());
        swap(jet_queue,tmp_jet_queue);
      }
         
      // make sure not merged
      if((unmerged_jets.find(jdpair.j1) == unmerged_jets.end()) || (jdpair.j2 != -1 && unmerged_jets.find(jdpair.j2) == unmerged_jets.end()))
        continue;
         
         
      if(jdpair.j2 == -1) // If closest distance is to beam, then merge with beam
        _merge_jet_with_beam(cs, jdpair, unmerged_jets);
      else // Otherwise, merge jets back together
        _merge_jets(cs, jdpair, jet_queue, unmerged_jets);
    }
  }

  
  // Add final jet to clust_seq.  No need to update jet_queue, since this jet has already been deleted.
  void VariableRPlugin::_merge_jet_with_beam(ClusterSequence & clust_seq, JetDistancePair & jdp, set<int>& unmerged_jets) const{
    clust_seq.plugin_record_iB_recombination(jdp.j1, jdp.distance);
    unmerged_jets.erase(jdp.j1);
  }
   
   
  // Add jet merging to clust_seq.  Here, we need to update the priority_queue since some of the measures have changes.
  void VariableRPlugin::_merge_jets(ClusterSequence & clust_seq,
                                    JetDistancePair & jdp,
                                    priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue,
                                    set<int>& unmerged_jets) const{
      
      
    int new_jet_num;
    clust_seq.plugin_record_ij_recombination(jdp.j1,jdp.j2,jdp.distance,new_jet_num);

    unmerged_jets.erase(jdp.j1);
    unmerged_jets.erase(jdp.j2);
      
    // take the resulting jet and recompute all distances with it
    for(set<int>::iterator it = unmerged_jets.begin(); it != unmerged_jets.end(); it++){
      JetDistancePair jpair;
      jpair.j1 = new_jet_num;
      jpair.j2 = (*it);
      jpair.distance = _get_JJ_distance_measure(clust_seq.jets()[*it],clust_seq.jets()[new_jet_num]);
      jet_queue.push(jpair);
    }
    unmerged_jets.insert(unmerged_jets.end(), new_jet_num);
      
    // also add the new distance to beam
    JetDistancePair jpair;
    jpair.j1  = new_jet_num;
    jpair.j2 = -1; // -1 is for the beam
    jpair.distance = _get_JB_distance_measure(clust_seq.jets()[new_jet_num]);
    jet_queue.push(jpair);
  }
   
  // Initial distance setup
  void VariableRPlugin::_setup_distance_measures(ClusterSequence & clust_seq,
                                                 vector<JetDistancePair> &jet_vec,
                                                 set<int> & unmerged_jets) const{
    JetDistancePair jpair;
      
    // Add the jet-jet distances
    for(set<int>::iterator it1 = unmerged_jets.begin(); it1 != unmerged_jets.end(); it1++){
      for(set<int>::iterator it2 = it1; it2 != unmerged_jets.end(); it2++){
        if((*it1) != (*it2)){
          jpair.j1 = (*it1);
          jpair.j2 = (*it2);
          jpair.distance = _get_JJ_distance_measure(clust_seq.jets()[*it1],clust_seq.jets()[*it2]);
          jet_vec.push_back(jpair);
        }
      }
      // Add the jet-beam distances, and set initial merge info
      jpair.j1  = (*it1);
      jpair.j2 = -1; // -1 is for the beam
      jpair.distance = _get_JB_distance_measure(clust_seq.jets()[*it1]);
      jet_vec.push_back(jpair);
    }
  }
   
  // get the dij between two jets
  // Different measures for AKTLIKE, CALIKE, and KTLIKE, as well as for GenKT
  double VariableRPlugin::_get_JJ_distance_measure(const PseudoJet& j1, const PseudoJet& j2) const {

    // new code since version 1.2 using generic p for GenKT
    // Note that this is written to avoid having to take the "pow" twice.
    // We also keep the +-1 cases separate (for the same reason)
    double ret;
    if (_clust_type == AKTLIKE){
      ret = min(1.0/j1.perp2(), 1.0/j2.perp2());
    } else if (_clust_type == CALIKE){
      ret = 1.0;
    } else if (_clust_type == KTLIKE){
      ret = min(j1.perp2(), j2.perp2());
    } else if (_clust_type>=0){ // (include = just to be sure)
      ret = pow(min(j1.perp2(), j2.perp2()), _clust_type);
    } else {
      ret = pow(min(1.0/j1.perp2(), 1.0/j2.perp2()), -_clust_type);
    }
      
    ret *= j1.squared_distance(j2);
    return ret;
  }
   
  // jet diB between jet and beam
  // Different measures for AKTLIKE, CALIKE, and KTLIKE, as well as for GenKT
  double VariableRPlugin::_get_JB_distance_measure(const PseudoJet& jet) const{
    double pre_factor = pow(jet.perp2(), _clust_type);
    double geom_factor = _rho2 / jet.perp2();

    // While the below could implement a simplification of the pt2 factors
    // explicitly for KTLIKE, the "Native" code will be eventually deprecated
    // (and KTVR is not recommended)
    if(geom_factor < _min_r2) return _min_r2 * pre_factor;
    if(geom_factor > _max_r2) return _max_r2 * pre_factor;
    return geom_factor * pre_factor;
  }
   
   
   
} // namespace contrib

FASTJET_END_NAMESPACE
