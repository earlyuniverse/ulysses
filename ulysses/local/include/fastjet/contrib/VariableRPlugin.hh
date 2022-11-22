//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2016
//  David Krohn, Gregory Soyez, Jesse Thaler, and Lian-Tao Wang
//
//  $Id$
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

#ifndef __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__
#define __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__

#include <fastjet/internal/base.hh>
#include <fastjet/config.h>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/LimitedWarning.hh>

#include <map>
#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
  
  ////////
  //
  //  Core VR Code
  //
  ////////
   
  // In version 1.1, this replaces CoreJetAlgorithm.
  // This acts like any fastjet plugin since it implements run_clustering
  class VariableRPlugin : public JetDefinition::Plugin {

  public:
    /// Type of clustering
    ///
    /// Since version 1.2.0 of VariableR, the clustering is treated as
    /// a generalised-kt algorithm and the previous "ClusterType"
    /// argument is replaced by the "p" parameter of the
    /// generalised-kt algorithm. The definitions below are shorthand
    /// for the antikt, C/A and kt algorithm which also allow for
    /// backwards compatibility.
    /// (These are initialized in VariableRPlugin.cc.)
    static const double CALIKE;  //  =  0.0;
    static const double KTLIKE;  //  =  1.0;
    static const double AKTLIKE; //  = -1.0;

    /// for backwards compatibility reasons, we also define ClusterType
    /// as "double"
    typedef double ClusterType;

    /// The strategy to be used with the clustering
    enum Strategy{
      Best,      ///< currently N2Tiled or N2Plain for FJ>3.2.0, Native for FastJet<3.2.0
      N2Tiled,   ///< the default (faster in most cases) [requires FastJet>=3.2.0]
      N2Plain,   ///< [requires FastJet>=3.2.0]
      NNH,       ///< slower but already available for FastJet<3.2.0
      Native     ///< original local implemtation of the clustering [the default for FastJet<3.2.0]
    };

    
    /// Constructor that sets VR algorithm parameters
    ///
    ///  \param rho         mass scale for effective radius (i.e. R ~ rho/pT)
    ///  \param min_r       minimum jet radius
    ///  \param max_r       maximum jet radius
    ///  \param clust_type  whether to use CA-like, kT-like, or anti-kT-like distance measure
    ///                     (this value is the same as the p exponent in generalized-kt, with
    ///                     anti-kt = -1.0, CA = 0.0, and kT = 1.0)
    ///  \param precluster  whether to use optional kT subjets (of size min_r) for preclustering
    ///                     (true is much faster, default=false). At the moment, the only option
    ///                     for preclustering is kT (use fastjet::NestedDefsPjugin otherwise)
    ///  \param strategy    decodes which algorithm to apply for the clustering
    ///
    /// Note that pre-clustering is deprecated and will likely be
    /// removed in a future releasse of this contrib. You can get the
    /// same behaviour using the NestedDefs fastjet plugin:
    ///
    /// \code
    ///    std::list<JetDefinition> jet_defs;
    ///    jet_defs.push_back(JetDefinition(kt_algorithm, min_r));
    ///    contrib::VariableRPlugin vr_plugin(rho,min_r,max_r,clust_type);
    ///    jet_defs.push_back(JetDefinition(&vr_plugin));
    ///
    ///    NestedDefsPlugin nd_plugin(jet_defs);
    ///    JetDefinition jet_def(&nd_plugin);
    /// \endcode
    ///
    /// For FastJet>=3.2.0, we will use by default the N2Tiled
    /// strategy unless pre-clustering is requested. Otherwise, we use
    /// the "old" native strategy.
    VariableRPlugin(double rho, double min_r, double max_r, ClusterType clust_type, bool precluster = false,
                    Strategy requested_strategy = Best);
      
    // virtual function from JetDefinition::Plugin that implements the actual VR algorithm
    void run_clustering(ClusterSequence & cs) const;
      
    // information string
    virtual string description() const;
      
    // TODO:  have this return a non-trivial answer.
    virtual double R() const{ return _max_r; }
      
  private:
      
    // parameters to define VR algorithm
    double _rho2, _min_r2, _max_r, _max_r2;
    ClusterType _clust_type;
    Strategy _requested_strategy;
      
    // For preclustering, can use kT algorithm to make subclusters of size min_r
    bool _precluster;
    JetDefinition _pre_jet_def;

    // warn about pre-clustering being deprecated
    static LimitedWarning _preclustering_deprecated_warning;
    
    // helper function to decide what strategy is best
    //
    // The optimal strategy will depend on the multiplicity and _max_r
    Strategy _best_strategy(unsigned int N) const;
    
    // helper function to apply kT preclustering if desired
    void _preclustering(ClusterSequence & cs, set<int>& unmerged_jets) const;

    // native implementation of the clustering
    void _native_clustering(ClusterSequence &cs) const;

    // implementation of the clustering using FastJet NN*** classes
    template<typename NN>
    void _NN_clustering(ClusterSequence &cs, NN &nn) const;

    // Helper struct to store two jets and a distance measure
    struct JetDistancePair{
      int j1,j2;
      double distance;
    };
      
    // Helper comparitor class for comparing JetDistancePairs
    class CompareJetDistancePair {
    public:
      CompareJetDistancePair(){};
      bool operator() (const JetDistancePair & lhs, const JetDistancePair &rhs) const {
        return (lhs.distance > rhs.distance);
      }
    };
      
    // helper function to merge jet with beam
    void _merge_jet_with_beam(ClusterSequence & clust_seq, JetDistancePair & jdp, set<int>& unmerged_jets) const;
      
    // helper function to merge two jets into a new pseudojet
    void _merge_jets(ClusterSequence & clust_seq,
                     JetDistancePair & jdp,
                     priority_queue< JetDistancePair, vector<JetDistancePair>, CompareJetDistancePair > &jet_queue,
                     set<int>& unmerged_jets) const;
      
    // use ClusterMode to determine jet-jet and jet-beam distance
    inline double _get_JJ_distance_measure(const PseudoJet& j1, const PseudoJet& j2) const;
    inline double _get_JB_distance_measure(const PseudoJet& jet) const;
      
    // helper function to establish measures
    void _setup_distance_measures(ClusterSequence & clust_seq,
                                  vector<JetDistancePair> &jet_vec,
                                  set<int>& unmerged_jets) const;

  };
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VARIABLERPLUGIN_HH__
