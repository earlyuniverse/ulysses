//  ClusteringVetoPlugin Package
//  Questions/Comments? liew@hep-th.phys.s.u-tokyo.ac.jp
//                      stoll@hep-th.phys.s.u-tokyo.ac.jp
//
//  Copyright (c) 2014-2015
//  Seng Pei Liew, Martin Stoll
//
//  $Id: ClusteringVetoPlugin.hh 792 2015-05-04 03:42:26Z martinstoll $
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

#ifndef __FASTJET_CONTRIB_CLUSTERINGVETOPLUGIN_HH__
#define __FASTJET_CONTRIB_CLUSTERINGVETOPLUGIN_HH__

#include <fastjet/internal/base.hh>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/LimitedWarning.hh>

#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

  //------------------------------------------------------------------------
  /// \class ClusteringVetoPlugin
  /// This class implements a terminating clustering veto in the sense that
  /// if a recombination step is vetoed, both jets do not further participate
  /// in jet clustering ("turn passive").
  ///
  /// The default veto function is the mass-jump veto:
  ///
  ///  - 2 veto parameters: mu, theta
  ///  - Maximum jet radius parameter: max_r
  ///  - Metric d_ab, d_aB: CA-like, kT-like, anti-kT-like
  ///
  /// Mass-jump veto function in a recombination step j_a+j_b -> j:
  ///
  ///  - If mass m(j) < mu, cluster.
  ///  - If mass m(j) > mu and theta*m(j) > max(m(j_a),m(j_b)), veto.
  ///  - If mass m(j) > mu but theta*m(j) < max(m(j_a),m(j_b)), no veto.
  ///    Check active-passive veto, i.e.
  ///     * find passive jet j_n with smallest d_an and d_an < d_nB
  ///     * if d_an < d_ab, check veto between j_a and j_n
  ///     * if a veto is called, j_a turns passive
  ///    Do the same for j_b.
  ///  - If no active-passive veto is called, cluster.
  ///
  /// A user-defined veto function can be passed, which takes the return values
  /// { CLUSTER, VETO, NOVETO }
  /// corresponding to the first three bullet points in the mass-jump function.
  class ClusteringVetoPlugin : public JetDefinition::Plugin {

  public:
    // Type of clustering
    enum ClusterType {
      CALIKE,
      KTLIKE,
      AKTLIKE
    };

    // Result of veto function
    enum VetoResult {
      CLUSTER,
      VETO,
      NOVETO
    };
      
    // Constructor
    // mu         : veto condition 1, m(i+j) > mu
    // theta      : veto condition 2, theta * m(i+j) > max( m_i, m_j )
    // max_r      : maximum jet radius
    // clust_type : whether to use CA-like, kT-like, or anti-kT-like distance measure
    ClusteringVetoPlugin(double mu, double theta, double max_r, ClusterType clust_type);
      
    // Virtual function from JetDefinition::Plugin that implements the algorithm
    void run_clustering(fastjet::ClusterSequence & cs) const;
      
    // Information string
    virtual string description() const;
      
    // NOTE: Required by JetDefinition::Plugin
    double R() const { return sqrt(_max_r2); }

    // Set user-defined mass-jump veto
    void set_veto_function(VetoResult (*f)(const PseudoJet& j1,
					   const PseudoJet& j2)) {
      _veto_function = f; }
      
  private:

    // Parameters of MJ clustering
    double _max_r2, _mu, _theta;
    ClusterType _clust_type;

    // pointer to CheckVeto_xxx function
    VetoResult (*_veto_function) (const PseudoJet& j1, const PseudoJet& j2);

  private:
      
    // Use ClusterMode to determine jet-jet and jet-beam distance
    inline double GetJJDistanceMeasure(const PseudoJet& j1, const PseudoJet& j2) const;
    inline double GetJBDistanceMeasure(const PseudoJet& jet) const;

    // Check veto condition
    VetoResult CheckVeto(const PseudoJet& j1, const PseudoJet& j2) const;
    // Pre-defined mass-jump veto
    VetoResult CheckVeto_MJ(const PseudoJet& j1, const PseudoJet& j2) const;

  };


  /// helper class to store some clustering parameters
  class ClusteringVetoJetInfo {
  public:
    ClusteringVetoPlugin::ClusterType clust_type;
    double max_r2;
  };
  /// helper class for NNH
  class ClusteringVetoJet {
  public:
    void init(const PseudoJet & jet,
	      const ClusteringVetoJetInfo * info) {
      ph = jet.phi();
      rp = jet.rap();
      max_r2 = info->max_r2;
      switch(info->clust_type) {
      case ClusteringVetoPlugin::AKTLIKE: perpfactor = 1./jet.perp2();
	break;
      case ClusteringVetoPlugin::CALIKE: perpfactor = 1.;
	break;
      case ClusteringVetoPlugin::KTLIKE: perpfactor = jet.perp2();
	break;
      default: assert(false);
      }
    }

    double distance(const ClusteringVetoJet * jet ) const {
      double dij = min ( perpfactor, jet->perpfactor );
      double dphi = fabs(ph-jet->ph);
      if (dphi > pi) {dphi = twopi - dphi;}
      dij *= (dphi*dphi+((rp-jet->rp)*(rp-jet->rp))) / max_r2;
      return dij;
    }

    double beam_distance() const {
      return perpfactor;
    }

  private:
    double ph, rp, perpfactor, max_r2;
  };
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_CLUSTERINGVETOPLUGIN_HH__
