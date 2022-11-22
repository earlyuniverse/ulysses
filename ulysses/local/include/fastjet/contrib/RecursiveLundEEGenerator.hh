// $Id: RecursiveLundEEGenerator.hh 1295 2021-11-26 12:36:26Z scyboz $
//
// Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
// Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
//
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

#ifndef __FASTJET_CONTRIB_RECURSIVELUNDEEGENERATOR_HH__
#define __FASTJET_CONTRIB_RECURSIVELUNDEEGENERATOR_HH__

#include "EEHelpers.hh"

#include <fastjet/internal/base.hh>
#include "fastjet/tools/Recluster.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include <string>
#include <vector>
#include <utility>
#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//----------------------------------------------------------------------
/// \class LundEEDeclustering
/// Contains the declustering variables associated with a single qnode
/// on the LundEE plane
class LundEEDeclustering {
public:

  /// return the pair PseudoJet, i.e. sum of the two subjets
  const PseudoJet & pair()  const {return pair_;}
  /// returns the subjet with larger transverse momentum
  const PseudoJet & harder() const {return harder_;}
  /// returns the subjet with smaller transverse momentum
  const PseudoJet & softer() const {return softer_;}


  /// returns pair().m() [cached]
  double m()         const {return m_;}

  /// returns the effective pseudorapidity of the emission [cached]
  double eta()       const {return eta_;}

  /// returns sin(theta) of the branching [cached]
  double sin_theta() const {return sin_theta_;}

  /// returns softer().modp() / (softer().modp() + harder().modp()) [cached]
  double z()         const {return z_;}

  /// returns softer().modp() * sin(theta()) [cached]
  double kt()        const {return kt_;}

  /// returns ln(softer().modp() * sin(theta())) [cached]
  double lnkt()      const {return lnkt_;}

  /// returns z() * Delta() [cached]
  double kappa()     const {return kappa_;}

  /// returns the index of the plane to which this branching belongs
  int iplane() const {return iplane_;}

  /// returns the depth of the plane on which this declustering
  /// occurred. 0 is the primary plane, 1 is the first set of leaves, etc. 
  int depth() const {return depth_;}
  
  /// returns iplane (plane index) of the leaf associated with the
  /// potential further declustering of the softer of the objects in
  /// this splitting
  int leaf_iplane() const {return leaf_iplane_;}

  /// Returns sign_s, indicating the initial parent jet index of this splitting
  int sign_s() const {return sign_s_;}
  
  /// (DEPRECATED)
  /// returns an azimuthal type angle between this declustering plane and the previous one
  /// Note: one should use psibar() instead, since we found that this definition of psi is
  /// not invariant under rotations of the event
  double psi()       const {return psi_;}

  /// update the azimuthal angle (deprecated)
  void set_psi(double psi) {psi_ = psi;}

  /// returns the azimuthal angle psibar between this declustering plane and the previous one
  double psibar()    const {return psibar_;}

  
  /// returns the coordinates in the Lund plane
  std::pair<double,double> const lund_coordinates() const {
    return std::pair<double,double>(eta_,lnkt_);
  }

  virtual ~LundEEDeclustering() {}

private:
  int iplane_;
  double psi_, psibar_, lnkt_, eta_;
  double m_, z_, kt_, kappa_, sin_theta_;
  PseudoJet pair_, harder_, softer_;
  int  depth_ = -1, leaf_iplane_ = -1;
  int sign_s_;

protected:
  /// the constructor is private, because users will not generally be
  /// constructing a LundEEDeclustering element themselves.
  LundEEDeclustering(const PseudoJet& pair,
		     const PseudoJet& j1, const PseudoJet& j2,
		     int iplane = -1, double psi = 0.0, double psibar = 0.0, int depth = -1, int leaf_iplane = -1, int sign_s = 1);

  friend class RecursiveLundEEGenerator;

};


/// Default comparison operator for LundEEDeclustering, using kt as the ordering.
/// Useful when including declusterings in structures like priority queues
inline bool operator<(const LundEEDeclustering& d1, const LundEEDeclustering& d2) {
  return d1.kt() < d2.kt();
}

//----------------------------------------------------------------------  
/// Class to carry out Lund declustering to get anything from the
/// primary Lund plane declusterings to the full Lund diagram with all
/// its leaves, etc.
class RecursiveLundEEGenerator {
 public:
  /// constructs a RecursiveLundEEGenerator with the specified depth.
  /// - depth = 0 means only primary declusterings are registered
  /// - depth = 1 means the first set of leaves are declustered
  /// - ...
  /// - depth < 0 means no limit, i.e. recurse through all leaves
  RecursiveLundEEGenerator(int max_depth = 0, bool dynamical_psi_ref = false) :
    max_depth_(max_depth), nx_(1,0,0,0), ny_(0,1,0,0), dynamical_psi_ref_(dynamical_psi_ref)
  {}

  /// destructor
  virtual ~RecursiveLundEEGenerator() {}

  /// This takes a cluster sequence with an e+e- C/A style algorithm, e.g.
  /// EECambridgePlugin(ycut=1.0).
  ///
  /// The output is a vector of LundEEDeclustering objects, ordered
  /// according to kt
  virtual std::vector<LundEEDeclustering> result(const ClusterSequence & cs) const {
    std::vector<PseudoJet> exclusive_jets = cs.exclusive_jets(2);
    assert(exclusive_jets.size() == 2);
    
    // order the two jets according to momentum along z axis
    if (exclusive_jets[0].pz() < exclusive_jets[1].pz()) {
      std::swap(exclusive_jets[0],exclusive_jets[1]);
    }

    PseudoJet d_ev = exclusive_jets[0] - exclusive_jets[1];
    Matrix3 rotmat = Matrix3::from_direction(d_ev);
    
    std::vector<LundEEDeclustering> declusterings;
    int depth = 0;
    int max_iplane_sofar = 1;
    for (unsigned ijet = 0; ijet < exclusive_jets.size(); ijet++) {

      // reference direction for psibar calculation
      PseudoJet axis = d_ev/sqrt(d_ev.modp2());
      PseudoJet ref_plane = axis;

      int sign_s = ijet == 0? +1 : -1;
      // We can pass a vector normal to a plane of reference for phi definitions
      append_to_vector(declusterings, exclusive_jets[ijet], depth, ijet, max_iplane_sofar,
                       rotmat, sign_s, exclusive_jets[0], exclusive_jets[1], ref_plane, 0., true);
    }
    
    // a typedef to save typing below
    typedef LundEEDeclustering LD;
    // sort so that declusterings are returned in order of decreasing
    // kt (if result of the lambda is true, then first object appears
    // before the second one in the final sorted list)
    sort(declusterings.begin(), declusterings.end(),
         [](const LD & d1, LD & d2){return d1.kt() > d2.kt();});

    return declusterings;
  }
  
 private:

  /// internal routine to recursively carry out the declusterings,
  /// adding each one to the declusterings vector; the primary
  /// ones are dealt with first (from large to small angle),
  /// and then secondary ones take place.
  void append_to_vector(std::vector<LundEEDeclustering> & declusterings,
                        const PseudoJet & jet, int depth,
                        int iplane, int & max_iplane_sofar,
                        const Matrix3 & rotmat, int sign_s,
                        const PseudoJet & harder,
                        const PseudoJet & softer,
                        const PseudoJet & psibar_ref_plane,
                        const double & last_psibar, bool first_time) const {
    PseudoJet j1, j2;
    if (!jet.has_parents(j1, j2)) return;
    if (j1.modp2() < j2.modp2()) std::swap(j1,j2);

    // calculation of azimuth psi
    Matrix3 new_rotmat;
    if (dynamical_psi_ref_) {
      new_rotmat = Matrix3::from_direction(rotmat.transpose()*(sign_s*jet)) * rotmat;
    } else {
      new_rotmat = rotmat;
    }
    PseudoJet rx = new_rotmat * nx_;
    PseudoJet ry = new_rotmat * ny_;
    PseudoJet u1 = j1/j1.modp(), u2 = j2/j2.modp();
    PseudoJet du = u2 - u1;
    double x = du.px() * rx.px() + du.py() * rx.py() + du.pz() * rx.pz();
    double y = du.px() * ry.px() + du.py() * ry.py() + du.pz() * ry.pz();
    double psi = atan2(y,x);

    // calculation of psibar
    double psibar = 0.;
    PseudoJet n1, n2;

    // First psibar for this jet
    if (first_time) {

      // Compute the angle between the planes spanned by (some axis,j1) and by (j1,j2)
      n1 = cross_product(psibar_ref_plane,j1);
      n2 = cross_product(j1,j2);

      double signed_angle = 0.;
      n2 /= n2.modp();
      if (n1.modp() != 0) {
        n1 /= n1.modp();
        signed_angle = signed_angle_between_planes(n1,n2,j1);
      }

      psibar = map_to_pi(j1.phi() + signed_angle);
    }
    // Else take the value of psibar_i and the plane from the last splitting to define psibar_{i+1}
    else {
      n2 = cross_product(j1,j2);
      n2 /= n2.modp();
      psibar = map_to_pi(last_psibar + signed_angle_between_planes(psibar_ref_plane, n2, j1));
    }

    int leaf_iplane = -1;
    // we will recurse into the softer "parent" only if the depth is
    // not yet at its limit or if there is no limit on the depth (max_depth<0)
    bool recurse_into_softer = (depth < max_depth_ || max_depth_ < 0);
    if (recurse_into_softer) {
      max_iplane_sofar += 1;
      leaf_iplane = max_iplane_sofar;
    }
    
    LundEEDeclustering declust(jet, j1, j2, iplane, psi, psibar, depth, leaf_iplane, sign_s);
    declusterings.push_back(declust);

    // now recurse
    // for the definition of psibar, we recursively pass the last splitting plane (normal to n2) and the last value
    // of psibar
    append_to_vector(declusterings, j1, depth, iplane, max_iplane_sofar, new_rotmat, sign_s, u1, u2, n2, psibar, false);
    if (recurse_into_softer) {
      append_to_vector(declusterings, j2, depth+1, leaf_iplane, max_iplane_sofar, new_rotmat, sign_s, u1, u2, n2, psibar, false);
    }
  }
  
  int max_depth_ = 0;
  /// vectors used to help define psi
  PseudoJet nx_;
  PseudoJet ny_;
  bool dynamical_psi_ref_;
};
  
} // namespace contrib

FASTJET_END_NAMESPACE

/// for output of declustering information
std::ostream & operator<<(std::ostream & ostr, const fastjet::contrib::LundEEDeclustering & d);

#endif  // __FASTJET_CONTRIB_RECURSIVELUNDEEGENERATOR_HH__
