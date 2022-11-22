// $Id: EEHelpers.hh 1292 2021-11-09 11:55:44Z scyboz $
//
// Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
// Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
//
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

#ifndef __FASTJET_CONTRIB_EEHELPERS_HH__
#define __FASTJET_CONTRIB_EEHELPERS_HH__

#include "fastjet/PseudoJet.hh"
#include <array>
#include <limits>

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//----------------------------------------------------------------------
/// Returns the 3-vector cross-product of p1 and p2. If lightlike is false
/// then the energy component is zero; if it's true the the energy component
/// is arranged so that the vector is lighlike
inline PseudoJet cross_product(const PseudoJet & p1, const PseudoJet & p2, bool lightlike=false) {
  double px = p1.py() * p2.pz() - p2.py() * p1.pz();
  double py = p1.pz() * p2.px() - p2.pz() * p1.px();
  double pz = p1.px() * p2.py() - p2.px() * p1.py();

  double E;
  if (lightlike) {
    E = sqrt(px*px + py*py + pz*pz);
  } else {
    E = 0.0;
  }
  return PseudoJet(px, py, pz, E);
}

/// Map angles to [-pi, pi]
inline const double map_to_pi(const double &phi) {
  if (phi < -M_PI)     return phi + 2 * M_PI;
  else if (phi > M_PI) return phi - 2 * M_PI;
  else                 return phi;
}

inline double dot_product_3d(const PseudoJet & a, const PseudoJet & b) {
  return a.px()*b.px() + a.py()*b.py() + a.pz()*b.pz();
}

/// Returns (1-cos theta) where theta is the angle between p1 and p2
inline double one_minus_costheta(const PseudoJet & p1, const PseudoJet & p2) {

  if (p1.m2() == 0 && p2.m2() == 0) {
    // use the 4-vector dot product. 
    // For massless particles it gives us E1*E2*(1-cos theta)
    return dot_product(p1,p2) / (p1.E() * p2.E());
  } else {
    double p1mod = p1.modp();
    double p2mod = p2.modp();
    double p1p2mod = p1mod*p2mod;
    double dot = dot_product_3d(p1,p2);

    if (dot > (1-std::numeric_limits<double>::epsilon()) * p1p2mod) {
      PseudoJet cross_result = cross_product(p1, p2, false);
      // the mass^2 of cross_result is equal to 
      // -(px^2 + py^2 + pz^2) = (p1mod*p2mod*sintheta_ab)^2
      // so we can get
      return -cross_result.m2()/(p1p2mod * (p1p2mod+dot));
    }

    return 1.0 - dot/p1p2mod;
    
  }
}

// Get the angle between two planes defined by normalized vectors
// n1, n2. The sign is decided by the direction of a vector n.
inline double signed_angle_between_planes(const PseudoJet& n1,
      const PseudoJet& n2, const PseudoJet& n) {

  // Two vectors passed as arguments should be normalised to 1.
  assert(fabs(n1.modp()-1) < sqrt(std::numeric_limits<double>::epsilon()) && fabs(n2.modp()-1) < sqrt(std::numeric_limits<double>::epsilon()));

  double omcost = one_minus_costheta(n1,n2);
  double theta;

  // If theta ~ pi, we return pi.
  if(fabs(omcost-2) < sqrt(std::numeric_limits<double>::epsilon())) {
    theta = M_PI;
  } else if (omcost > sqrt(std::numeric_limits<double>::epsilon())) {
    double cos_theta = 1.0 - omcost;
    theta = acos(cos_theta);
  } else {
    // we are at small angles, so use small-angle formulas
    theta = sqrt(2. * omcost);
  }

  PseudoJet cp = cross_product(n1,n2);
  double sign = dot_product_3d(cp,n);

  if (sign > 0) return theta;
  else          return -theta;
}

class Matrix3 {
public:
  /// constructs an empty matrix
  Matrix3() : matrix_({{{{0,0,0}},   {{0,0,0}},   {{0,0,0}}}}) {}

  /// constructs a diagonal matrix with "unit" along each diagonal entry
  Matrix3(double unit) : matrix_({{{{unit,0,0}},   {{0,unit,0}},   {{0,0,unit}}}}) {}

  /// constructs a matrix from the array<array<...,3>,3> object
  Matrix3(const std::array<std::array<double,3>,3> & mat) : matrix_(mat) {}

  /// returns the entry at row i, column j
  inline double operator()(int i, int j) const {
    return matrix_[i][j];
  }

  /// returns a matrix for carrying out azimuthal rotation
  /// around the z direction by an angle phi
  static Matrix3 azimuthal_rotation(double phi) {
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    Matrix3 phi_rot( {{ {{cos_phi, sin_phi, 0}},
	             {{-sin_phi, cos_phi, 0}},
	             {{0,0,1}}}});
    return phi_rot;
  }

  /// returns a matrix for carrying out a polar-angle
  /// rotation, in the z-x plane, by an angle theta
  static Matrix3 polar_rotation(double theta) {
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    Matrix3 theta_rot( {{ {{cos_theta, 0, sin_theta}},
	               {{0,1,0}},
	               {{-sin_theta, 0, cos_theta}}}});
    return theta_rot;
  }

  /// This provides a rotation matrix that takes the z axis to the
  /// direction of p. With skip_pre_rotation = false (the default), it
  /// has the characteristic that if p is close to the z axis then the
  /// azimuthal angle of anything at much larger angle is conserved.
  ///
  /// If skip_pre_rotation is true, then the azimuthal angles are not
  /// when p is close to the z axis.
  template<class T>
  static Matrix3 from_direction(const T & p, bool skip_pre_rotation = false) {
    double pt = p.pt();
    double modp = p.modp();
    double cos_theta = p.pz() / modp;
    double sin_theta = pt / modp;
    double cos_phi, sin_phi;
    if (pt > 0.0) {
      cos_phi = p.px()/pt;
      sin_phi = p.py()/pt;
    } else {
      cos_phi = 1.0;
      sin_phi = 0.0;
    }

    Matrix3 phi_rot({{ {{ cos_phi,-sin_phi, 0 }},
	                     {{ sin_phi, cos_phi, 0 }},
	                     {{       0,       0, 1 }} }});
    Matrix3 theta_rot( {{ {{ cos_theta, 0, sin_theta }},
	                        {{         0, 1,         0 }},
	                        {{-sin_theta, 0, cos_theta }} }});

    // since we have orthogonal matrices, the inverse and transpose
    // are identical; we use the transpose for the frontmost rotation
    // because 
    if (skip_pre_rotation) {
      return phi_rot * theta_rot;
    } else {
      return phi_rot * (theta_rot * phi_rot.transpose());
    }
  }

  template<class T>
  static Matrix3 from_direction_no_pre_rotn(const T & p) {
    return from_direction(p,true);
  }



  /// returns the transposed matrix
  Matrix3 transpose() const {
    // 00 01 02
    // 10 11 12
    // 20 21 22
    Matrix3 result = *this;
    std::swap(result.matrix_[0][1],result.matrix_[1][0]);
    std::swap(result.matrix_[0][2],result.matrix_[2][0]);
    std::swap(result.matrix_[1][2],result.matrix_[2][1]);
    return result;
  }

  // returns the product with another matrix
  Matrix3 operator*(const Matrix3 & other) const {
    Matrix3 result;
    // r_{ij} = sum_k this_{ik} & other_{kj}
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          result.matrix_[i][j] += this->matrix_[i][k] * other.matrix_[k][j];
        }
      }
    }
    return result;
  }
  
  friend std::ostream & operator<<(std::ostream & ostr, const Matrix3 & mat);
private:
  std::array<std::array<double,3>,3> matrix_;
};

inline std::ostream & operator<<(std::ostream & ostr, const Matrix3 & mat) {
  ostr << mat.matrix_[0][0] << " " << mat.matrix_[0][1] << " " << mat.matrix_[0][2] << std::endl;
  ostr << mat.matrix_[1][0] << " " << mat.matrix_[1][1] << " " << mat.matrix_[1][2] << std::endl;
  ostr << mat.matrix_[2][0] << " " << mat.matrix_[2][1] << " " << mat.matrix_[2][2] << std::endl;
  return ostr;
}

/// returns the project of this matrix with the PseudoJet,
/// maintaining the 4th component of the PseudoJet unchanged
inline PseudoJet operator*(const Matrix3 & mat, const PseudoJet & p) {
  // r_{i} = m_{ij} p_j
  std::array<double,3> res3{{0,0,0}};
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) {
      res3[i] += mat(i,j) * p[j];
    }
  }
  // return a jet that maintains all internal pointers by
  // initialising the result from the input jet and
  // then resetting the momentum.
  PseudoJet result(p);
  // maintain the energy component as it was
  result.reset_momentum(res3[0], res3[1], res3[2], p[3]);
  return result;
}

} // namespace contrib

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_EEHELPERS_HH__

