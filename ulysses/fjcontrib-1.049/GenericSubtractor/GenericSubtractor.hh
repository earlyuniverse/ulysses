// $Id: GenericSubtractor.hh 861 2015-09-21 10:35:22Z gsoyez $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of the GenericSubtractor package of FastJet Contrib.
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

#ifndef __FASTJET_CONTRIB_GENERIC_SUBTRACTOR_HH__
#define __FASTJET_CONTRIB_GENERIC_SUBTRACTOR_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"

#include "ShapeWithComponents.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

// forward declaration of assistional class(es) defined below.
class GenericSubtractorInfo;


//------------------------------------------------------------------------
/// \class GenericSubtractor
///
/// A class to perform subtraction of background (eg pileup or
/// underlying event) from a jet shape.
///
/// This class is a tool that allows one to subtract pileup from jet shapes
/// (i.e. a FunctionOfPseudoJet<double>). It implements the method
/// described in:
/// Gregory Soyez, Gavin P. Salam, Jihun Kim, Souvik Dutta and Matteo Cacciari,
/// Phys. Rev. Lett. 110, 162001 (2013) [arXiv:1211.2811]
///
/// The basic usage of this class is as follows:
///
///   GenericSubtractor gensub(& some_background_estimator);
///   FunctionOfPseudoJet<double> shape;
///   double subtracted_shape_value = gensub(shape, jet);
///
/// By default, this only subtracts the transverse momentum density
/// rho (obtained from the background estimator provided to the
/// class).
///
/// Two options allow also for the inclusion of the "particle mass"
/// contribution to the background density (the "rho_m" term): one can
/// either instruct the GenericSubtractor class to compute rho_m from
/// the same background estimator using
/// 
///   gensub.set_common_bge_for_rho_and_rhom();                   (1)
/// 
/// or explicitly construct gensub with two background estimators
/// 
///   GenericSubtractor gensub(& background_estimator_for_rho,
///                            & background_estimator_for_rhom);  (2)
///
/// Note that since FastJet 3.1, the option (1) will work with any
/// background estimator that internally computes rho_m (and has that
/// computation enabled). For FastJet 3.0, the first option is only
/// available for JetMedianBackgroundEstimator.
///
/// For the option (2), 'gensub' will obtain rho_m using
///
///   background_estimator_for_rhom->rho()   [NOT rho_m()!]  
///
/// unless one calls gensub.set_use_bge_rhom_rhom() (it requires
/// FJ>=3.1) in which case it will be obtained using
///
///   background_estimator_for_rhom->rhom()
///
/// The right choice between these two procedures for option (2)
/// depends on the details of 'background_estimator_for_rhom'.
///
///
/// If the background transverse momentum density rho (and optionally the 
/// density of sqrt(pt^2+m^2) - pt, denoted by rhom) are known, the
/// subtractor can be constructed directly as
///
///   GenericSubtractor gensub(rho,rhom);
///
/// Extra information can be retrieved using
///
///   GenericSubtractorInfo info;
///   double subtracted_shape_value = gensub(shape, jet, info);
///
class GenericSubtractor{
public:
  /// default constructor
  /// leaves the object unusable
  GenericSubtractor() : 
    _bge_rho(0), _bge_rhom(0), _jet_pt_fraction(0.01),
    _common_bge(false), _rhom_from_bge_rhom(false), 
    _rho(_invalid_rho), _externally_supplied_rho_rhom(false) {}

  /// Constructor that takes a pointer to a background estimator for
  /// rho and optionally a pointer to a background estimator for
  /// rho_m.  If the latter is not supplied, rho_m is assumed to
  /// always be zero (this behaviour can be changed by calling
  /// set_common_bge_for_rho_and_rhom).
  /// See also the discussion above (in particular for the use of
  /// bge_rhom)
  GenericSubtractor(BackgroundEstimatorBase *bge_rho,
		    BackgroundEstimatorBase *bge_rhom=0) :
    _bge_rho(bge_rho), _bge_rhom(bge_rhom), _jet_pt_fraction(0.01),
    _common_bge(false), _rhom_from_bge_rhom(false), 
     _rho(_invalid_rho), _externally_supplied_rho_rhom(false)  {}

  /// Constructor that takes an externally supplied value for rho and,
  /// optionally, for rho_m. The latter defaults to zero.
  GenericSubtractor(double rho, double rhom=0);
  

  /// destructor
  ~GenericSubtractor(){}

  /// a description of what this class does
  std::string description() const;

  /// returns the estimate of the supplied shape, for the given jet,
  /// after background subtraction.
  double operator()(const FunctionOfPseudoJet<double> &shape,
		    const PseudoJet &) const;

  /// returns the estimate of the supplied shape, for the given jet,
  /// after background subtraction. It also sets the "info" variable
  /// with information about the details of the subtraction (e.g. the
  /// shape's derivatives wrt to the amount of pileup, and the background
  /// densities used).
  double operator()(const FunctionOfPseudoJet<double> & shape,
		    const PseudoJet & jet, GenericSubtractorInfo & info) const;

  /// when only one background estimator ('bge_rho') is specified,
  /// calling this routine with argument "true", causes rho_m to be be
  /// calculated from the same background estimator as rho, instead of
  /// being set to zero. 
  ///
  /// Since FastJet 3.1, this will use the rho_m calculation from the
  /// background estimator if it is available [the default for both
  /// GridMedianBackgroundEstimator and JetMedianBackgroundEstimator].
  /// In all other cases, the use of a common background estimator for
  /// rho and rho_m only works if the estimator is a
  /// JetMedianBackgroundEstimator (or derived from it), [since it
  /// makes use of that class's set_jet_density_class(...) facility].
  /// 
  void set_common_bge_for_rho_and_rhom(bool value=true);

  /// equivalent to 'set_common_bge_for_rho_and_rhom'
  ///
  /// NOTE: this is DEPRECATED and kept only for backwards
  /// compatibility reasons. Use 'set_common_bge_for_rho_and_rhom'
  /// instead.
  void use_common_bge_for_rho_and_rhom(bool value=true);

  /// returns true if it uses a common background estimator for both
  /// rho and rho_m (see set_common_bge_for_rho_and_rhom for
  /// details)
  bool common_bge_for_rho_and_rhom() const { return _common_bge; }
  
  /// when the GenericSubtractor has been created with two background
  /// estimators (one for rho, the second for rho_m), setting this to
  /// true will result in rho_m being estimated using
  /// bge_rhom->rho_m() instead of bge_rhom->rho()
  /// See also the discussion at the top of the class.
  void set_use_bge_rhom_rhom(bool value=true);

  /// returns true if rho_m is being estimated using
  /// bge_rhom->rho_m(). (s
  bool use_bge_rhom_rhom() const { return _rhom_from_bge_rhom; }
  
  /// sets the fraction of the jet pt that will be used in the
  /// stability condition when computing the derivatives
  void set_jet_pt_fraction_for_stability(double jet_pt_fraction){
    _jet_pt_fraction=jet_pt_fraction;
  }



protected:
  // tools to help compute the derivatives
  //----------------------------------------------------------------------
  /// do the computation of the various derivatives
  /// 
  /// rhoA_pt_fraction is rho_p/(rho_p + rho_m) used to know along
  /// which trajectory in the (rho_p, rho_m) plane one must estimate
  /// the derivative.
  void _compute_derivatives(const FunctionOfPseudoJet<double> &shape,
                            const PseudoJet &jet,
                            const double original_ghost_scale,
                            const double ghost_area,
                            const double f0, 
			    const double rho_pt_fraction,
			    GenericSubtractorInfo &info) const;

  /// make a sweep in stepsize to determine which step is the most precise
  ///
  /// Note that this returns the values of the function at the points
  /// needed to compute the derivatives
  double _optimize_step(const FunctionOfPseudoJet<double> &shape,
                        const PseudoJet &jet,
                        const double original_ghost_scale,
                        const double ghost_area,
                        const double x_fraction,
                        const double f0, 
			double cached_functions[4],
                        const double max_step) const;

  /// the function that does the rescaling
  double _shape_with_rescaled_ghosts(const FunctionOfPseudoJet<double> &shape,
				     const PseudoJet &jet,
				     const double original_ghost_scale,
				     const double new_ghost_scale, 
				     const double new_dmass_scale) const;

  /// if the shape is a ShapeWithComponents, apply the subtraction on
  /// each of the components
  double _component_subtraction(const ShapeWithComponents * shape_ptr, 
				const PseudoJet & jet,
				GenericSubtractorInfo &info) const;

  BackgroundEstimatorBase *_bge_rho, *_bge_rhom;
  double _jet_pt_fraction;
  bool _common_bge, _rhom_from_bge_rhom;
  double _rho, _rhom;
  bool _externally_supplied_rho_rhom;
  /// a value of rho that is used as a default to label that the stored
  /// rho is not valid for subtraction. 
  //
  // NB: there are two reasons for not having the value written here:
  // 1) that it caused problems on karnak with g++ 4.0.1 and 2) that
  // we anyway like -infinity as a default, and since that's a function,
  // that's not allowed in an include file.
  static const double _invalid_rho;

  /// deprecated "use_common_bge_for_rho_and_rhom"
  static LimitedWarning _warning_depracated_use_common_bge;
  static LimitedWarning _warning_unused_rhom;
};

//------------------------------------------------------------------------
/// \class GenericSubtractorInfo
///
/// Helper class that allows to get extra information about the shape
/// subtraction
class GenericSubtractorInfo{
public:
  /// returns the unsubtracted shape
  double unsubtracted() const{ return _unsubtracted;}

  /// returns the result of applying the subtraction including only the
  /// first-order correction
  double first_order_subtracted() const{ return _first_order_subtracted;}
  
  /// returns the result of applying the subtraction including the
  /// first and second-order corrections (this is the default returned
  /// by the subtractions).
  double second_order_subtracted() const{ return _second_order_subtracted;}
  
  /// returns an estimate of the 3rd order correction. Note that the
  /// third derivative is quite likely to be poorly evaluated in some
  /// cases. For this reason, it is not used by default.
  double third_order_subtracted() const{ return _third_order_subtracted;}

  /// returns the first derivative (wrt to rho+rhom). Derivatives
  /// are taken assuming a constant value for the ratio of rhom/rho, as
  /// determined from the background estimates for the jet.
  ///
  /// Note: derivatives are not evaluated for shapes with components,
  /// and therefore set to zero.
  double first_derivative() const{ return _first_derivative;}

  /// returns the second derivative (wrt to rho+rhom). It is not evaluated
  /// for shapes with components, and therefore set to zero.
  double second_derivative() const{ return _second_derivative;}

  /// returns an estimate of the 3rd derivative (wrt to rho+rhom);
  /// note that the third derivative is quite likely to be poorly
  /// evaluated in some cases. It is not used by default. It is not
  /// evaluated for shapes with components, and therefore set to zero.
  double third_derivative() const{ return _third_derivative;}

  /// returns the ghost scale actually used for the computation
  double ghost_scale_used() const{ return _ghost_scale_used;}

  /// return the value of rho and rhom used in the subtraction
  double rho() const  { return _rho;  }
  double rhom() const { return _rhom; }
  
protected:
  double _unsubtracted;
  double _first_order_subtracted;
  double _second_order_subtracted;
  double _third_order_subtracted;
  double _first_derivative, _second_derivative, _third_derivative;
  double _ghost_scale_used;
  double _rho, _rhom;

  friend class GenericSubtractor;
};


}

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_GENERIC_SUBTRACTION_HH__
