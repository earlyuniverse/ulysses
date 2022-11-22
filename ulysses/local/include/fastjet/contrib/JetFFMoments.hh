// $Id: JetFFMoments.hh 3602 2012-09-25 13:03:36Z salam $
//
// Copyright (c) 2012-, Matteo Cacciari, Paloma Quiroga-Arias, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_CONTRIB_JET_FF_MOMENTS_HH__
#define __FASTJET_CONTRIB_JET_FF_MOMENTS_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/config.h"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// \class JetFFMoments 
/// class that computes the moments of a jet fragmentation function,
/// optionally including background subtraction as described in
/// arXiv:1209.6086.
///
/// This class computes the moments of the jet fragmentation function
/// \f[
///   M_N = \frac{\sum_i p_{t,i}^N}{({\rm norm})^N}
/// \f]
/// where the sum runs over the jet constituents and several options
/// exist for the normalisation:
///
///  - by default, 'norm' is the scalar pt sum of the jet constituents.
///
///  - calling use_scalar_sum(false) causes the transverse momentum of
///    the jet to be used, instead of the scalar pt sum.
//
///  - calling set_return_numerator(true) causes the numerator only to be
///    returned i.e. norm=1.
///
///  - calling set_denominator(norm) with 'norm' positive, causes that
///    value to be used as the normalisation. This is particularly
///    useful e.g. in photon+jet events where the transverse momentum of the
///    photon is used to normalise the moments.
///
/// The values of 'N' at which one wants the computation to be performed
/// are passed as a constructor argument, and the result of this class
/// will be a vector containing the jet fragmentation moments computed
/// at these values of N.
///
/// When an optional JetMedianBackgroundEstimator is passed to the
/// constructor, the returned moments are subtracted according to the
/// method described in arXiv:1209.6086. The subtraction can include an "improvement" 
/// corrects for the effect of fluctuations when the jet spectrum
/// is some steeply function (.
///
/// Additional remarks:
///
///  - jets are required to have constituents. If subtraction is to be
///    performed, jets are also required to have area.
///
///  - FastJet version compatibility: this class is meant to be
///    compatible with FastJet 3.0 and higher. However it uses features 
///    that are still preliminary in FastJet 3.0.
///     . For usage with FastJet 3.0, please make sure that both the
///       jets and the background estimator have (or do not have)
///       explicit ghosts if you want to compute multiplicities (M_0).
///     . The interface may evolve in future versions of FastJet.
///
class JetFFMoments : public FunctionOfPseudoJet<std::vector<double> >{
public:
  /// constructor from a vector of N values
  ///  \param ns   the vector of N values
  ///  \param bge  an optional background estimator
  JetFFMoments(const std::vector<double> & ns, 
               JetMedianBackgroundEstimator *bge=0);

  /// constructor using regularly-spaced values of N
  ///  \param nmin the minimum N value 
  ///  \param nmax the maximum N value (must be larger than nmin)
  ///  \param nn   the number of N values (at least 1, nmin used if nn==1)
  ///  \param bge  an optional background estimator
  JetFFMoments(double nmin, double nmax, unsigned int nn, 
               JetMedianBackgroundEstimator *bge=0);

  /// default destructor (virtual to allow safe polymorphism)
  virtual ~JetFFMoments(){}

  // configuration handles
  //--------------------------------------------------

  /// when 'return_numerator' is set to true, just compute (and
  /// return) the numerator of the moments
  inline void set_return_numerator(bool return_numerator){
    _return_numerator=return_numerator;
  }

  /// Set a fixed value ('norm') for the denominator instead of
  /// calculating it from the jet.
  ///
  /// When 'norm' is specified (and positive) and subtraction is
  /// requested, no subtraction is applied to the denominator.
  ///
  /// A negative value for 'norm' means that the denominator is
  /// calculated internally from the jet, typically as its (scalar)
  /// transverse momentum. This is the default behaviour.
  inline void set_denominator(double norm){ _norm = norm;}

  /// Use a scalar sum of the pt of the jet consituents to calculate
  /// the denominator (when appropriate). This is the default.
  /// The alternative is to use the transverse component of the full 
  /// jet momentum. 
  inline void use_scalar_sum(bool use_scalar_sum = true){
    _use_scalar_sum = use_scalar_sum;
  }

#if (FASTJET_VERSION_NUMBER >= 30100)
  /// Turn on "improved" subtraction and provide the information that
  /// it requires.
  ///
  /// - mu is the \mu of d\sigma/dpt parametrised as sigma_0
  ///   exp(-pt/mu) in the neighbourhood of the pt of the jet (see
  ///   Section 5 of arXiv:1209.6086).
  ///
  /// If mu is negative, improved subtraction is turned off (the default).
  ///
  /// NB: with versions of FJ prior to FJ3.1, the improved subtraction
  /// also needs explicit information about the jets used for
  /// background estimation; this can be supplied with the alternative
  /// set_improved_subtraction(...) given below.
  inline void set_improved_subtraction(double mu){
    _mu = mu;
    _jets_for_improved_sub.clear();
  }
#endif

  
  /// Turn on "improved" subtraction and provide the information that
  /// it requires.
  ///
  /// - mu is the \mu of d\sigma/dpt parametrised as sigma_0
  ///   exp(-pt/mu) in the neighbourhood of the pt of the jet (see
  ///   Section 5 of arXiv:1209.6086).
  ///
  /// If mu is negative, improved subtraction is turned off (the default).
  ///
  /// The other inputs are needed so that the class can establish
  /// which jets were used in estimating the background subtraction,
  /// so that they can be used to determine the correlation between
  /// numerator and denominator of the FF moment.
  ///
  /// From FJ v 3.1 upwards, the information about the jets used for
  /// background estimation is available automatically from the bge
  /// provided in the constructor, and so this version of
  /// set_improved_subtraction is usually redundant.
  void set_improved_subtraction(double mu,
				const Selector &rho_range,
				const ClusterSequenceAreaBase &csa);


  /// Turn on "improved" subtraction and provide the information that
  /// it requires.
  ///
  /// - mu is the \mu of d\sigma/dpt parametrised as sigma_0
  ///   exp(-pt/mu) in the neighbourhood of the pt of the jet (see
  ///   Section 5 of arXiv:1209.6086).
  ///
  /// If mu is negative, improved subtraction is turned off (the default).
  ///
  /// The other parameters are needed so that the class can establish
  /// which jets were used in estimating the background subtraction,
  /// so that they can be used to determine the correlation between
  /// numerator and denominator of the FF moment.
  ///
  /// If correcting FF's for multiple jets in the same event, the
  /// version that takes a rho_range and csa will be more efficient.
  ///
  /// From FJ v 3.1 upwards, the information about the jets used for
  /// background estimation is available automatically from the bge
  /// provided in the constructor, and so this version of
  /// set_improved_subtraction is usually redundant.
  void set_improved_subtraction(double mu,
				const Selector & rho_range,
				const std::vector<PseudoJet> & particles,
				const JetDefinition  & rho_jet_def,
				const AreaDefinition & rho_area_def);

  // basic methods from FunctionOfPseudoJet (+ variants)
  //--------------------------------------------------

  /// description of the class
  virtual std::string description() const;

  // pre-declaration of subclass Info
  class Info;

  /// the computation for a given jet
  ///
  /// Requirement: the jet 'jet' needs to have constituents. If
  /// subtraction is requested (bge!=0) it also needs to have an area
  virtual std::vector<double> result(const PseudoJet &jet) const {
    Info dummy_info;
    return (*this)(jet, dummy_info);
  }

  virtual std::vector<double> operator()(const PseudoJet &jet, Info & info) const;

  // because we've introduced a new operator() with an "info" argument, we
  // also need to redefine the other operator() members
  virtual std::vector<double> operator()(const PseudoJet &jet) const {return result(jet);}
  virtual std::vector<std::vector<double> > operator()(const std::vector<PseudoJet> & jets) const {
    return FunctionOfPseudoJet<std::vector<double> >::operator()(jets);
  }

  // access to extra information
  //--------------------------------------------------
  /// return the moment N value for the i^{th} moment
  inline double N(unsigned int i) const{
    if (i>=_Ns.size())
      throw Error("JetFFMoments: out-of-range value of n requested");
    return _Ns[i];
  }

  /// return the vector of the values of N
  inline std::vector<double> Ns() const{ return _Ns;}

  // internal class to store bits of info that might be of "diagnostic"
  // interest
  class Info {
  public:
    const std::vector<double> & rhoNs()        const {return _rhoNs;}
    const std::vector<double> & sigmaNs()      const {return _sigmaNs;}
    const std::vector<double> & rNs()          const {return _rNs;}
    const std::vector<double> & ptshift_term() const {return _ptshift_term;};
    const std::vector<double> & correl_term()  const {return _correl_term ;};
    
    double rho()   const {return _rho;}
    double sigma() const {return _sigma;}

    void resize(unsigned int n) {
      _rhoNs.resize(n);
      _sigmaNs.resize(n);
      _rNs.resize(n, 0.0);
      _ptshift_term.resize(n, 0.0);
      _correl_term.resize(n, 0.0);
    }

  protected:
    std::vector<double> _rhoNs;
    std::vector<double> _sigmaNs;
    std::vector<double> _rNs;
    std::vector<double> _ptshift_term;
    std::vector<double> _correl_term;
    double _rho, _sigma;

    friend class JetFFMoments;
  };

private:
  std::vector<double> _Ns;
  JetMedianBackgroundEstimator *_bge;
  bool _return_numerator; ///< just compute the numerator
  double _norm;         ///< use a user-provided (unsubtracted) denominator
  bool _use_scalar_sum;   ///< use scalar pt sum for the denominator

  double _mu;             ///< use this value of mu for the improved_sub correction
  std::vector<PseudoJet> _jets_for_improved_sub;
  mutable Selector _rho_range_for_improved_sub; // mutable as it can take a reference

  static LimitedWarning _warnings_negative_pt;

  /// initialisation of the internal behaviour
  void _initialise();

  /// This function returns the transverse momentum to be used in the denominator
  /// of the definition of the moments of a fragmentation function for a jet, e.g.
  /// in eq.(3) of arXiv:1209.6086
  /// This can be 
  ///   1. the (subtracted) jet pt
  ///   2. the (subtracted) scalar sum of the pt of the constituents of the jet [this is the default]
  ///   3. just the number 1, in case one is interested only in the numerator.
  double _compute_normalisation(const PseudoJet &jet, 
                                const std::vector<PseudoJet> &constituents,
                                double &rho, double &sigma) const;
};

} // contrib namespace

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_JET_FF_MOMENTS_HH__
