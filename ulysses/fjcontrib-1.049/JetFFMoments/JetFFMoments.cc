// $Id: JetFFMoments.cc 3602 2012-09-25 13:03:36Z salam $
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

#include "fastjet/ClusterSequenceArea.hh"
#include "JetFFMoments.hh"
#include "fastjet/config.h"
#include <ostream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  using namespace std;

#if (FASTJET_VERSION_NUMBER < 30100)
  /// \class BackgroundJetScalarPtDensity
  /// Class that implements (scalar pt sum of jet)/(scalar area of jet)
  /// for background estimation <i>(this is a preliminary class)</i>.
  ///
  /// Optionally it can return a quantity based on the sum of pt^n,
  /// e.g. for use in subtracting fragementation function moments.
  class BackgroundJetScalarPtMomentDensity : public FunctionOfPseudoJet<double> {
  public:
    /// Constructor to provide background estimation based on 
    /// \f$ sum_{i\in jet} p_{ti}^{n} \f$
    BackgroundJetScalarPtMomentDensity(double n=1) : _pt_power(n) {}
    
    virtual double result(const PseudoJet & jet) const{
      // do not include the ghosts in the list of constituents to have a
      // correct behaviour when _pt_power is <= 0
      vector<PseudoJet> constituents = (!SelectorIsPureGhost())(jet.constituents());
      double scalar_pt = 0;
      for (unsigned i = 0; i < constituents.size(); i++) scalar_pt += pow(constituents[i].perp(), _pt_power);
      return scalar_pt / jet.area();
    }
    
    virtual std::string description() const {return "BackgroundScalarJetPtMomentDensity";}
    
  private:
    double _pt_power;
  };
#endif

  LimitedWarning JetFFMoments::_warnings_negative_pt;

  // ctor from a vector of n values
  //  - ns   the vector of n values (all non-negative)
  //  - bge  an optional background estimator
  JetFFMoments::JetFFMoments(const vector<double> & ns, 
                             JetMedianBackgroundEstimator *bge){
    _Ns = ns;
    _bge = bge;
    _initialise();
  }

  // ctor using regularly-spaced values of n
  //  - nmin the minimal n value (must be non-negative)
  //  - nmax the maximal n value (must be larger than nmin)
  //  - nn   the number of n values (at least 1, nmin used if nn==1)
  //  - bge  an optional background estimator
  JetFFMoments::JetFFMoments(double nmin, double nmax, unsigned int nn, 
                             JetMedianBackgroundEstimator *bge){

    if (nn==0){
      throw Error("JetFFMoments should be constructed with at least one element");
    }

    _Ns.resize(nn);
    if (nn==1){ _Ns[0] = nmin;}
    else {
      for (unsigned int in=0; in<nn; in++){
        _Ns[in] = nmin + in*(nmax-nmin)/(nn-1);
      }
    }

    _bge = bge;
    _initialise();
  }

  // configuration handles
  //--------------------------------------------------
  void JetFFMoments::set_improved_subtraction(
 	           double mu,
                   const Selector &rho_range,
                   const std::vector<PseudoJet> & particles,
                   const JetDefinition &jet_def,
                   const AreaDefinition &area_def){
    _mu = mu;
    ClusterSequenceArea *csa = new ClusterSequenceArea(particles, jet_def, area_def);
    _jets_for_improved_sub = csa->inclusive_jets();
    _rho_range_for_improved_sub = rho_range;
    csa->delete_self_when_unused();
  }

  void JetFFMoments::set_improved_subtraction(
                   double mu,
                   const Selector &rho_range,
                   const ClusterSequenceAreaBase &csa){
    _mu = mu;
    _jets_for_improved_sub = csa.inclusive_jets();
    _rho_range_for_improved_sub = rho_range;
  }

  // basic methods from FunctionOfPseudoJet
  //--------------------------------------------------

  /// description of the class
  std::string JetFFMoments::description() const{
    ostringstream oss;
    if (_return_numerator) oss << "Numerator of the ";
    oss << "Jet fragmentation function moments calculated";
    if (!_return_numerator){
      if (_norm>0){
        oss << " with a fixed denominator";
      } else if (_use_scalar_sum){
        oss << " using the scalar pt sum as denominator";
      } else {
        oss << " using the pt of the jet as denominator";
      }
    }
    if (_bge){
      oss << ", with background subtracted using the estimator " << _bge->description();
    }
    if ( _mu > 0 ) {
#if (FASTJET_VERSION_NUMBER >= 30100)
      if (_jets_for_improved_sub.size() == 0)
	oss << ", subtraction improved using jets from the background estimator and mu = " << _mu; 
      else 
#endif
      oss << ", subtraction improved using jets in the range " << _rho_range_for_improved_sub.description() << " and mu = " << _mu;
    }
    oss << ".";
    return oss.str();
  }

  /// the computation for a given jet
  ///
  /// Requirement: the jet 'jet' need to have constituents. If
  /// subtraction is requested (bge!=0) it also need to have an area
  vector<double> JetFFMoments::operator()(const PseudoJet &jet, Info & info) const{
    // check if we have constituents
    if (! jet.has_constituents())  //Q: to jets or on jets?
      throw Error("JetFFMoments can only be applied to jets having constituents");

    // if subtraction is requested, check if we have an area
    if ((_bge) && (! jet.has_area()))
      throw Error("JetFFMoments with background subtraction can only be applied to jets having an area");

    vector<double> ffm(_Ns.size(), 0.0);     // the vector for the results

    // set up size of our complementary info
    info.resize(_Ns.size());

    //   When areas are evaluated with explicit ghosts, these contribute to the
    //   multiplicity. This is fine as long as both the jet and the bge
    //   have (or do not have) explicit ghosts, but would lead to wrong
    //   results (for the multiplicity) if only one of the 2 has them.
    //   A possible way out is to discard the ghosts from all
    //   sums over constituents. This must be done in the line below
    //   as well as in BackgroundJetScalarPtDensity in FastJet.
#if (FASTJET_VERSION_NUMBER >= 30100)
    vector<PseudoJet> constituents = (!SelectorIsPureGhost())(jet.constituents());
#else // In FJ versions prior to 3.1, ghosts are not discarded in the density
      // class and so should be retained here too to avoid a mismatch in their
      // subtraction from the multiplicity.
    vector<PseudoJet> constituents = jet.constituents();
#endif

    // Note: the code below traverses the constituents in their
    // original order. For better precision (at large N) it may be
    // useful to browse them from softest to hardest

    double rho, sigma;
    // This is the (subtracted) denominator that goes in the
    // definition of the moments.  It's called S1 to follow the
    // notation of eq.(10) of 1209.6086, but _compute_normalisation()
    // may actually return something different (e.g. a constant 1)
    double S1 = _compute_normalisation(jet, constituents, rho, sigma);
    info._rho = rho;
    info._sigma = sigma;

    if (S1<=0){
      _warnings_negative_pt.warn("JetFFMoments: Negative or zero (subtracted) denominator. Returning 1 for all moments.");
      for (unsigned int i=0;i<ffm.size();i++) ffm[i]=1.0;
      return ffm;
    }
      
    // now do compute the moments (numerator only)
    for (unsigned int i=0; i<constituents.size(); i++){
      double pti = constituents[i].pt();
      for (unsigned int in=0; in<_Ns.size(); in++){
        ffm[in] += pow(pti, _Ns[in]);
      }
    }

    // apply subtraction to the moments (numerator only), if requested
    const FunctionOfPseudoJet<double> * old_density_class = 0;
    if (_bge){
      old_density_class = _bge->jet_density_class();
      for (unsigned int in=0; in<_Ns.size(); in++){
#if (FASTJET_VERSION_NUMBER >= 30100)
        BackgroundJetScalarPtDensity scalar_density(_Ns[in]);
        _bge->set_jet_density_class(&scalar_density);
#else
        BackgroundJetScalarPtMomentDensity scalar_density(_Ns[in]);
        _bge->set_jet_density_class(&scalar_density);
#endif
        info._rhoNs[in] = _bge->rho(jet);
        info._sigmaNs[in] = _bge->sigma(jet);
        ffm[in] -= info._rhoNs[in] * jet.area();
      }
    }

    // apply the normalisation to the moments
    for (unsigned int in=0; in<_Ns.size(); in++){
      ffm[in] /= pow(S1, _Ns[in]);
    }

    // apply the improved_sub correction if requested
    // See Eqs. (18) and (25) in arXiv:1209.6086
    if (_mu>0){
      if ((_norm>0) || (_return_numerator))
        throw Error("JetFFMoments: improved subtraction (mu>0) is not available for a fixed denominator");

      // get the jets that have been used for the computation of rho
      vector<PseudoJet> jets_used;
      if (_jets_for_improved_sub.size() == 0) {
#if (FASTJET_VERSION_NUMBER >= 30100)
        jets_used = _bge->jets_used();
#else 
        throw Error("jets_for_improved_sub was empty; this is inconsistent when improved subtraction is requested (mu>0)");
#endif       
      } else {
        if (_rho_range_for_improved_sub.takes_reference()) _rho_range_for_improved_sub.set_reference(jet); 
        jets_used = _rho_range_for_improved_sub(_jets_for_improved_sub);
      }

      // compute the relevant sums
      double sum_qt=0.0;
      double sum_qt2=0.0;
      vector<double> sum_Qn(_Ns.size(), 0.0);
      vector<double> sum_QnQn(_Ns.size(), 0.0);
      vector<double> sum_Qq(_Ns.size(), 0.0);

      // loop over jets used for estimating the background
      for (unsigned int ijet=0; ijet<jets_used.size(); ijet++){
        const PseudoJet & current_jet = jets_used[ijet];
        double qt = 0.0;
        vector<PseudoJet> constits = (!SelectorIsPureGhost())(current_jet.constituents());
        if (!_use_scalar_sum){
          qt = (current_jet - rho*current_jet.area_4vector()).pt();

          // we need to make sure we have both the +ve and -ve
          // "subtracted" pt. The prescription we use here is to flip
          // the sign whenever the transverse component of the
          // subtracted 4-vector is larger than the transverse
          // component of the original jet.
          if (current_jet.pt2() < rho*rho*current_jet.area_4vector().pt2()) { qt *= -1.0; }
        } else {
          for (unsigned int ic=0; ic<constits.size(); ic++) qt+=constits[ic].pt();
          qt -= rho * current_jet.area();
        }

        sum_qt  += qt;
        sum_qt2 += qt*qt;

        // sum_Qn =  ( sum_(i in jet) constit_i.pt()^N ) - rho_N * jet.area
        for (unsigned int in=0; in<_Ns.size(); in++){
          double Qn=0.0;
          for (unsigned int ic=0; ic<constits.size(); ic++) {
            Qn += pow(constits[ic].pt(), _Ns[in]);
	  }  
          Qn -= info._rhoNs[in] * current_jet.area();
          sum_Qn[in] += Qn;
          sum_QnQn[in] += Qn*Qn;
          sum_Qq[in] += Qn*qt;
        } // end loop over N
      } // end loop over jets


      // the average expected qt fluctuation in the jet, given
      // the knowledge of the background fluctuations (characterised by sigma) and 
      // and the hard jet spectrum (characterised by a steeply falling exp(-mu pt) shape)
      double expected_jet_qt = sigma*sigma*jet.area()/_mu; // eq.(18)
      // the average qt of the jets used to estimate rho (it should be close to zer).
      double avg_qt = sum_qt/jets_used.size(); // 

      // now we have all the sums, apply the improved_sub
      for (unsigned int in=0; in<_Ns.size(); in++){
        // here we work out the correlation coefficient rN = Cov(qt,QN)/sqrt(Var(qt)Var(Qn))
        info._rNs[in] = (sum_Qq[in]-avg_qt*sum_Qn[in])/
          sqrt((sum_qt2-avg_qt*sum_qt) * (sum_QnQn[in] - pow(sum_Qn[in],2)/jets_used.size()) );

        // make a record of the different contributions
        info._ptshift_term[in] = ffm[in] *_Ns[in]/S1*expected_jet_qt; // first term in eq.(21)
	// 2nd term in eq.(21), evaluated dirrect from rN
        info._correl_term[in]  = info._rNs[in]* sigma*info._sigmaNs[in]*jet.area()/(_mu* pow(S1, _Ns[in]));

        // improved_sub correction
	ffm[in] *= 1 + _Ns[in]/S1*expected_jet_qt; // first term in eq.(21)  
	ffm[in] -= info._correl_term[in];          // second term in eq.(21) 
      }
    } // end of improved_sub correction

    // now that we've done everything needed with the _bge, if relevant reset its density class
    if (_bge) _bge->set_jet_density_class(old_density_class);

    return ffm;
  }

  // initialisation of the internal behaviour
  void JetFFMoments::_initialise(){
    use_scalar_sum      (true);
    set_return_numerator(false);
    set_denominator     (-1.0);  // negative means denominator calculated (and subtracted) automatically
    _mu = -1.0;                 // negative means no improved_sub done
    _jets_for_improved_sub.clear();
  }

  // This function returns the transverse momentum to be used in the denominator
  // of the definition of the moments of a fragmentation function for a jet, e.g.
  // in eq.(3) of arXiv:1209.6086
  // This can be 
  //   1. the (subtracted) jet pt
  //   2. the (subtracted) scalar sum of the pt of the constituents of the jet [this is the default]
  //   3. just the number 1, in case one is interested only in the numerator.
  double JetFFMoments::_compute_normalisation(const PseudoJet & jet, const vector<PseudoJet> &constituents, double &rho, double &sigma) const{
    rho = sigma = 0.0; // in case there is no subtraction;

    if (_return_numerator){ return 1.0;}
    if (_norm>0){ return _norm;}
    if (!_use_scalar_sum){  // use the jet pt
      // if a bge is defined, subtract the appropriate quantity
      if (_bge){
        rho = _bge->rho(jet);
        sigma = _bge->sigma(jet);
        PseudoJet subtracted_jet = jet;
        PseudoJet to_subtract = rho * jet.area_4vector();
        if (to_subtract.pt2() >= jet.pt2()){
          return -1.0;
        }
        return  (jet-to_subtract).pt();
      }

      // no BGE
      return jet.pt();
    }

    // use the scalar pt sum
    double scalar_sum=0.0;
    for (unsigned int i=0; i<constituents.size(); i++){
      scalar_sum += constituents[i].pt();
    }

    // if a bge is defined, subtract the appropriate quantity
    if (_bge){
      const FunctionOfPseudoJet<double> * old_density_class = _bge->jet_density_class();
      BackgroundJetScalarPtDensity scalar_density;
      _bge->set_jet_density_class(&scalar_density);
      rho = _bge->rho(jet);
      sigma = _bge->sigma(jet);
      scalar_sum -= rho * jet.area();
      _bge->set_jet_density_class(old_density_class);
    }
    return scalar_sum;
  }


} // contrib namespace

FASTJET_END_NAMESPACE
