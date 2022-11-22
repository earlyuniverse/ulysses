// ConstituentSubtractor package                                                                                                                        
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz                   
//                                                                                                                                                     
// Copyright (c) 2014-, Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner                                                                   
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

#ifndef __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__
#define __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__



#include <fastjet/FunctionOfPseudoJet.hh>
#include <iostream>
#include <vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh                                                                                     

namespace contrib{


  template<class T>
  class BackgroundRescalingYFromRoot : public FunctionOfPseudoJet<double> {
  public:
    ///
    /// construct a background rescaling function using ROOT TH1 histogram bin contents (rapidity binning)
    BackgroundRescalingYFromRoot(): _hist(0) {}
    BackgroundRescalingYFromRoot(T* hist=0) {_hist = hist;}

    // return the rescaling factor associated with this jet  
    virtual double result(const PseudoJet & particle) const {
      if (!_hist){
	throw Error("BackgroundRescalingYFromRoot (from ConstituentSubtractor)  The histogram for rescaling not defined! ");
      }
      double rap = particle.rap();
      if (rap<_hist->GetXaxis()->GetBinLowEdge(1)) return _hist->GetBinContent(1);
      if (rap>=_hist->GetXaxis()->GetBinUpEdge(_hist->GetNbinsX())) return _hist->GetBinContent(_hist->GetNbinsX());
      int bin=_hist->FindBin(rap);
      return _hist->GetBinContent(bin);
    }

  private:
    T* _hist;
  };



  template<class T>
  class BackgroundRescalingYPhiFromRoot : public FunctionOfPseudoJet<double> {
  public:
    ///
    /// construct a background rescaling function using ROOT TH2 histogram bin contents (rapidity vs azimuth binning)
    BackgroundRescalingYPhiFromRoot(): _hist(0) {}
    BackgroundRescalingYPhiFromRoot(T* hist=0) {_hist = hist;}

    // return the rescaling factor associated with this jet  
    virtual double result(const PseudoJet & particle) const {
      if (!_hist){
	throw Error("BackgroundRescalingYPhiFromRoot (from ConstituentSubtractor)  The histogram for rescaling not defined! ");
      }
      double rap = particle.rap();
      double phi = particle.phi();
      int xbin=1;
      if (rap<_hist->GetXaxis()->GetBinLowEdge(1)) xbin=1;
      else if (rap>=_hist->GetXaxis()->GetBinUpEdge(_hist->GetNbinsX())) xbin=_hist->GetNbinsX();
      else xbin=_hist->GetXaxis()->FindBin(rap);
      int ybin=1;
      if (phi<_hist->GetYaxis()->GetBinLowEdge(1) || phi>_hist->GetYaxis()->GetBinUpEdge(_hist->GetNbinsY())){
	throw Error("BackgroundRescalingYPhiFromRoot (from ConstituentSubtractor)  The phi range of the histogram does not correspond to the phi range of the particles! Change the phi range of the histogram.");
      }
      else ybin=_hist->GetYaxis()->FindBin(phi);
      return _hist->GetBinContent(xbin,ybin);
    }

  private:
    T* _hist;
  };





  template<class T>
  class BackgroundRescalingYFromRootPhi : public FunctionOfPseudoJet<double> {
  public:
    ///  Construct background rescaling function in rapidity and azimuth using ROOT TH1 histogram bin contents for the rapidity dependence and this parametrization for the azimuth:

    ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))
    ///  with four parameters v2, v3, v4, and psi.

    ///  This product of the TH1 histogram and function is used to rescale the background which is subtracted such that one can correctly account
    ///  for the modulation of the UE due to rapidity dependence of the particle production
    ///  and/or due to the modulation in the azimuthal angle which is characteristic for heavy ion collisions.
    ///  The overall normalization of the rescaling function is arbitrary since it divides out in the calculation of position dependent rho (background is first demodulated to obtain unbiased position independent rho, and then it is modulated to obtain position dependent rho, see fastjet classes GridMedianBackgroundEstimator and JetMedianBackgroundEstimator for detailed calculation).

    BackgroundRescalingYFromRootPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _use_rap(false), _use_phi(false), _hist(0) {}
    BackgroundRescalingYFromRootPhi(double v2, double v3, double v4, double psi, T* hist=0){
      _v2=v2;
      _v3=v3;
      _v4=v4;
      _psi=psi;
      _hist=hist;
      _use_phi=true;
      if (!_hist){
	std::cout << std::endl << std::endl << "ConstituentSubtractor::BackgroundRescalingYFromRootPhi WARNING: The histogram for rapidity rescaling is not defined!!! Not performing rapidity rescaling." << std::endl << std::endl << std::endl;
	_use_rap=false;
      }
      else _use_rap=true;
    }

    ///
    /// Turn on or off the rapidity rescaling. Throwing in case true is set and no histogram is provided.
    void use_rap_term(bool use_rap){
      _use_rap=use_rap;
      if (!_hist && _use_rap){
	throw Error("BackgroundRescalingYFromRootPhi (from ConstituentSubtractor)  Requested rapidity rescaling, but the histogram for rescaling is not defined!");
      }
    }

    ///
    /// Turn on or off the azimuth rescaling.
    void use_phi_term(bool use_phi){
      _use_phi=use_phi;
    }
    
    ///
    /// Return the rescaling factor associated with this particle
    virtual double result(const PseudoJet & particle) const{
      double phi_term=1;
      if (_use_phi){
	double phi=particle.phi();
	phi_term=1 + 2*_v2*_v2*cos(2*(phi-_psi)) + 2*_v3*_v3*cos(3*(phi-_psi)) +  2*_v4*_v4*cos(4*(phi-_psi));
      }
      double rap_term=1;
      if (_use_rap){
	double y=particle.rap();
	int bin=_hist->FindBin(y);
	rap_term=_hist->GetBinContent(bin);
      }

      return phi_term*rap_term;
    }

  private:
    double _v2, _v3, _v4, _psi;
    bool _use_rap, _use_phi;
    T* _hist;
  };
  





  class BackgroundRescalingYPhi : public FunctionOfPseudoJet<double> {
  public:
    ///
    ///  Construct background rescaling function in rapidity and azimuth using this parameterization:
    ///    
    ///  f(y,phi) = phi_term(phi) * rap_term(y)
    ///  where
    ///
    ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))
    ///  with four parameters v2, v3, v4, and psi.
    ///
    ///  rap_term(y) = a1*exp(-pow(y,2)/(2*sigma1^2)) + a2*exp(-pow(y,2)/(2*sigma2^2))
    ///  with four parameters sigma1, sigma2, a1, and a2. 
    ///
    ///  This function is used to rescale the background which is subtracted such that one can correctly account
    ///  for the modulation of the UE due to rapidity dependence of the particle production
    ///  and/or due to the modulation in the azimuthal angle which is characteristic for heavy ion collisions.
    ///  The overall normalization of function f is arbitrary since it divides out in the calculation of position dependent rho (background is first demodulated to obtain unbiased position independent rho, and then it is modulated to obtain position dependent rho, see fastjet classes GridMedianBackgroundEstimator and JetMedianBackgroundEstimator for detailed calculation).
    
    BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
    BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2);
    
    void use_rap_term(bool use_rap);
    void use_phi_term(bool use_phi);
    
    /// return the rescaling factor associated with this jet  
    virtual double result(const PseudoJet & particle) const;
  private:
    double _v2, _v3, _v4, _psi, _a1, _sigma1, _a2, _sigma2;
    bool _use_rap, _use_phi;
  };

  

  class BackgroundRescalingYPhiUsingVectorForY : public FunctionOfPseudoJet<double> {
  public:
    ///
    ///  Construct background rescaling function in rapidity and azimuth using this parameterization:
    ///    
    ///  f(y,phi) = phi_term(phi) * rap_term(y)
    ///  where
    ///
    ///  phi_term(phi) = 1 + 2 * v2^2 * cos(2*(phi-psi)) + 2 * v3^2 * cos(3*(phi-psi)) +  2 * v4^2 * cos(4*(phi-psi))
    ///  with four parameters v2, v3, v4, and psi.
    ///
    ///  rap_term(y) = provided in a vector. 
    ///
    /// The size of the input vector "values" for rapidity dependence is N bins and the corresponding binning should be specified in a separate input vector "rap_binning" of size (N+1). The bin boundaries must be in increasing order.
    ///
    ///  This function is used to rescale the background which is subtracted such that one can correctly account
    ///  for the modulation of the UE due to rapidity dependence of the particle production
    ///  and/or due to the modulation in the azimuthal angle which is characteristic for heavy ion collisions.
    ///  The overall normalization of function f is arbitrary since it divides out in the calculation of position dependent rho (background is first demodulated to obtain unbiased position independent rho, and then it is modulated to obtain position dependent rho, see fastjet classes GridMedianBackgroundEstimator and JetMedianBackgroundEstimator for detailed calculation).
    
    BackgroundRescalingYPhiUsingVectorForY(): _v2(0), _v3(0), _v4(0), _psi(0), _values(0), _rap_binning(0), _use_rap(false), _use_phi(false) {}
    BackgroundRescalingYPhiUsingVectorForY(double v2, double v3, double v4, double psi, std::vector<double> values, std::vector<double> rap_binning);
    
    void use_rap_term(bool use_rap);
    void use_phi_term(bool use_phi);
    
    /// return the rescaling factor associated with this jet  
    virtual double result(const PseudoJet & particle) const;
  private:
    double _v2, _v3, _v4, _psi;
    std::vector<double> _values;
    std::vector<double> _rap_binning;
    bool _use_rap, _use_phi;
  };

  

  class BackgroundRescalingYPhiUsingVectors : public FunctionOfPseudoJet<double> {
  public:
    ///
    ///  Construct background rescaling function in rapidity and azimuth using dependence recorded in input object vector<vector<double> > values. Its size is N bins for rapidity and M bins for azimuth. The binning of the rapidity should be specified in a separate vector "rap_binning" of size (N+1), and similarly the binning of the azimuth  should be specified in a separate vector "phi_binning" of size (M+1). The bin boundaries must be in increasing order.  
    ///

    BackgroundRescalingYPhiUsingVectors(): _values(0), _rap_binning(0), _phi_binning(0) {}
    BackgroundRescalingYPhiUsingVectors(std::vector<std::vector<double> > values, std::vector<double> rap_binning, std::vector<double> phi_binning);
    
    void use_rap_term(bool use_rap);
    void use_phi_term(bool use_phi);
    
    /// return the rescaling factor associate
    virtual double result(const PseudoJet & particle) const;
  private:
    std::vector<std::vector<double> > _values;
    std::vector<double> _rap_binning;
    std::vector<double> _phi_binning;
    bool _use_rap, _use_phi;
  };

  

} // namespace contrib                                                                                                                                  

FASTJET_END_NAMESPACE


#endif   //__FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_RESCALINGCLASSES_HH__
