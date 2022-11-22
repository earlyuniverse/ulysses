// ConstituentSubtractor package                                                                                                                       
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz                   
                                                                                                                                                      
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
//                                                                                                                                                     // You should have received a copy of the GNU General Public License                                                                                   
// along with this code. If not, see <http://www.gnu.org/licenses/>.                                                                                   
//----------------------------------------------------------------------                                                                               


#include "RescalingClasses.hh"
 
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh                                                                                    


namespace contrib{

  ///-----------------------------------------------------
  /// BackgroundRescalingYPhi
  BackgroundRescalingYPhi::BackgroundRescalingYPhi(double v2, double v3, double v4, double psi, double a1, double sigma1, double a2, double sigma2){
    _v2=v2;
    _v3=v3;
    _v4=v4;
    _psi=psi;
    _a1=a1;
    _sigma1=sigma1;
    _a2=a2;
    _sigma2=sigma2;
    _use_rap=true;
    _use_phi=true;
  }

  void BackgroundRescalingYPhi::use_rap_term(bool use_rap){
    _use_rap=use_rap;
  }

  void BackgroundRescalingYPhi::use_phi_term(bool use_phi){
    _use_phi=use_phi;
  }

  double BackgroundRescalingYPhi::result(const PseudoJet & particle) const {
    double phi_term=1;
    if (_use_phi){
      double phi=particle.phi();
      phi_term=1 + 2*_v2*_v2*cos(2*(phi-_psi)) + 2*_v3*_v3*cos(3*(phi-_psi)) +  2*_v4*_v4*cos(4*(phi-_psi));
    }
    double rap_term=1;
    if (_use_rap){
      double y=particle.rap();
      rap_term=_a1*exp(-y*y/(2*_sigma1*_sigma1)) + _a2*exp(-y*y/(2*_sigma2*_sigma2));
    }

    return phi_term*rap_term;
  }






  ///----------------------------------------------------
  /// BackgroundRescalingYPhiUsingVectorForY
  BackgroundRescalingYPhiUsingVectorForY::BackgroundRescalingYPhiUsingVectorForY(double v2, double v3, double v4, double psi, std::vector<double> values, std::vector<double> rap_binning){
    _v2=v2;
    _v3=v3;
    _v4=v4;
    _psi=psi;
    _values=values;
    _rap_binning=rap_binning;
    _use_phi=true;
    if (_rap_binning.size()>=2){
      _use_rap=true;
      if (_values.size()!=_rap_binning.size()-1) throw Error("BackgroundRescalingYPhiUsingVectorForY (from ConstituentSubtractor) The input vectors have wrong dimension. The vector with binning shuld have the size by one higher than the vector with values.");
    }
    else _use_rap=false;
  }

  void BackgroundRescalingYPhiUsingVectorForY::use_rap_term(bool use_rap){
    _use_rap=use_rap;
    if (use_rap && _rap_binning.size()<2) throw Error("BackgroundRescalingYPhiUsingVectorForY (from ConstituentSubtractor)  Requested rapidity rescaling, but the vector with binning has less than two elements!");
  }

  void BackgroundRescalingYPhiUsingVectorForY::use_phi_term(bool use_phi){
    _use_phi=use_phi;
  }

  double BackgroundRescalingYPhiUsingVectorForY::result(const PseudoJet & particle) const {
    double phi_term=1;
    if (_use_phi){
      double phi=particle.phi();
      phi_term=1 + 2*_v2*_v2*cos(2*(phi-_psi)) + 2*_v3*_v3*cos(3*(phi-_psi)) +  2*_v4*_v4*cos(4*(phi-_psi));
    }
    double rap_term=1;
    if (_use_rap){
      int rap_index=0;
      double rap=particle.rap();
      if (rap<_rap_binning[0]) rap_index=0;  // take the lowest rapidity bin in this case
      else if (rap>=_rap_binning[_rap_binning.size()-1]) rap_index=_rap_binning.size()-2;  // take the highest rapidity bin in this case
      else{
	for (unsigned int i=1;i<_rap_binning.size();++i){
	  if (rap<_rap_binning[i]){
	    rap_index=i-1;
	    break;
	  }
	}
      }
      rap_term=_values[rap_index];
    }

    return phi_term*rap_term;
  }




  ///---------------------------------------------------------
  ///BackgroundRescalingYPhiUsingVectors
  BackgroundRescalingYPhiUsingVectors::BackgroundRescalingYPhiUsingVectors(std::vector<std::vector<double> > values, std::vector<double> rap_binning, std::vector<double> phi_binning){
    _values=values;
    _rap_binning=rap_binning;
    _phi_binning=phi_binning;
    if (_rap_binning.size()>=2) _use_rap=true;
    else _use_rap=false;
    if (_phi_binning.size()>=2) _use_phi=true;
    else _use_phi=false;
  }

  void BackgroundRescalingYPhiUsingVectors::use_rap_term(bool use_rap){
    _use_rap=use_rap;
    if (use_rap && _rap_binning.size()<2) throw Error("BackgroundRescalingYPhiUsingVectors (from ConstituentSubtractor)  Requested rapidity rescaling, but the vector with binning has less than two elements!");
  }

  void BackgroundRescalingYPhiUsingVectors::use_phi_term(bool use_phi){
    _use_phi=use_phi;
    if (use_phi && _phi_binning.size()<2) throw Error("BackgroundRescalingYPhiUsingVectors (from ConstituentSubtractor)  Requested azimuth rescaling, but the vector with binning has less than two elements!");
  }

  double BackgroundRescalingYPhiUsingVectors::result(const PseudoJet & particle) const {
    unsigned int phi_index=0;
    if (_use_phi){
      double phi=particle.phi();
      if (phi<_phi_binning[0] || phi>=_phi_binning[_phi_binning.size()-1]) throw Error("BackgroundRescalingYPhiUsingVectors (from ConstituentSubtractor) The phi binning does not correspond to the phi binning of the particles.");
      for (unsigned int i=1;i<_phi_binning.size();++i){
	if (phi<_phi_binning[i]){
	  phi_index=i-1;
	  break;
	}
      }
    }
    unsigned int rap_index=0;
    if (_use_rap){
      double rap=particle.rap();
      if (rap<_rap_binning[0]) rap_index=0;  // take the lowest rapidity bin in this case
      else if (rap>=_rap_binning[_rap_binning.size()-1]) rap_index=_rap_binning.size()-2;  // take the highest rapidity bin in this case
      else{
	for (unsigned int i=1;i<_rap_binning.size();++i){
	  if (rap<_rap_binning[i]){
	    rap_index=i-1;
	    break;
	  }
	}
      }
    }
    if (_values.size()<=rap_index) throw Error("BackgroundRescalingYPhiUsingVectors (from ConstituentSubtractor) The input vector<vector<double> > with values has wrong size.");
    else if (_values[rap_index].size()<=phi_index) throw Error("BackgroundRescalingYPhiUsingVectors (from ConstituentSubtractor) The input vector<vector<double> > with values has wrong size.");
 
    return _values[rap_index][phi_index];
  }


} // namespace contrib                                                                                                                                  


FASTJET_END_NAMESPACE



