// $Id: LundGenerator.hh 1289 2021-11-09 11:53:53Z scyboz $
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

#ifndef __FASTJET_CONTRIB_LUNDGENERATOR_HH__
#define __FASTJET_CONTRIB_LUNDGENERATOR_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/tools/Recluster.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include <string>
#include <vector>
#include <utility>

// TODO:
// - add interface to write declusterings to json files
//   [gps, possibly as a separate header, in order to factorise the json.hh dependence]
// - something for pileup subtraction?
// - do we want to update json.hh to latest? And handle
//   the precision issue more elegantly than the current 
//   hack of editing json.hh
// - what do we do about the fact that json.hh is c++11?

FASTJET_BEGIN_NAMESPACE

namespace contrib{

class LundGenerator;

//----------------------------------------------------------------------
/// \class LundDeclustering
/// Contains the declustering variables associated with a single qnode
/// on the Lund plane
class LundDeclustering {
public:

  /// return the pair PseudoJet, i.e. sum of the two subjets
  const PseudoJet & pair()  const {return pair_;}
  /// returns the subjet with larger transverse momentum
  const PseudoJet & harder() const {return harder_;}
  /// returns the subjet with smaller transverse momentum
  const PseudoJet & softer() const {return softer_;}


  /// returns pair().m() [cached]
  double m()         const {return m_;}

  /// returns the rapidity-azimuth separation of the pair of subjets [cached]
  double Delta()     const {return Delta_;}

  /// returns softer().pt() / (softer().pt() + harder().pt()) [cached]
  double z()         const {return z_;}

  /// returns softer().pt() * Delta() [cached]
  double kt()        const {return kt_;}

  /// returns z() * Delta() [cached]
  double kappa()     const {return kappa_;}
  
  /// returns an azimuthal type angle of softer() around harder()
  double psi()       const {return psi_;}
  
  /// returns the x,y coordinates that are used in the Lund-plane plots
  /// of arXiv:1807.04758: ln(1/Delta()), and ln(kt()) respectively
  std::pair<double,double> const lund_coordinates() const {
    return std::pair<double,double>(std::log(1.0/Delta()),std::log(kt()));
  }

  virtual ~LundDeclustering() {}

private:
  double m_, Delta_, z_, kt_, kappa_, psi_;
  PseudoJet pair_, harder_, softer_;

protected:
  /// the constructor is private, because users will not generally be
  /// constructing a LundDeclustering element themselves.
  LundDeclustering(const PseudoJet& pair,
		   const PseudoJet& j1, const PseudoJet& j2);

  friend class LundGenerator;

};
  

//----------------------------------------------------------------------
/// \class LundGenerator
/// Generates vector of LundDeclustering for a given jet
/// corresponding to its Lund plane.
class LundGenerator : public FunctionOfPseudoJet< std::vector<LundDeclustering> > {
public:
  /// LundGenerator constructor
  LundGenerator(JetAlgorithm jet_alg = cambridge_algorithm)
    : recluster_(JetDefinition(jet_alg, JetDefinition::max_allowable_R)){}
    
  /// LundGenerator constructor with jet definition
  LundGenerator(const JetDefinition & jet_def) : recluster_(jet_def){}
  
  /// destructor
  virtual ~LundGenerator() {}

  /// obtain the declusterings of the primary plane of the jet
  virtual std::vector<LundDeclustering> result(const PseudoJet& jet) const;

  /// description of the class
  virtual std::string description() const;
  
private:
  /// recluster definition
  Recluster recluster_;
};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_LUNDGENERATOR_HH__

