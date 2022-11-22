// $Id: ExampleShapes.hh 859 2015-09-21 10:11:32Z gsalam $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_CONTRIB_EXAMPLE_SHAPES_HH__
#define __FASTJET_CONTRIB_EXAMPLE_SHAPES_HH__

#include "fastjet/FunctionOfPseudoJet.hh"
#include "ShapeWithPartition.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//------------------------------------------------------------------------
// a series of simple shapes:
//  - Pt
//  - ScalarPt
//  - Mt
//  - Mass
//  - MassSquare
//------------------------------------------------------------------------
/// the jet transverse momentum
class Pt : public FunctionOfPseudoJet<double>{
public:
  virtual std::string description() const{return "jet pt";}
  virtual double result(const PseudoJet &jet) const{ return jet.pt();}
};

/// the jet scalar transverse momentum
class ScalarPt : public FunctionOfPseudoJet<double>{
public:
  virtual std::string description() const{return "jet scalar pt";}
  virtual double result(const PseudoJet &jet) const;
};

/// the transverse mass 
class Mt : public FunctionOfPseudoJet<double>{
public:
  virtual std::string description() const{return "jet transverse mass";}
  virtual double result(const PseudoJet &jet) const{ return jet.mt();}
};

/// the jet mass
class Mass : public FunctionOfPseudoJet<double>{
public:
  virtual std::string description() const{return "jet mass";}
  virtual double result(const PseudoJet &jet) const{ return jet.m();}
};

/// the jet mass squared
class MassSquare : public FunctionOfPseudoJet<double>{
public:
  virtual std::string description() const{return "jet mass squared";}
  virtual double result(const PseudoJet &jet) const{ return jet.m2();}
};


//------------------------------------------------------------------------
/// the of the kt distance between the 2 kt subjets
class KtDij : public ShapeWithPartition{
public:
  virtual std::string description() const{return "kt distance between the 2 kt subjets";}
  virtual PseudoJet partition(const PseudoJet &jet) const;
  virtual double result_from_partition(const PseudoJet &partit) const;
};

//------------------------------------------------------------------------
/// angularities
///
/// This is defined as 
/// \f[
///   \theta^{(\alpha)} = 1/{\tilde p_t} \sum_i p_{ti} \Delta R_{i,jet}^{2-\alpha}
/// \f]
class Angularity : public FunctionOfPseudoJet<double>{
public:
  /// default ctor
  Angularity(double alpha=1.0) : _alpha(alpha){}

  /// description
  virtual std::string description() const;

  /// compute the function
  virtual double result(const PseudoJet &jet) const;

protected:
  double _alpha;
};


/// angularities: numerator
///
/// Angularities are defined as
/// \f[
///   \theta^{(\alpha)} = 1/pttilde \sum_i p_{ti} \Delta R_{i,jet}^{2-\alpha}
/// \f]
/// this calss computes just the numerator
class AngularityNumerator : public FunctionOfPseudoJet<double>{
public:
  /// default ctor
  AngularityNumerator(double alpha=1.0) : _alpha(alpha){}

  /// description
  virtual std::string description() const;

  /// compute the function
  virtual double result(const PseudoJet &jet) const;

protected:
  double _alpha;
};


/// Energy-energy correlators
///
/// Defined as 
/// \f[
///   \tau_EEC = \sum_{i,j} p_ti p_tj \theta_{ij}^\beta / p_t^2
/// \f]
class TauEEC : public FunctionOfPseudoJet<double>{
public:
  /// default ctor
  TauEEC(double i_beta=1.0) : _beta(i_beta){}

  /// description
  virtual std::string description() const;

  double beta() const {return _beta;}

  /// compute the function
  virtual double result(const PseudoJet &jet) const;

protected:
  double _beta;
};


//------------------------------------------------------------------------
/// Energy-energy correlators: numerator
///
/// Defined as 
/// \f[
///   \tau_EEC = \sum_{i,j} p_ti p_tj \theta_{ij}^\beta / p_t^2
/// \f]
class TauEECNumerator : public ShapeWithPartition{
public:
  /// default ctor
  TauEECNumerator(double beta=1.0) : _beta(beta){}

  /// description
  virtual std::string description() const;

  /// compute the function
  virtual double result(const PseudoJet &jet) const;

protected:
  double _beta;
};

// a quick coding of the "numerator" of N-subjettiness
class NSubjettinessNumerator : public ShapeWithPartition{
public:
  NSubjettinessNumerator(unsigned int N) : _N(N){}
  virtual std::string description() const{ return "N-subjettiness numerator";}

  /// compute the partition associated with a given jet
  virtual PseudoJet partition(const PseudoJet &jet) const;

  /// compute the shape from the partition
  virtual double result_from_partition(const PseudoJet &partit) const;

protected:
  const unsigned int _N;
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_EXAMPLE_SHAPES_HH__
