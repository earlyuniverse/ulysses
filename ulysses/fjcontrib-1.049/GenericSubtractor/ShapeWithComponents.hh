// $Id: ShapeWithComponents.hh 859 2015-09-21 10:11:32Z gsalam $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of the GenericSubtractor package of FastJet
// Contrib.
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

#ifndef __FASTJET_CONTRIB_SHAPE_WITH_COMPONENTS_HH__
#define __FASTJET_CONTRIB_SHAPE_WITH_COMPONENTS_HH__

#include "fastjet/FunctionOfPseudoJet.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

/// \class ShapeWithComponents
/// base class for jet shapes that provide extra information about
/// about separate components from which they are composed.
/// 
/// Many jet shapes involve ratios or more complicated combinations of
/// different components. E.g. a jet broadening is given by a ratio
///
///     B = (\sum_i p_ti DeltaR_{i,jet}) / (\sum_i p_ti)
///
/// The ShapeWithComponents class provides an interface that allows
/// access not just to the shape but also to the individual
/// components, as well as a mechanism for calculating the components
/// individually and combining them back into the overall shape.  This
/// can be of interest, for example, in pileup subtraction, providing
/// a structure to allow subtraction to be performed automatically on
/// the invidividual components.
///
class ShapeWithComponents : public FunctionOfPseudoJet<double> {
public:
  //------------------------------------------------------------------
  //------- first, a reminder of the standard functions for a FoPJ ---

  /// the result of computing the shape from a given PseudoJet
  virtual double result(const PseudoJet &jet) const{
    return result_from_components(components(jet));
  }

  /// the description that has to be defined in derived classes
  virtual std::string description() const = 0;

  //------------------------------------------------------------------
  //------- below, the functions that provide "hinting" --------------

  /// returns the number of components 
  virtual unsigned int n_components() const = 0;

  /// Returns a vector containing the different components which go
  /// into the calculation of the shape. The vector will be of size
  /// n_components().
  virtual std::vector<double> components(const PseudoJet &jet) const = 0;

  /// Return component i of the shape for the given jet. The numbering
  /// of components is 0...(n_components()-1).
  ///
  /// The default implementation of this function evaluates all
  /// components and then returns the requested one. Derived classes
  /// may wish to proceed differently.
  virtual double component(int i, const PseudoJet &jet) const {
    assert(i < int(n_components()));
    return components(jet)[i];
  }

  /// given a vector of components, determine the result of the event
  /// shape
  virtual double result_from_components(const std::vector <double> &) const = 0;

  /// return a pointer to a new instance of a ComponentShape that
  /// calculates component i of this shape. It is the caller's
  /// responsibility to delete the pointer when it is no longer
  /// needed.
  ///
  /// Note: the generic implementation included below should be adequate
  /// for most uses.
  virtual const FunctionOfPseudoJet<double> * component_shape(unsigned index) const;


  /// a helper class intended to contain a single component of a
  /// composite shape.
  class ComponentShape : public FunctionOfPseudoJet<double> {
  public:
    ComponentShape(const ShapeWithComponents * shape, unsigned index) : 
      _shape(shape), _index(index) {}

    double result(const PseudoJet & jet) const {
      return _shape->component(_index, jet);
    }
  private:
    const ShapeWithComponents * _shape;
    const unsigned _index;
  };


};
  
/// generic implementation of component_shape -- note, this
/// implementation will hide properties of component shapes, e.g. if
/// the underlying shape class derives from ShapeWithPartition, the
/// result of component_shape will not; If that's a problem, you
/// should reimplement component_shape.
inline const FunctionOfPseudoJet<double> * ShapeWithComponents::component_shape(unsigned index) const {
  return new ComponentShape(this, index);
}


} // namespace contrib

FASTJET_END_NAMESPACE


#endif // __FASTJET_CONTRIB_SHAPE_WITH_COMPONENTS_HH__
