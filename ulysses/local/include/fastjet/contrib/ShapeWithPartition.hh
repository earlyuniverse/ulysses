// $Id: ShapeWithPartition.hh 859 2015-09-21 10:11:32Z gsalam $
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

#ifndef __FASTJET_CONTRIB_SHAPE_WITH_PARTITION_HH__
#define __FASTJET_CONTRIB_SHAPE_WITH_PARTITION_HH__

#include "fastjet/FunctionOfPseudoJet.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//------------------------------------------------------------------------
/// \class ShapeWithPartition
///
/// This class describes a shape where the partitioning has been
/// disconnected from the result of the shape (obtained from that
/// pre-computed partition)
///
/// This is meant to be a helper for operations in which we want to
/// recompute the result several times while keeping the partition
/// fixed (such as GenericSubtraction for which this would avoid any
/// back-reaction effect)
class ShapeWithPartition : public FunctionOfPseudoJet<double>{
public:
  /// just a reminder that one has to define this too in derived classes
  virtual std::string description() const = 0;

  /// compute the partition associated with a given jet.
  /// It is returned under the form of a composite jet.
  virtual PseudoJet partition(const PseudoJet &jet) const = 0;

  /// return the result of the shape computed from a given partition
  virtual double result_from_partition(const PseudoJet &partit) const = 0;

  /// the overall result of the shape. This is nothing that the
  /// "result_from_partition" obtained from the partition of the given
  /// jet.
  virtual double result(const PseudoJet &jet) const{
    return result_from_partition(partition(jet));
  }
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif // __FASTJET_CONTRIB_SHAPE_WITH_PARTITION_HH__
