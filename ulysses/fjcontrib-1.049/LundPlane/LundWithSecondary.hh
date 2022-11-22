// $Id: LundWithSecondary.hh 1289 2021-11-09 11:53:53Z scyboz $
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

#ifndef __FASTJET_CONTRIB_LUNDWITHSECONDARY_HH__
#define __FASTJET_CONTRIB_LUNDWITHSECONDARY_HH__

#include "LundGenerator.hh"
#include "SecondaryLund.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//----------------------------------------------------------------------
/// \class LundWithSecondary
/// Define a primary and secondary Lund plane
///
///  \param secondary_def definition used for the leading emission.
//
class LundWithSecondary {
public:
  /// LundWithSecondary constructor
  LundWithSecondary(SecondaryLund * secondary_def = 0)
    : secondary_def_(secondary_def) {}
  
  /// LundWithSecondary constructor with jet alg
  LundWithSecondary(JetAlgorithm jet_alg,
		    SecondaryLund * secondary_def = 0)
    : lund_gen_(jet_alg), secondary_def_(secondary_def) {}
  
  /// LundWithSecondary constructor with jet def
  LundWithSecondary(const JetDefinition & jet_def,
		    SecondaryLund * secondary_def = 0)
    : lund_gen_(jet_def), secondary_def_(secondary_def) {}
  
  /// destructor
  virtual ~LundWithSecondary() {}
  
  /// primary Lund declustering
  std::vector<LundDeclustering> primary(const PseudoJet& jet) const;
  
  /// secondary Lund declustering (slow)
  std::vector<LundDeclustering> secondary(const PseudoJet& jet) const;
  
  /// secondary Lund declustering with primary sequence as input
  std::vector<LundDeclustering> secondary(
			 const std::vector<LundDeclustering> & declusts) const;

  /// return the index associated of the primary declustering that is to be
  /// used for the secondary plane.
  int secondary_index(const std::vector<LundDeclustering> & declusts) const;
  
  /// description of the class
  std::string description() const;
  
private:
  /// lund generator
  LundGenerator lund_gen_;
  
  /// secondary definition
  SecondaryLund * secondary_def_;
};

  
// //----------------------------------------------------------------------
// /// \class LundWithSecondaryAndTertiary
// /// Define a primary, secondary and tertiary Lund plane
// class LundWithSecondaryAndTertiary : public LundWithSecondary {
// public:
//   /// LundWithSecondaryAndTertiary constructor
//   LundWithSecondaryAndTertiary(SecondaryLund * secondary_def = 0,
// 			       SecondaryLund * tertiary_def = 0)
//     : LundWithSecondary(secondary_def), tertiary_def_(tertiary_def) {}
  
//   /// LundWithSecondaryAndTertiary constructor with jet alg
//   LundWithSecondaryAndTertiary(JetAlgorithm jet_alg,
// 		    SecondaryLund * secondary_def = 0,
// 		    SecondaryLund * tertiary_def = 0)
//     : LundWithSecondary(jet_alg, secondary_def), tertiary_def_(tertiary_def) {}
  
//   /// LundWithSecondaryAndTertiary constructor with jet def
//   LundWithSecondaryAndTertiary(const JetDefinition & jet_def,
// 			       SecondaryLund * secondary_def = 0,
// 			       SecondaryLund * tertiary_def = 0)
//     : LundWithSecondary(jet_def, secondary_def), tertiary_def_(tertiary_def) {}
  
//   /// destructor
//   virtual ~LundWithSecondaryAndTertiary() {}
  
//   /// tertiary Lund declustering
//   virtual std::vector<LundDeclustering> tertiary(const PseudoJet& jet) const;

//   /// description of the class
//   virtual std::string description() const;
  
// private:
//   /// tertiary definition
//   SecondaryLund * tertiary_def_;
// };


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_LUNDWITHSECONDARY_HH__

