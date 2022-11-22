// $Id$
//
// Copyright (c) 2013, Oxford University.
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib ScJet.
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

#ifndef __FASTJET_CONTRIB_SCJET_HH__
#define __FASTJET_CONTRIB_SCJET_HH__

#include <fastjet/internal/base.hh>
#include <fastjet/JetDefinition.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

class ClusterSequence;

namespace contrib{

//------------------------------------------------------------------------
/// \class ScJet
/// Jet algorithm based on semi-classical jet model.
///
/// Implementation of an NNH-based sequential recombination algorithm
/// using distance measures derived from a semi-classical jet model.
/// Documented in J Tseng, H Evans,
/// "Semi-classical approach to sequential recombination algorithms
/// for jet clustering", arXiv:1304.1025 (2013).
class ScJet : public JetDefinition::Plugin {
public:

  /// enumeration to specify types of energy input.
  /// mt and pt are longitudinal boost-invariant, et is not.
  typedef enum { use_mt, use_pt, use_et } energyModeType;

  /// print energy mode
  static std::string energyModeString(energyModeType mode) {
    switch (mode) {
      case use_mt : return std::string("Mt");
      case use_pt : return std::string("Pt");
      case use_et : return std::string("Et");
      default :     return std::string("Undefined");
    }
  }

  /// default ctor
  ScJet(double Rin, energyModeType mode = use_mt, int Rexp = 3) :
    _R(Rin), _Rexp(Rexp), _energyMode(mode) {}

  /// copy constructor
  ScJet(const ScJet& plugin) :
    _R(plugin._R), _Rexp(plugin._Rexp), _energyMode(plugin._energyMode) {}

  /// default dtor
  ~ScJet(){}

  /// required by base class
  virtual std::string description() const;

  /// required by base class
  virtual void run_clustering(ClusterSequence&) const;

  /// plugin mechanism's standard way of accessing the jet radius
  virtual double R() const { return _R; }

  /// return radius exponent parameter
  int Rexp() const { return _Rexp; }

  /// return energy mode
  energyModeType energyMode() const { return _energyMode; }

  /// return energy mode string (non-static version)
  std::string energyModeString() const { return energyModeString(_energyMode); }

  /// avoid warning when the user requests exclusive jets
  virtual bool exclusive_sequence_meaningful() const { return true; }

private:

  double _R; // jet radius parameter
  int _Rexp; // radius exponent parameter
  energyModeType _energyMode; // energy calculation type

};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SCJET_HH__
