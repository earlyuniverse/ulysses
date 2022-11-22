// $Id: FlavorCone.hh 1056 2017-09-06 20:59:06Z philten $
//
// Copyright (c) 2017, Philip Ilten, Nicholas Rodd, Jesse Thaler, 
// and Michael Williams
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

#ifndef __FASTJET_CONTRIB_FLAVORCONE_HH__
#define __FASTJET_CONTRIB_FLAVORCONE_HH__

#include <fastjet/internal/base.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/LimitedWarning.hh>

FASTJET_BEGIN_NAMESPACE // defined in fastjet/internal/base.hh

namespace contrib{

  /** \mainpage FlavorCone contrib
   
   The FlavorCone contrib is a lightweight tool to cluster jets of
   fixed radius around a given number of seeds, where overlapping jets
   are partitioned by nearest neighbor.  The contrib consists of the
   following files:
   - contrib::FlavorConePlugin
   - contrib::FlavorConePlugin::Extras
   - example.cc provides usage examples
   */
  
  
//----------------------------------------------------------------------
/// \class FlavorConePlugin
///
/// Plugin to take input seeds and particles and create FlavorCone
/// jets.
///
/// On construction, one must supply seeds and an rcut value. To obtain
/// the jets call ClusterSequence::inclusive_jets(). The seed for
/// each jet can be accessed via ClusterSequence::extras() with the
/// method FlavorConeExtras::seed(jet), where jet is from the returned
/// jet vector.
class FlavorConePlugin : public JetDefinition::Plugin {
public:
  //--------------------------------------------------------------------
  /// \class Extras
  ///
  /// Hold extras for FlavorConePlugin.
  ///
  /// Access the seed for a given jet constructed with the
  /// FlavorConePlugin.
  class Extras : public ClusterSequence::Extras {
  public:
    /// Main constructor for the FlavorCone Extras class.
    Extras();
    /// Return the seed for a jet. A seed returned with an energy of
    /// -1 indicates the jet has no valid seed.
    const PseudoJet& seed(const PseudoJet &jet) const;
    
  private:
    /// Map which holds the seed for each jet.
    std::map<int, PseudoJet> _seeds;
    /// Stored jet to return when no seed is found.
    const PseudoJet _invalid_seed;
    /// Warning when no seed is found.
    static LimitedWarning _warn_seed;
    friend class FlavorConePlugin;
  };

  /// Main constructor for the FlavorCone Plugin class.
  FlavorConePlugin(const std::vector<PseudoJet> &seeds, double rcut);

  /// Copy constructor.
  FlavorConePlugin(const FlavorConePlugin & plugin);

  // Methods required by base class.
  /// Brief description of FlavorCone for self-documentation (gives
  /// number of seeds and rcut).
  virtual std::string description() const;
  /// Run the FlavorCone clustering.
  virtual void run_clustering(ClusterSequence &cs) const;
  /// Returns rcut value.
  virtual double R() const;
  /// Exclusive sequence is not meaningful for this algorithm.
  virtual bool exclusive_sequence_meaningful() const;
  /// At the moment, FlavorCone only works for pp, not e+e-.
  virtual bool is_spherical() const;

private:
  double                 _rcut;  ///< Jet radius cut.
  std::vector<PseudoJet> _seeds; ///< Vector of seeds.
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_FLAVORCONE_HH__
