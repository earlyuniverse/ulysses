//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#ifndef __PXCONEPLUGIN_HH__
#define __PXCONEPLUGIN_HH__

#include "fastjet/JetDefinition.hh"
#include "fastjet/internal/thread_safety_helpers.hh"  // helpers to write transparent code w&wo C++11 features

// questionable whether this should be in fastjet namespace or not...

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class PxConePlugin
/// Implementation of the PxCone algorithm (plugin for fastjet v2.1 upwards)
///
/// PxConePlugin is a plugin for fastjet (v2.1 upwards) that provides
/// an interface to the fortran pxcone iterative cone algorithm with
/// midpoint seeds.
///
/// Pxcone was written by Luis del Pozo and Michael H. Seymour. It is
/// not a "supported" program, so if you encounter problems, you are
/// on your own...
///
/// Note that pxcone sometimes encounters non-stable iterations; in
/// such cases it returns an error -- the plugin propagates this by
/// throwing a fastjet::Error exception; if the user wishes to have
/// robust code, they should catch this exception.
///
/// Pxcone has a hard-coded limit (by default 4000) on the maximum
/// number of particles and protojets; if the number of particles or
/// protojets exceeds this, again a fastjet::Error exception will be
/// thrown.
///
/// The functionality of pxcone is described at 
/// http://www.hep.man.ac.uk/u/wplano/ConeJet.ps
//
//----------------------------------------------------------------------
class PxConePlugin : public JetDefinition::Plugin {
public:

  /// constructor for the PxConePlugin, whose arguments have the
  /// following meaning:
  ///
  ///   - the cone_radius is as usual in cone algorithms
  ///
  ///   - stables cones (protojets) below min_jet_energy are discarded
  ///     before calling the splitting procedure to resolve overlaps
  ///     (called epslon in pxcone).
  ///
  ///   - when two protojets overlap, if
  ///       (overlapping_Et)/(Et_of_softer_protojet) < overlap_threshold
  ///     the overlapping energy is split between the two protojets;
  ///     otherwise the less energetic protojet is discarded. Called
  ///     ovlim in pxcone.
  ///
  ///   - pxcone carries out p-scheme recombination, and the resulting 
  ///     jets are massless; setting E_scheme_jets = true (default
  ///     false) doesn't change the jet composition, but the final
  ///     momentum sum for the jets is carried out by direct
  ///     four-vector addition instead of p-scheme recombination.
  ///
  ///   - mode: set to 1 for the e+e- version
  ///           set to 2 for the hadron-hadron version (the default)
  ///
  PxConePlugin (double  cone_radius_in, 
		double  min_jet_energy_in = 5.0, 
		double  overlap_threshold_in = 0.5,
                bool    E_scheme_jets_in = false,
                int     mode = 2) : 
    _cone_radius        (cone_radius_in      ),
    _min_jet_energy     (min_jet_energy_in   ),
    _overlap_threshold  (overlap_threshold_in),
    _E_scheme_jets      (E_scheme_jets_in    ),
    _mode               (mode                ){}


  // some functions to return info about parameters ----------------

  /// the cone radius
  double cone_radius        () const {return _cone_radius        ;}

  /// minimum jet energy (protojets below this are thrown own before
  /// merging/splitting) -- called epslon in pxcone
  double min_jet_energy     () const {return _min_jet_energy     ;}

  /// Maximum fraction of overlap energy in a jet -- called ovlim in pxcone.
  double overlap_threshold  () const {return _overlap_threshold  ;}

  /// if true then the final jets are returned as the E-scheme recombination
  /// of the particle momenta (by default, pxcone returns massless jets with
  /// a mean phi,eta type of recombination); regardless of what is
  /// returned, the internal pxcone jet-finding procedure is
  /// unaffected.
  bool E_scheme_jets()         const {return _E_scheme_jets      ;}

  int mode()                   const {return _mode               ;}

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return cone_radius();}

private:

  double _cone_radius       ;
  double _min_jet_energy    ;
  double _overlap_threshold ;

  bool _E_scheme_jets;

  static thread_safety_helpers::FirstTimeTrue _first_time;
  int _mode;  // 1 = e+e-, 2 = hh (default)

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __PXCONEPLUGIN_HH__
