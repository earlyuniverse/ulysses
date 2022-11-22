// $Id$
//
// Copyright (c)-, 
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

#ifndef __FASTJET_CONTRIB_SUBJETCOUNTING_HH__
#define __FASTJET_CONTRIB_SUBJETCOUNTING_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include <vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
//
//   Implementation of the two subjet counting algorithms used in 
//   "Learning How to Count: A High Multiplicity Search for the LHC"
//   Sonia El Hedri, Anson Hook, Martin Jankowiak, Jay G. Wacker
//   JHEP 1308:136,2013
//   http://arxiv.org/abs/1302.1870
//
//   for usage example see example.cc 
//
//   for questions contact jankowiak@gmail.com
//
//   tested with Fastjet 3.0.3
//
//------------------------------------------------------------------------

class SubjetCountingKt : public FunctionOfPseudoJet<unsigned int> {

 // **************************************************************************************
 //
 // constructor for subjet counting with the exclusive Kt algorithm
 //
 // ~ input:  f_Kt (or rather f_Kt*total_jet_pt) defines the scale at which the exclusive Kt algorithm
 //           stops recombination; typical values might be f_Kt ~ 0.04 - 0.10 
 // ~ input:  pt_cut is the minimum pt required of the identified subjets included in n_Kt;
 //           typical values might be 20-70 GeV
 //
 // ~ return: a non-negative integer (accessed through the parent class's () operator as 
 //           stipulated by the FuntionOfPseudoJet paradigm: see example.cc);
 //           alternatively the subjets can be accessed directly via getSubjets
 // 
 //   see original reference for details and discussion
 //
 // **************************************************************************************

public:
  SubjetCountingKt(double f_Kt, double pt_cut) : _f_Kt(f_Kt), _pt_cut(pt_cut) {}

  /// default dtor
  virtual ~SubjetCountingKt() {}

  /// returns the value of the subjet count for a jet's
  /// constituents. (Normally accessed by the parent class's operator()).
  unsigned int result(const PseudoJet& jet) const;

  /// get the actual subjets identified by the algorithm
  std::vector<fastjet::PseudoJet> getSubjets(const fastjet::PseudoJet &jet) const;

  /// returns the description of the class
  std::string description() const; 

private: 
double _f_Kt, _pt_cut; // parameters defining n_Kt

// internal function that computes the result
unsigned int _n_Kt(const fastjet::PseudoJet &jet) const;

}; //end class SubjectCountingKt

class SubjetCountingCA : public FunctionOfPseudoJet<unsigned int> {
 // *************************************************************************************
 //
 // constructor for subjet counting with the CA algorithm
 //
 // ~ input:  mass_cut_off is the mass scale down to which the fat jet is declustered;
 //           typical values might be 20-70 GeV 
 // ~ input:  ycut is the asymmetry parameter that determines how hard to cut on asymmetric
 //           splittings; typical values might be ycut ~ 0.10 - 0.15 
 // ~ input:  R_min defines a lower angular cutoff to how far the declustering proceeds;
 //           this parameter may not make sense or be useful depending on the source of the
 //           constituent objects making up the fat jet (topoclusters, etc.);
 //           can be made inoperative by being set to zero.
 //           typical values might be R_min ~ 0.10 - 0.20 
 // ~ input:  pt_cut is the minimum pt required of the identified subjets included in n_CA;
 //           typical values might be 20-70 GeV
 //
 // ~ return: a non-negative integer (accessed through the parent class's () operator as 
 //           stipulated by the FuntionOfPseudoJet paradigm: see example.cc);
 //           alternatively the subjets can be accessed directly via getSubjets
 //
 //   see original reference for details and discussion
 //
 // *************************************************************************************

public:
  SubjetCountingCA(double mass_cut_off, double ycut, double R_min, double pt_cut) : 
                  _mass_cut_off(mass_cut_off), _ycut(ycut), _R_min(R_min), _pt_cut(pt_cut) {}

  /// default dtor
  virtual ~SubjetCountingCA(){}

  /// returns the value of the subjet count for a jet's
  /// constituents. (Normally accessed by the parent class's operator()).
  unsigned int result(const PseudoJet& jet) const;

  /// get the actual subjets identified by the algorithm
  std::vector<fastjet::PseudoJet> getSubjets(const fastjet::PseudoJet &jet) const;

  /// returns the description of the class
  std::string description() const; 

private: 
double _mass_cut_off, _ycut, _R_min,  _pt_cut; // parameters defining n_CA

// internal function that calculates the result
unsigned int _n_CA(const fastjet::PseudoJet &jet) const;

// internal function called recursively by _n_CA
void _FindHardSubst(const fastjet::PseudoJet & this_jet, 
                   std::vector<fastjet::PseudoJet> & t_parts) const;
}; // end class SubjectCountingCA

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SUBJETCOUNTING_HH__
