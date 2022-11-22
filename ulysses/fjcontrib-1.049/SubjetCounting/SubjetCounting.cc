// $Id$
//
// Copyright (c) -, 
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

#include <sstream>
#include "SubjetCounting.hh"

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
//------------------------------------------------------------------------

 // **************************************************************************************
 //
 // subjet counting with the exclusive Kt algorithm
 //
 // **************************************************************************************

 // the actual n_Kt function; called internally
 unsigned int SubjetCountingKt::_n_Kt(const fastjet::PseudoJet &jet) const
 {
 return (unsigned int)SubjetCountingKt::getSubjets(jet).size();
 }

 // get the subjets identified by the Kt subjet counting algorithm
 std::vector<fastjet::PseudoJet> SubjetCountingKt::getSubjets(const fastjet::PseudoJet &jet) const
 {
 JetAlgorithm algorithm = kt_algorithm;
 // this choice of jet_rad ensures that no beam jets are identified
 double jet_rad = fastjet::JetDefinition::max_allowable_R; 
 JetDefinition jetDef = JetDefinition(algorithm, jet_rad, E_scheme, Best);
 ClusterSequence clust_seq(jet.constituents(), jetDef);
 double total_jet_pt=jet.perp();
 double Kt_scale = total_jet_pt*total_jet_pt*_f_Kt*_f_Kt;
 // Kt_scale has to be scaled down to compensate for the choice of jet_rad
 Kt_scale /= fastjet::JetDefinition::max_allowable_R*fastjet::JetDefinition::max_allowable_R;
 std::vector<PseudoJet> Kt_jets = sorted_by_pt(clust_seq.exclusive_jets(Kt_scale));
 std::vector<PseudoJet> return_subjets;

 for (int k=0; k<(int)Kt_jets.size(); k++) 
   if (Kt_jets[k].perp()>_pt_cut) 
     return_subjets.push_back(Kt_jets[k]);

 return return_subjets;
 }

 // *************************************************************************************
 //
 // subjet counting with the CA algorithm
 //
 // *************************************************************************************

 // internal function called recursively by the CA variant of getSubjets (see below)
 void SubjetCountingCA::_FindHardSubst(const fastjet::PseudoJet & this_jet, 
                                     std::vector<fastjet::PseudoJet> & t_parts) const
 {
  fastjet::PseudoJet parent1(0.0,0.0,0.0,0.0), parent2(0.0,0.0,0.0,0.0);
  bool had_parents=this_jet.validated_cs()->has_parents(this_jet,parent1,parent2);

  if ( (this_jet.m() < _mass_cut_off) || (!had_parents) )
  {
  t_parts.push_back(this_jet);  return;
  } // stop recursion on this branch

  if (had_parents && parent1.plain_distance(parent2) < (_R_min*_R_min) )
  {
  t_parts.push_back(this_jet);  return;
  } // stop recursion on this branch

  if (parent1.perp() < parent2.perp()) std::swap(parent1,parent2);

  double pt1=parent1.perp();
  double pt2=parent2.perp();
  double totalpt=pt1+pt2;

  if (pt2>_ycut*totalpt)   
  {
    _FindHardSubst(parent1, t_parts);
    _FindHardSubst(parent2, t_parts);
  } // continue recursion on both branches

  else
    _FindHardSubst(parent1, t_parts);
    // continue recursion on harder branch, discarding softer branch

  return;
 }

 // the actual n_CA function, which is used internally; 
 // most of the hard work done by _FindHardSubst which is called by getSubjets 
 unsigned int SubjetCountingCA::_n_CA(const fastjet::PseudoJet &jet) const
 {
 return (unsigned int)SubjetCountingCA::getSubjets(jet).size();
 }

 // get the subjets identified by the CA subjet counting algorithm
 std::vector<fastjet::PseudoJet> SubjetCountingCA::getSubjets(const PseudoJet& jet) const
 {
 JetAlgorithm algorithm = cambridge_algorithm;
 double jet_rad = fastjet::JetDefinition::max_allowable_R; 
 // want to make sure we capture all the constituents in a single jet
 JetDefinition jetDef = JetDefinition(algorithm, jet_rad, E_scheme, Best);
 ClusterSequence clust_seq(jet.constituents(), jetDef);
 std::vector<PseudoJet> ca_jets = sorted_by_pt(clust_seq.inclusive_jets());
 
 std::vector<fastjet::PseudoJet> empty_parts, return_subjets;
 _FindHardSubst(ca_jets[0], empty_parts);

 for (int k=0; k<(int)empty_parts.size(); k++) 
   if (empty_parts[k].perp() > _pt_cut) 
     return_subjets.push_back(empty_parts[k]);

 return return_subjets;
 }

 //---------------------------------
 unsigned int SubjetCountingKt::result(const PseudoJet& jet) const {
   
   // if jet does not have constituents, throw error
   if (!jet.has_constituents()) throw Error("SubjetCountingKt called on jet with no constituents.");
   
   return _n_Kt(jet);
 }

 unsigned int SubjetCountingCA::result(const PseudoJet& jet) const {
   
   // if jet does not have constituents, throw error
   if (!jet.has_constituents()) throw Error("SubjetCountingCA called on jet with no constituents.");
   
   return _n_CA(jet);
 }

 //---------------------------------
 std::string SubjetCountingCA::description() const {
  std::ostringstream oss;
  oss << "SubjetCountingCA using ";
    oss << "parameters mass_cutoff = " << _mass_cut_off << 
                ", ycut = " << _ycut << ", Rmin = " << _R_min << " and pt_cut = " << _pt_cut; 
  return oss.str();
 }

 std::string SubjetCountingKt::description() const {
  std::ostringstream oss;
  oss << "SubjetCountingKt using ";
  oss << "parameters f_Kt = " << _f_Kt << " and pt_cut = " << _pt_cut;
  return oss.str();
 }
 
} // namespace contrib

FASTJET_END_NAMESPACE
