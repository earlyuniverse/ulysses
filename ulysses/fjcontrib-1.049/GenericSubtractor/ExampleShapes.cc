// $Id: ExampleShapes.cc 859 2015-09-21 10:11:32Z gsalam $
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

#include "fastjet/ClusterSequenceArea.hh"
#include "ExampleShapes.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//------------------------------------------------------------------------
double ScalarPt::result(const PseudoJet &jet) const{
  // check the jet is appropriate for computation
  if (!jet.has_constituents())
    throw Error("ScalarPt can only be applied on jets for which the constituents are known.");
  
  vector<PseudoJet> constits = jet.constituents();
  double ptsum = 0.0;
  for (unsigned int i=0; i<constits.size(); i++) ptsum += constits[i].pt();
  return ptsum;
}


//------------------------------------------------------------------------
PseudoJet KtDij::partition(const PseudoJet &jet) const{
  // check the jet is appropriate for computation
  if (!jet.has_constituents()){
    throw Error("KtDij can only be applied on jets for which the constituents are known.");
  }
  
  vector<PseudoJet> constits = jet.constituents();
  if (constits.size()<2){
    return join(vector<PseudoJet>(2,PtYPhiM(1,0,0,0)));
  }  

  // compute the observable
  // recluster the constituents using kt
  vector<PseudoJet> particles, ghosts;
  SelectorIsPureGhost().sift(constits, ghosts, particles);
  double ghost_area = ghosts.size() ? ghosts[0].area() : 0.01;
  ClusterSequenceActiveAreaExplicitGhosts *cs
    = new ClusterSequenceActiveAreaExplicitGhosts(particles,
						  JetDefinition(kt_algorithm, 999.9),
						  ghosts, ghost_area);
  PseudoJet ktjet = SelectorNHardest(1)(cs->inclusive_jets())[0];

  // undo the last clustering and use that as a partition  
  PseudoJet parent1, parent2;
  ktjet.has_parents(parent1, parent2);
  cs->delete_self_when_unused();

  return join(parent1, parent2);
}

double KtDij::result_from_partition(const PseudoJet &partit) const{
  // ensure that we have a composite jet with 2 pieces
  if (!partit.has_pieces())
    throw Error("KtDij::result_from_partition can only be computed for composite jets");
 
  vector<PseudoJet> pieces = partit.pieces();
  if (pieces.size() != 2)
    throw Error("KtDij::result_from_partition can only be computed for composite jets made of 2 pieces");

  return pieces[0].kt_distance(pieces[1]);
}

//------------------------------------------------------------------------
string Angularity::description() const{
  ostringstream oss;
  oss << "Angularity with alpha=" << _alpha;
  return oss.str();
}

double Angularity::result(const PseudoJet &jet) const{
  // check the jet is appropriate for computation
  if (!jet.has_constituents())
    throw Error("Angularities can only be applied on jets for which the constituents are known.");
  
  vector<PseudoJet> constits = jet.constituents();
  double num=0.0, den=0.0;
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++){
    double pt = ci->pt();
    num += pt * pow(ci->squared_distance(jet), 1-_alpha/2);
    den += pt;
  }
  return num/den;
}


string AngularityNumerator::description() const{
  ostringstream oss;
  oss << "Angularity numerator with alpha=" << _alpha;
  return oss.str();
}

double AngularityNumerator::result(const PseudoJet &jet) const{
  // check the jet is appropriate for computation
  if (!jet.has_constituents())
    throw Error("Angularities can only be applied on jets for which the constituents are known.");

  vector<PseudoJet> constits = jet.constituents();
  double num=0.0;
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++)
    num += ci->pt() * pow(ci->squared_distance(jet), 1-_alpha/2);
  return num;
}


//------------------------------------------------------------------------
string TauEEC::description() const{
  ostringstream oss;
  oss << "Energy-energy correlator with beta=" << _beta;
  return oss.str();
}

double TauEEC::result(const PseudoJet &jet) const{
  vector<PseudoJet> constits = jet.constituents();
  double num=0.0, den=0.0;
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++){
    vector<PseudoJet>::iterator cj = constits.begin();
    double pti = ci->pt();
    while (cj!=ci){
      num += pti * cj->pt() * pow(ci->squared_distance(*cj), 0.5*_beta);
      cj++;
    }
    den += pti;
  }
  return num/(den*den);
}

string TauEECNumerator::description() const{
  ostringstream oss;
  oss << "Numerator of Energy-energy correlator with beta=" << _beta;
  return oss.str();
}
  
double TauEECNumerator::result(const PseudoJet &jet) const{
  vector<PseudoJet> constits = jet.constituents();
  double t=0.0;
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++){
    vector<PseudoJet>::iterator cj = constits.begin();
    while (cj!=ci){
      t += sqrt(ci->perp2() * cj->perp2()) * pow(ci->squared_distance(*cj), 0.5*_beta);
      cj++;
    }
  }
  return t;
}


//------------------------------------------------------------------------
PseudoJet NSubjettinessNumerator::partition(const PseudoJet &jet) const{
  if (! jet.has_constituents())
    throw Error("N-subjettiness can only be computed for jets with available constituents");

  // get the constituents  
  vector<PseudoJet> constituents = jet.constituents();
  
  // recluster them with a kt algorithm
  JetDefinition subdef(kt_algorithm, 999.0);
  vector<PseudoJet> particles, ghosts;
  SelectorIsPureGhost().sift(constituents, ghosts, particles);
  double ghost_area = ghosts.size() ? ghosts[0].area() : 0.01;
  ClusterSequenceActiveAreaExplicitGhosts *cs
    = new ClusterSequenceActiveAreaExplicitGhosts(particles, subdef,
						  ghosts, ghost_area);

  vector<PseudoJet> axes = cs->exclusive_jets_up_to(_N);
  cs->delete_self_when_unused();

  // take the axes as the N exclusive subjets
  return join(axes);
}

double NSubjettinessNumerator::result_from_partition(const PseudoJet &partit) const{
  // ensure that we have a composite jet with pieces
  if (!partit.has_pieces())
    throw Error("NSubjettinessNumerator::result_from_partition can only be computed for composite jets");
 
  vector<PseudoJet> axes = partit.pieces();
  if (axes.size() < _N) return 0.0;
  if (axes.size() > _N)
    throw Error("NSubjettinessNumerator::result_from_partition can only be computed for composite jets made of N pieces");
  
  // get the constituents  
  vector<PseudoJet> constituents = partit.constituents();
  
  // compute N-subjettiness from these axes and constituents
  double sum=0.0;
  for (unsigned int ic=0; ic<constituents.size(); ic++){
    const PseudoJet &c = constituents[ic];
    double mind2=std::numeric_limits<double>::max();
    for (unsigned int ia=0; ia<axes.size(); ia++)
      mind2 = std::min(mind2, c.squared_distance(axes[ia]));
    sum += sqrt(c.pt2() * mind2);
  }
  
  return sum;
}


} // namespace contrib

FASTJET_END_NAMESPACE
