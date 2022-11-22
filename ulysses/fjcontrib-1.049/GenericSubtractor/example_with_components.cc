//----------------------------------------------------------------------
// Example on how to use this contribution
//
// run it with
//  ./example_with_components < ../data/Pythia-Zp2jets-lhc-pileup-1ev.dat
//----------------------------------------------------------------------

// $Id: example_with_components.cc 859 2015-09-21 10:11:32Z gsalam $
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

#include <iostream>
#include <iomanip>
#include <limits>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"

#include "ShapeWithComponents.hh"
#include "ExampleShapes.hh"
#include "GenericSubtractor.hh"

using namespace std;
using namespace fastjet;

// an example of a ShapeWithComponents
//
// We define tau_N/tau_{N-1} as the explicit ratio of (the numerator
// of) tau_N and (the numerator of) tau_{N-1}. When this shape will be
// subtracted, tau_N and tau_{N-1} will be subtracted independently
class NSubjettinessRatio : public contrib::ShapeWithComponents{
public:
  NSubjettinessRatio(int N) : _N(N){
    assert(_N>1);
  }

  // a (rather loosy) description
  virtual std::string description() const{ return "N-subjettiness ratio from components";}


  /// returns the number of components 
  virtual unsigned int n_components() const { return 2;}

  /// computes individually tau_N and tau_{N-1}
  virtual std::vector<double> components(const PseudoJet &jet) const{
    vector<double> comp(n_components());
    comp[0] = contrib::NSubjettinessNumerator(_N)(jet);
    comp[1] = contrib::NSubjettinessNumerator(_N-1)(jet);
    return comp;
  }

  /// given the components, determine the result of the event shape
  virtual double result_from_components(const std::vector <double> &components) const{
    return components[0]/components[1];
  }

  /// since the components are shapes with artition, take the
  /// advantage of it.
  ///
  /// This is done as follows (and is not necessary if no partition is
  /// needed) : the overloaded method below return a shape that
  /// computes the ith component.
  virtual FunctionOfPseudoJet<double> * component_shape(unsigned int index) const{
    return new contrib::NSubjettinessNumerator(_N-index);
  }

protected:
  const int _N;
};

// fwd declaration
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event);

//----------------------------------------------------------------------
int main(){
  // read in input particles
  //
  // since we use here simulated data we can split the hard event
  // from the full (i.e. with pileup added) one
  //
  // (see also example 07 in FastJet)
  //----------------------------------------------------------
  vector<PseudoJet> hard_event, full_event;
  read_event(hard_event, full_event);

  // keep the particles up to 4 units in rapidity
  hard_event = SelectorAbsRapMax(4.0)(hard_event);
  full_event = SelectorAbsRapMax(4.0)(full_event);
  
  // create what we need for the clustering
  //----------------------------------------------------------
  JetDefinition jet_def(antikt_algorithm, 0.7);
  JetDefinition jet_def_for_rho(kt_algorithm, 0.4);
  AreaDefinition area_def(active_area_explicit_ghosts,
                          GhostedAreaSpec(SelectorAbsRapMax(4.0)));
  Selector rho_range =  SelectorAbsRapMax(3.0);
  JetMedianBackgroundEstimator bge(rho_range, jet_def_for_rho, area_def);
  Subtractor subtractor(&bge);

  Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);

  // do the clustering
  //----------------------------------------------------------
  ClusterSequenceArea clust_seq_hard(hard_event, jet_def, area_def);
  ClusterSequenceArea clust_seq_full(full_event, jet_def, area_def);

  vector<PseudoJet> hard_jets = sel_jets(clust_seq_hard.inclusive_jets());
  vector<PseudoJet> full_jets = sel_jets(clust_seq_full.inclusive_jets());

  // the shape part
  //----------------------------------------------------------
  NSubjettinessRatio tau21(2);
  contrib::GenericSubtractor gen_sub(&bge);
  bge.set_particles(full_event);

  cout << "# original hard jets" << endl;
  for (unsigned int i=0; i<hard_jets.size(); i++){
    const PseudoJet &jet = hard_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", tau21 = " << tau21(jet) << endl;
  }
  cout << endl;

  cout << "# unsubtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];
    cout << "pt = " << jet.pt()
	 << ", rap = " << jet.rap()
	 << ", tau21 = " << tau21(jet) << endl;
  }
  cout << endl;

  cout << "# subtracted full jets" << endl;
  for (unsigned int i=0; i<full_jets.size(); i++){
    const PseudoJet &jet = full_jets[i];

    // subtract the jet and the shape value
    PseudoJet subtracted_jet = subtractor(jet);
    double subtracted_shape = gen_sub(tau21, jet);

    cout << "pt = " << subtracted_jet.pt()
	 << ", rap = " << subtracted_jet.rap()
	 << ", tau21 = " << subtracted_shape << endl;
  }
  cout << endl;

  return 0;
}

//------------------------------------------------------------------------
// read the event with and without pileup
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event){
  string line;
  int  nsub  = 0; // counter to keep track of which sub-event we're reading
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {break;}
    if (line.substr(0,9) == "#SUBSTART") {
      // if more sub events follow, make copy of first one (the hard one) here
      if (nsub == 1) hard_event = full_event;
      nsub += 1;
    }
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    full_event.push_back(particle);
  }

  // if we have read in only one event, copy it across here...
  if (nsub == 1) hard_event = full_event;

  // if there was nothing in the event 
  if (nsub == 0) {
    cerr << "Error: read empty event\n";
    exit(-1);
  }

  cout << "# " << nsub-1 << " pileup events on top of the hard event" << endl;
}
