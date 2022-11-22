//----------------------------------------------------------------------
/// \file
/// \page Example14 14 - unified use of transformers
///
/// fastjet example program to illustrate the use of the
/// fastjet::Filter and fastjet::Pruner classes in a unified way
/// through their derivation from fastjet::Transformer
///
/// The two hardest jets in a boosted top event, clustered with an
/// (abnormally) large R, are then groomed using different tools. One
/// notes the reduction in the mass of the jets after grooming.
///
/// run it with    : ./14-groomers < data/boosted_top_event.dat
///
/// Source code: 14-groomers.cc
//----------------------------------------------------------------------

//STARTHEADER
// $Id$
//
// Copyright (c) 2005-2018, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include <iostream>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <cstdio>   // needed for io

using namespace fastjet;
using namespace std;

/// an example program showing how to use Filter and Pruner in FastJet
int main(){
  // read in input particles
  //----------------------------------------------------------
  vector<PseudoJet> input_particles;
  
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a PseudoJet with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(PseudoJet(px,py,pz,E)); 
  }
 
  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  JetDefinition jet_def(cambridge_algorithm, 1.5);
  ClusterSequence clust_seq(input_particles, jet_def);
  vector<PseudoJet> inclusive_jets = 
                             sorted_by_pt(clust_seq.inclusive_jets(5.0));

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt", "mass");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp(),inclusive_jets[i].m());
  }

  // simple test to avoid that the example below crashes:
  // make sure there is at least 3 jets above our 5 GeV
  if (inclusive_jets.size()<2){
    cout << "Please provide an event with at least 2 jets above 5 GeV" << endl;
    return 1;
  }

  // We will groom the two hardest jets of the event
  //----------------------------------------------------------
  vector<PseudoJet> candidates = SelectorNHardest(2)(inclusive_jets);

  // create 3 groomers
  //----------------------------------------------------------
  vector<Transformer *> groomers;
  
  // 1.
  // the Cambridge/Aachen filter with Rfilt=0.3 
  // (simplified version of arXiv:0802.2470)
  double Rfilt = 0.3;
  unsigned int nfilt = 3;
  groomers.push_back(new Filter(JetDefinition(cambridge_algorithm, Rfilt), 
                                SelectorNHardest(nfilt) ) );

  // 2.
  // Filtering with a pt cut as for trimming (arXiv:0912.1342)
  double Rtrim = 0.2;
  double ptfrac = 0.03;
  groomers.push_back(new Filter(JetDefinition(kt_algorithm, Rtrim), 
                                SelectorPtFractionMin(ptfrac) ) );

  // 3.
  // Pruning (arXiv:0903.5081)
  double zcut = 0.1;
  double rcut_factor = 0.5;
  groomers.push_back(new Pruner(cambridge_algorithm, zcut, rcut_factor));

  // apply the various groomers to the test PseudoJet's
  // and show the result
  //----------------------------------------------------------

  // print out original jet candidates
  cout << "\nOriginal jets that will be grooomed: " << endl;
  for (vector<PseudoJet>::iterator jit=candidates.begin(); jit!=candidates.end(); jit++){
    const PseudoJet & c = *jit;
    cout << "  rap = " << c.rap() << ", phi = " << c.phi() << ", pt = " << c.perp() 
         << ", mass = " << c.m() 
         << "  [" << c.description() << "]" <<  endl;
  }

  // loop on groomers
  for (unsigned int i=0; i < groomers.size(); i++){
    const Transformer & f = *groomers[i];
    cout << "\nUsing groomer: " << f.description() << endl;
    
    // loop on jet candidates
    for (vector<PseudoJet>::iterator jit=candidates.begin(); jit!=candidates.end(); jit++){
      const PseudoJet & c = *jit;
      
      // apply groomer f to jet c      
      PseudoJet j = f(c);
      
      // access properties specific to the given transformer
      //
      // We first make sure that the jet indeed has a structure
      // compatible with the result of a Filter or Pruner (using
      // has_structure_of()), and then retrieve the pieces rejected by the
      // groomer (using structure_of())
      int n_rejected;         
      if (j.has_structure_of<Filter>()) {
         const Filter::StructureType & fj_struct = j.structure_of<Filter>(); 
	 n_rejected = fj_struct.rejected().size();
      }else { 
         assert(j.has_structure_of<Pruner>());  // make sure
         const Pruner::StructureType & fj_struct = j.structure_of<Pruner>();
	 n_rejected = fj_struct.rejected().size();
      }
      
      // write out result
      cout << "  rap = " << j.rap() << ", phi = " << j.phi() << ", pt = " << j.perp()
           << " mass = " << j.m() 
           << "  [kept: " << j.pieces().size() 
           << ", rejected: " << n_rejected;
      if (j.has_structure_of<Pruner>()) {
        cout << ", Rcut: " << j.structure_of<Pruner>().Rcut();
      }
      cout << "]" << endl;
    }
  }

  // a bit of memory cleaning
  for (unsigned int i=0; i < groomers.size(); i++) delete groomers[i];

  return 0;
}
