//----------------------------------------------------------------------
/// \file
/// \page Example11 11 - use of filtering
///
/// fastjet example program to illustrate the use of the
/// fastjet::Filter class
///
/// We apply different filter examples to either the hardest jet of
/// the given event, or to the composition of the two hardest jets:
///
///   - two examples of a filter keeping a fixed number of subjets (as 
///     in arXiv:0802.2470)
///   - a "trimmer" i.e. a filter keeping subjets carrying at least a given 
///     fraction of the pt of the jet (arXiv:0912.1342).
///   - two examples of filter in combination with background subtraction
///
/// run it with    : ./11-filter < data/single-event.dat
///
/// Source code: 11-filter.cc
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
// the following includes are only needed when combining filtering with subtraction
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"

#include <cstdio>   // needed for io

using namespace fastjet;
using namespace std;

// a function returning
//   min(Rmax, deltaR_factor * deltaR(j1,j2))
// where j1 and j2 are the 2 subjets of j
// if the jet does not have 2 exactly pieces, Rmax is used.
class DynamicRfilt : public FunctionOfPseudoJet<double>{
public:
  // default ctor 
  DynamicRfilt(double Rmax, double deltaR_factor) : _Rmax(Rmax), _deltaR_factor(deltaR_factor){}

  // action of the function
  double result(const PseudoJet &j) const{
    if (! j.has_pieces()) return _Rmax;

    vector<PseudoJet> pieces = j.pieces();
    if (pieces.size() != 2) return _Rmax;

    double deltaR = pieces[0].delta_R(pieces[1]);
    return min(_Rmax, _deltaR_factor * deltaR);
  }

private:
  double _Rmax, _deltaR_factor;
};

/// an example program showing how to use Filter in FastJet
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
  JetDefinition jet_def(cambridge_algorithm, 1.2);
  // the use of a ClusterSequenceArea (instead of a plain ClusterSequence)
  // is only needed because we will later combine filtering with area-based
  // subtraction
  ClusterSequenceArea clust_seq(input_particles, jet_def, 
                                AreaDefinition(active_area_explicit_ghosts));
  vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(5.0));

  // label the columns
  printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    printf("%5u %15.8f %15.8f %15.8f\n",
	   i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
	   inclusive_jets[i].perp());
  }

  // simple test to avoid that the example below crashes:
  // make sure there is at least 3 jets above our 5 GeV
  if (inclusive_jets.size()<3){
    cout << "Please provide an event with at least 3 jets above 5 GeV" << endl;
    return 1;
  }

  // the sample PseudoJet that we will filter
  //  - the hardest jet of the event
  //  - the composition of the second and third hardest jets 
  ///   (this shows that the Filter can also be applied to a composite jet)
  //----------------------------------------------------------
  vector<PseudoJet> candidates;
  candidates.push_back(inclusive_jets[0]);
  candidates.push_back(join(inclusive_jets[1],inclusive_jets[2]));


  // create 5 filters
  //----------------------------------------------------------
  vector<Filter> filters;
  
  // 1.
  // the Cambridge/Aachen filter with Rfilt=0.3 (simpliefied version of arXiv:0802.2470)
  filters.push_back(Filter(JetDefinition(cambridge_algorithm, 0.3), SelectorNHardest(3)));

  // 2.
  // the Cambridge/Aachen filter with Rfilt=min(0.3, 0.5*Rbb) as in arXiv:0802.2470
  SharedPtr<DynamicRfilt> dynamic_Rfilt(new DynamicRfilt(0.3, 0.5));
  filters.push_back(Filter(dynamic_Rfilt.get(), SelectorNHardest(3)));

  // 3.
  // Filtering with a pt cut as for trimming (arXiv:0912.1342)
  filters.push_back(Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03)));

  // 4.
  // First example of filtering with subtraction of the background: provide rho
  // First, estimate the background for the given event
  GridMedianBackgroundEstimator bkgd(4.5, 0.55); // uses particles up to |y|=4.5
  bkgd.set_particles(input_particles);
  double rho = bkgd.rho();
  // Then, define the filter
  filters.push_back(Filter(JetDefinition(cambridge_algorithm, 0.3), SelectorNHardest(3), rho));

  // 5.
  // Second example of filtering with subtraction of the background: set a subtractor
  // First, define a subtractor from a background estimator
  Subtractor subtractor(&bkgd);
  // Then, define the filter
  Filter filt(JetDefinition(cambridge_algorithm, 0.3), SelectorNHardest(3));
  // Finally, tell the filter about the subtractor
  filt.set_subtractor(&subtractor);
  filters.push_back(filt);


  // apply the various filters to the test PseudoJet
  // and show the result
  //----------------------------------------------------------

  // print out original jet candidates
  cout << "\nOriginal jets that will be filtered: " << endl;
  for (vector<PseudoJet>::iterator jit=candidates.begin(); jit!=candidates.end(); jit++){
    const PseudoJet & c = *jit;
    cout << "  rap = " << c.rap() << ", phi = " << c.phi() << ", pt = " << c.perp() 
         << "  [" << c.description() << "]" <<  endl;
  }

  // loop on filters
  for (vector<Filter>::iterator it=filters.begin(); it!=filters.end(); it++){
    const Filter & f = *it;
    cout << "\nUsing filter: " << f.description() << endl;
    
    // loop on jet candidates
    for (vector<PseudoJet>::iterator jit=candidates.begin(); jit!=candidates.end(); jit++){
      const PseudoJet & c = *jit;
      
      // apply filter f to jet c      
      PseudoJet j = f(c);
      
      // access properties specific to the Filter
      //
      // We first make sure that the jet indeed has a structure
      // compatible with the result of a Filter (using
      // has_structure_of()), and then retrieve the pieces rejected by the
      // filter (using structure_of())
      assert(j.has_structure_of<Filter>());
      const Filter::StructureType & fj_struct = j.structure_of<Filter>();
      
      // write out result
      cout << "  rap = " << j.rap() << ", phi = " << j.phi() << ", pt = " << j.perp() 
           << "  [kept: " << j.pieces().size() << ", rejected: "
	   << fj_struct.rejected().size() << "]" << endl;
    }
  }

  return 0;
}
