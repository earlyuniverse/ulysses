//----------------------------------------------------------------------
/// \file
/// \page Example12 12 - boosted Higgs tagging
///
/// fastjet example program, illustration of carrying out boosted
/// Higgs subjet ID analysis
///
/// It illustrates two kinds of functionality: 
///
///  - using a boosted higgs tagger
///  - following information on a b-tag through the jet
///
/// This kind of functionality was used in arXiv:0802.2470
/// (Butterworth, Davison, Rubin & Salam) for boosted Higgs searches,
/// and related functionality was used in arXiv:0806.0848 (Kaplan,
/// Rehermann, Schwartz & Tweedie) in searching for boosted tops
/// (without b-tag assumptions).
///
/// run it with    : ./12-boosted_higgs < data/HZ-event-Hmass115.dat
///
/// Source code: 12-boosted_higgs.cc
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

#include <iostream> // needed for io
#include <sstream>  // needed for internal io
#include <iomanip>  
#include <cmath>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/tools/MassDropTagger.hh>
#include <fastjet/tools/Filter.hh>

using namespace std;
using namespace fastjet;


//----------------------------------------------------------------------
// set up a class to give standard (by default E-scheme)
// recombination, with additional tracking of flavour information in
// the user_index. 
//
// b-tagged particles are assumed to have their user_index set to 1,
// and other particles should have user_index to 0.
//
// Watch out however that, by default, the user_index of a particle is
// set to -1 and you may not have control over that (e.g. if you
// compute the jet area using explicit ghosts, the ghosts will have a
// default user_index of -1). For that reason, if one of the particle
// being combined has a user index of -1, we assume it is not b-tagged
// (i.e. we count it as 0 in the recombination)
//
// This will work for native algorithms, but not for all plugins
//----------------------------------------------------------------------
typedef JetDefinition::DefaultRecombiner DefRecomb;

class FlavourRecombiner : public  DefRecomb {
public:
  FlavourRecombiner(RecombinationScheme recomb_scheme = E_scheme) : 
    DefRecomb(recomb_scheme) {};

  virtual std::string description() const {
    return DefRecomb::description()+" (with user index addition)";}

  /// recombine pa and pb and put result into pab
  virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                         PseudoJet & pab) const {
    DefRecomb::recombine(pa,pb,pab);
    // Note: see the above discussion for the fact that we consider
    // negative user indices as "0"
    pab.set_user_index(max(pa.user_index(),0) + max(pb.user_index(),0));
  }
};


//----------------------------------------------------------------------
// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);


//----------------------------------------------------------------------
// core of the program
//----------------------------------------------------------------------
int main(){

  vector<PseudoJet> particles;

  // read in data in format px py pz E b-tag [last of these is optional]
  // lines starting with "#" are considered as comments and discarded
  //----------------------------------------------------------

  string line;
  while (getline(cin,line)) {
    if (line.substr(0,1) == "#") {continue;}
    istringstream linestream(line);
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;

    // optionally read in btag information
    int    btag;
    if (! (linestream >> btag)) btag = 0;

    // construct the particle
    PseudoJet particle(px,py,pz,E);
    particle.set_user_index(btag); // btag info goes in user index, for flavour tracking
    particles.push_back(particle);
  }


  // set up the jet finding
  //
  // This also shows how to use the "FlavourRecombiner" user-defined
  // recombiner 
  // ----------------------------------------------------------
  double R = 1.2;
  FlavourRecombiner flav_recombiner; // for tracking flavour
  JetDefinition jet_def(cambridge_algorithm, R, &flav_recombiner);

  
  // run the jet finding; find the hardest jet
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  cout << "Ran: " << jet_def.description() << endl << endl;
  cout << "Hardest jet: " << jets[0] << endl << endl;

  // now do jet tagging using a mass drop tagger
  //
  // Note: if you prefer, you may as well use a CASubJetTagger
  //    CASubJetTagger ca_tagger;
  //    PseudoJet tagged = ca_tagger(jets[0]);
  // This requires including fastjet/tools/CASubJetTagger.hh 
  // You also need to adapt the 2 lines below accessing
  // the extra structural information provided by the tagger
  //----------------------------------------------------------
  MassDropTagger md_tagger(0.667, 0.09);
  PseudoJet tagged = md_tagger(jets[0]);

  if (tagged == 0){
    cout << "No substructure found" << endl;
    return 0;
  }

  PseudoJet parent1 = tagged.pieces()[0];
  PseudoJet parent2 = tagged.pieces()[1];
  cout << "Found suitable pair of subjets: " << endl;
  cout << " " << parent1 << endl;
  cout << " " << parent2 << endl;
  cout << "Total = " << endl;
  cout << " " << tagged << endl;
  cout << "(mass drop = " << tagged.structure_of<MassDropTagger>().mu()
       << ", y = " << tagged.structure_of<MassDropTagger>().y() << ")" 
       << endl << endl;

  // next we "filter" it, to remove UE & pileup contamination
  //----------------------------------------------------------
  //
  // [there are two ways of doing this; here we directly use the
  // existing cluster sequence and find the exclusive subjets of
  // this_jet (i.e. work backwards within the cs starting from
  // this_jet); alternatively one can recluster just the
  // constituents of the jet]
  //
  // first get separation between the subjets (called Rbb -- assuming
  // it's a Higgs!)
  //
  // See example 11-filter.cc for another way of implementing the dynamic
  // Rfilt used below
  double   Rbb = parent1.delta_R(parent2);
  double   Rfilt = min(Rbb/2, 0.3); // somewhat arbitrary choice
  unsigned nfilt = 3;               // number of pieces we'll take
  cout << "Subjet separation (Rbb) = " << Rbb << ", Rfilt = " << Rfilt << endl;

  Filter filter(JetDefinition(cambridge_algorithm, Rfilt, &flav_recombiner),
		SelectorNHardest(nfilt));
  PseudoJet filtered = filter(tagged);

  // now print out the filtered jets and reconstruct total 
  // at the same time
  const vector<PseudoJet> & filtered_pieces = filtered.pieces();
  cout << "Filtered pieces are " << endl;
  for (unsigned i = 0; i < nfilt && i < filtered_pieces.size(); i++) {
    cout << " " << filtered_pieces[i] << endl;
  }
  cout << "Filtered total is " << endl;
  cout << " " << filtered << endl;

}


//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "pt, y, phi =" 
       << " " << setw(10) << jet.perp() 
       << " " << setw(6) <<  jet.rap()  
       << " " << setw(6) <<  jet.phi()  
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
}
