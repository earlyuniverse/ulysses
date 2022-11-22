//----------------------------------------------------------------------
/// \file
/// \page Example09 09 - adding user information to a fastjet::PseudoJet
///
/// This example illustrates how it is possible to associate
/// user-defined information to a fastjet::PseudoJet (beyond the
/// simple user index), using a class derived from
/// fastjet::UserInfoBase.
///
/// Note that in this example we have chosen to use this
/// user-defined information to obtain properties of the constituents
/// of the reconstructed jet (e.g. if the event is made of a hard
/// interaction and pileup, what part of the reconstructed jet comes
/// from the hard interaction). To do that, we also show how to
/// introduce a user-defined fastjet::Selector. For some applications,
/// it might also be useful to define new recombination schemes using
/// the extra information.
///
/// run it with    : ./09-user_info < data/Pythia-dijet-ptmin100-lhc-pileup-1ev.dat
///
/// (Note that this event consists of many sub-events, the first one
/// being the "hard" interaction and the following being minbias
/// events composing the pileup. It has the specificity that it also
/// contains the PDG id of the particles)
///
/// Source code: 09-user_info.cc
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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include <iostream> // needed for io
#include <sstream>  // needed for io
#include <cstdio>   // needed for io

using namespace std;
using namespace fastjet;

//------------------------------------------------------------------------
// the user information
// 
// To associate extra information to a PseudoJet, one first has to
// create a class, derived from UserInfoBase, that contains
// that information.
//
// In our simple example, we shall use 2 informations
//  - the PDG id associated with the particle
//  - the "vertex number" associated with the particle
class MyUserInfo : public PseudoJet::UserInfoBase{
public:
  // default ctor
  //  - pdg_id        the PDG id of the particle
  //  - vertex_number theid of the vertex it originates from
  MyUserInfo(const int & pdg_id_in, const int & vertex_number_in) :
    _pdg_id(pdg_id_in), _vertex_number(vertex_number_in){}

  /// access to the PDG id
  int pdg_id() const { return _pdg_id;}
  
  /// access to the vertex number
  int vertex_number() const { return _vertex_number;}
  
protected:
  int _pdg_id;         // the associated pdg id
  int _vertex_number;  // the associated vertex number
};


//------------------------------------------------------------------------
// Select pi0 and photons
//
// This shows how we can build a Selector that uses the user-defined
// information to select particles that are either pi0's or photons
// (we choose this purely for simplicity).
// 
// To create a user-defined Selector, the first step is to
// create its associated "worker" class, i.e. to derive a class from
// SelectorWorker. Then (see below), we just write a function
// (SelectorIsPi0Gamma()) that creates a Selector with the
// appropriate worker class.
class SW_IsPi0Gamma : public SelectorWorker{
public:
  // default ctor
  SW_IsPi0Gamma(){}

  // the selector's description
  string description() const{
    return "neutral pions or photons"; 
  }

  // keeps the ones that have the pdg id of the pi0
  bool pass(const PseudoJet &p) const{
    // This is how we access the extra information associated with the
    // particles we test.
    //   p.user_info<MyUserInfo>()
    // returns a reference to the user-defined information (of type
    // MyUserInfo, as mentioned explicitly). It ensures automatically
    // that there is an associated user info compatible with the
    // requested type (and throws an error if it is not the case)
    //
    // We can then access the "pdg_id" member of MyUserInfo to
    // extract the targeted information.
    const int & pdgid = p.user_info<MyUserInfo>().pdg_id();
    return (pdgid == 111) || (pdgid == 22);
  }
};

// the function that allows to write simply
//    Selector sel = SelectorIsPi0Gamma();
Selector SelectorIsPi0Gamma(){
  return Selector(new SW_IsPi0Gamma());
}

//------------------------------------------------------------------------
// Select particles from a given vertex number
//
// This is the same kind of construct as just above except that we
// select on particles that are originated from a given vertex. The
// test event has been structured as a superposition of sub-events
// (the 0th being the hard interaction) and each particle will be
// associated a vertex number. This Selector allows to select
// particles corresponding to a given vertex number.
//
// As in the previous case, we start with the worker class and then
// write a function for the Selector itself.
class SW_VertexNumber : public SelectorWorker{
public:
  // ctor from the vertex we want to keep
  SW_VertexNumber(const int & vertex_number) : _vertex_number(vertex_number){}

  // the selector's description
  string description() const{
    ostringstream oss;
    oss << "vertex number " << _vertex_number;
    return oss.str();
  }

  // keeps the ones that have the correct vertex number
  bool pass(const PseudoJet &p) const{
    // This is how we access the extra information associated with the
    // particles we test.
    //   p.user_info<MyUserInfo>()
    // returns a reference to the user-defined information (of type
    // MyUserInfo, as mentioned explicitly). It ensures automatically
    // that there is an associated user info compatible with the
    // requested type (and throws an error if it is not the case)
    //
    // We can then access the "vertex_number" member of MyUserInfo to
    // extract the targeted information.
    return p.user_info<MyUserInfo>().vertex_number()==_vertex_number;
  }

private:
  int _vertex_number;  // the vertex number we're selecting
};

// The function that allows to write e.g.
//  Selector sel = !SelectorVertexNumber(0);
// to select particles from all vertices except the 0th.
Selector SelectorVertexNumber(const int & vertex_number){
  return Selector(new SW_VertexNumber(vertex_number));
}


//------------------------------------------------------------------------
// The example code associating user-info to the particles in the event
int main(){
  // read in input particles
  //----------------------------------------------------------
  vector<PseudoJet> input_particles;
  
  double px, py , pz, E;
  string str;
  int vertex_number=-1;
  int pdg_id = 21;
  while (getline(cin, str)){
    // if the line does not start with #, it's a particle
    // read its momentum and pdg id
    if (str[0] != '#'){
      istringstream iss(str);
      if (! (iss >> px >> py >> pz >> E >> pdg_id)){
	cerr << "Wrong file format: particles must be specified with" << endl;
	cerr << "  px py pz E pdg_id" << endl;
      }

      // first create a PseudoJet with the correct momentum
      PseudoJet p(px,py,pz,E);

      // associate to that our user-defined extra information
      // which is done using 
      //   PseudoJet::set_user_info()
      //
      // IMPORTANT NOTE: set_user_info(...) takes a pointer as an
      // argument. It will "own" that pointer i.e. will delete it when
      // all the PseudoJet's using it will be deleted.
      //
      // NB: once you've done p.set_user_info(my_user_info_ptr), you must
      // not call p2.set_user_info(my_user_info_ptr) with the same pointer
      // because p and p2 will both attempt to delete it when they go out
      // of scope causing a double-free corruption error. Instead do
      // p2.user_info_shared_ptr() = p.user_info_shared_ptr();
      p.set_user_info(new MyUserInfo(pdg_id, vertex_number));
      PseudoJet p2; // defined only to make the above documentation consistent!

      input_particles.push_back(p);
      continue;
    }

    // check if we start a new sub-event
    if ((str.length()>=9) && (str.compare(0,9,"#SUBSTART")==0)){
      vertex_number++;
    }
  }
  

  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  //----------------------------------------------------------
  double R = 0.6;
  JetDefinition jet_def(antikt_algorithm, R);


  // run the jet clustering with the above jet definition
  //----------------------------------------------------------
  ClusterSequence clust_seq(input_particles, jet_def);


  // get the resulting jets ordered in pt
  //----------------------------------------------------------
  double ptmin = 25.0;
  vector<PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));


  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt}
  //----------------------------------------------------------
  cout << "Ran " << jet_def.description() << endl;

  // label the columns
  printf("%5s %15s %15s %15s %15s %15s\n","jet #",
	 "rapidity", "phi", "pt",
	 "pt_hard", "pt_pi0+gamma");

  // a selection on the 1st vertex
  Selector sel_vtx0 = SelectorVertexNumber(0);

  // a selection on the pi0
  Selector sel_pi0gamma = SelectorIsPi0Gamma();

  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
    const PseudoJet & full = inclusive_jets[i];
    const vector<PseudoJet> constituents = full.constituents();

    // get the contribution from the 1st vertex
    PseudoJet hard = join(sel_vtx0(constituents));
    
    // get the contribution from the pi0's
    PseudoJet pi0gamma = join(sel_pi0gamma(constituents));
    
    // print the result
    printf("%5u %15.8f %15.8f %15.8f %15.8f %15.8f\n", i,
	   full.rap(), full.phi(), full.perp(),
	   hard.perp(), pi0gamma.perp());
  }

  return 0;
}
