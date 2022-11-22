//----------------------------------------------------------------------
/// \file example_advanced_usage.cc
///
/// This example program is meant to illustrate some advanced features
/// of the fastjet::contrib::SoftDrop class is used.
///
/// Run this example with
///
/// \verbatim
///     ./example_advanced_usage < ../data/single-event.dat
/// \endverbatim
//----------------------------------------------------------------------

// $Id: example_advanced_usage.cc 1016 2017-04-20 16:51:52Z knzhou $
//
// Copyright (c) 2014, Gavin P. Salam
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
#include <sstream>

#include <sstream>
#include <iomanip>
#include <cmath>
#include "fastjet/ClusterSequence.hh"
#include "SoftDrop.hh" // In external code, this should be fastjet/contrib/SoftDrop.hh

#define MY_INF std::numeric_limits<double>::infinity()

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
ostream & operator<<(ostream &, const PseudoJet &);

// Simple class to store SoftDrop objects along with display information
class SoftDropStruct {

private:
  string _name;
  double _beta;
  double _z_cut;
  SoftDrop::SymmetryMeasure _symmetry_measure;
  double _R0;
  double _mu;  //mass drop
  SoftDrop::RecursionChoice  _recursion_choice;
  JetAlgorithm _recluster_algorithm;
  bool _tagging_mode;
  SoftDrop _soft_drop;
  Recluster _reclusterer;
  
public:
  SoftDropStruct(string name,
                 double beta,
                 double z_cut,
                 SoftDrop::SymmetryMeasure symmetry_measure,
                 double R0,
                 double mu,
                 SoftDrop::RecursionChoice recursion_choice,
                 JetAlgorithm recluster_algorithm,
                 bool tagging_mode = false)
  :  _name(name),
    _beta(beta),
    _z_cut(z_cut),
    _symmetry_measure(symmetry_measure),
    _R0(R0),
    _mu(mu),
    _recursion_choice(recursion_choice),
    _recluster_algorithm(recluster_algorithm),
    _tagging_mode(tagging_mode),
    _soft_drop(beta,z_cut,symmetry_measure,R0,mu,recursion_choice),
    _reclusterer(recluster_algorithm,JetDefinition::max_allowable_R)
  {
    // no need to recluser if already CA algorithm
    if (recluster_algorithm != cambridge_algorithm) {
      _soft_drop.set_reclustering(true,&_reclusterer);
    }
    
    // if beta is negative, typically want to use in tagging mode
    // MMDT behavior is also tagging mode
    // set this option here
    if (tagging_mode) {
      _soft_drop.set_tagging_mode();
    }
    
    //turn verbose structure on (off by default)
    _soft_drop.set_verbose_structure();
  }
  
  const SoftDrop& soft_drop() const { return _soft_drop;}
  string name() const { return _name;}
  double beta() const { return _beta;}
  double z_cut() const { return _z_cut;}
  string symmetry_measure_name() const {
    switch (_symmetry_measure) {
      case SoftDrop::scalar_z:
        return "scalar_z";
      case SoftDrop::vector_z:
        return "vector_z";
      case SoftDrop::y:
        return "y";
      default:
        return "unknown";
    }
  }
  double R0() const { return _R0;}
  double mu() const { return _mu;}

  string recursion_choice_name() const {
    switch (_recursion_choice) {
      case SoftDrop::larger_pt:
        return "larger_pt";
      case SoftDrop::larger_mt:
        return "larger_mt";
      case SoftDrop::larger_m:
        return "larger_m";
      default:
        return "unknown";
    }
  }

  string reclustering_name() const {
    switch (_recluster_algorithm) {
      case kt_algorithm:
        return "KT";
      case cambridge_algorithm:
        return "CA";
      default:
        return "unknown";
    }
  }
  
  string tag_or_groom() const {
    return (_tagging_mode ? "tag" : "groom");
  }
  
};


//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;
  
  // first get some anti-kt jets
  double R = 1.0, ptmin = 20.0;
  JetDefinition jet_def(antikt_algorithm, R);
  ClusterSequence cs(event, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptmin));

  
  //----------------------------------------------------------
  // Make vector of structs to store a bunch of different SoftDrop options
  // This is a vector of pointers because SoftDropStruct doesn't
  // have a valid copy constructor
  vector<SoftDropStruct *> sd_vec;
  
  // make some standard SoftDrop
  sd_vec.push_back(new SoftDropStruct("beta=2.0 zcut=.1",  2.0,0.10,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("beta=1.0 zcut=.1",  1.0,0.10,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("beta=0.5 zcut=.1",  0.5,0.10,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("beta=2.0 zcut=.2",  2.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("beta=1.0 zcut=.2",  1.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("beta=0.5 zcut=.2",  0.5,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  
  // make a mMDT-like tagger using SoftDrop
  sd_vec.push_back(new SoftDropStruct("MMDT-like zcut=.1", 0.0,0.10,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  sd_vec.push_back(new SoftDropStruct("MMDT-like zcut=.2", 0.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  sd_vec.push_back(new SoftDropStruct("MMDT-like zcut=.3", 0.0,0.30,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  sd_vec.push_back(new SoftDropStruct("MMDT-like zcut=.4", 0.0,0.40,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  
  // make some tagging SoftDrop (negative beta)
  sd_vec.push_back(new SoftDropStruct("beta=-2.0 zcut=.05", -2.0,0.05,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  sd_vec.push_back(new SoftDropStruct("beta=-1.0 zcut=.05", -1.0,0.05,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  sd_vec.push_back(new SoftDropStruct("beta=-0.5 zcut=.05", -0.5,0.05,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm,true));
  
  // make a SoftDrop with R0 parameter
  sd_vec.push_back(new SoftDropStruct("b=.5 z=.3 R0=1.0",     0.5,0.30,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=.5 z=.3 R0=0.5",     0.5,0.30,SoftDrop::scalar_z,0.5,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=.5 z=.3 R0=0.2",     0.5,0.30,SoftDrop::scalar_z,0.2,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));

  // make a SoftDrop with different symmetry measure
  sd_vec.push_back(new SoftDropStruct("b=2 z=.4 scalar_z",   2.0,0.4,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=2 z=.4 vector_z",   2.0,0.4,SoftDrop::vector_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=2 z=.4 y",          2.0,0.4,SoftDrop::y,       1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));

  // make a SoftDrop with different recursion choice
  sd_vec.push_back(new SoftDropStruct("b=3 z=.2 larger_pt",   3.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=3 z=.2 larger_mt",   3.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_mt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=3 z=.2 larger_m",    3.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_m,cambridge_algorithm));

  // make a SoftDrop with mass drop
  sd_vec.push_back(new SoftDropStruct("b=2 z=.1 mu=1.0",    2.0,0.10,SoftDrop::scalar_z,1.0,1.0,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=2 z=.1 mu=0.8",    2.0,0.10,SoftDrop::scalar_z,1.0,0.8,SoftDrop::larger_pt,cambridge_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=2 z=.1 mu=0.5",    2.0,0.10,SoftDrop::scalar_z,1.0,0.5,SoftDrop::larger_pt,cambridge_algorithm));
  
  // make a SoftDrop with a different clustering scheme (kT instead of default CA)
  sd_vec.push_back(new SoftDropStruct("b=2.0 z=.2 kT",        2.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=1.0 z=.2 kT",        1.0,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=0.5 z=.2 kT",        0.5,0.20,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=2.0 z=.4 kT",        2.0,0.40,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=1.0 z=.4 kT",        1.0,0.40,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  sd_vec.push_back(new SoftDropStruct("b=0.5 z=.4 kT",        0.5,0.40,SoftDrop::scalar_z,1.0,MY_INF,SoftDrop::larger_pt,kt_algorithm));
  
  //----------------------------------------------------------
  // Output information about the various Soft Drop algorithms

  
  // header lines
  cout << "---------------------------------------------------------------------------------------------" << endl;
  cout << "Soft Drops to be tested:" << endl;
  cout << "---------------------------------------------------------------------------------------------" << endl;
  cout << std::setw(18) << "name"
    << std::setw(8) << "beta"
    << std::setw(8) << "z_cut"
    << std::setw(9) << "sym"
    << std::setw(8) << "R0"
    << std::setw(8) << "mu"
    << std::setw(10)<< "recurse"
    << std::setw(8) << "reclust"
    << std::setw(8) << "mode"
    << endl;
  
  // set precision for display
  cout << setprecision(3) << fixed;
  
  // line for each SoftDrop
  for (unsigned i_sd = 0; i_sd < sd_vec.size(); i_sd++) {
    cout << std::setw(18) << sd_vec[i_sd]->name()
    << std::setw(8) << sd_vec[i_sd]->beta()
    << std::setw(8) << sd_vec[i_sd]->z_cut()
    << std::setw(9) << sd_vec[i_sd]->symmetry_measure_name()
    << std::setw(8) << sd_vec[i_sd]->R0()
    << std::setw(8) << sd_vec[i_sd]->mu()
    << std::setw(10)<< sd_vec[i_sd]->recursion_choice_name()
    << std::setw(8) << sd_vec[i_sd]->reclustering_name()
    << std::setw(8) << sd_vec[i_sd]->tag_or_groom()
    << endl;
  }
  cout << "---------------------------------------------------------------------------------------------" << endl;



  for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
  
      cout << "---------------------------------------------------------------------------------------------" << endl;
      cout << "Analyzing Jet " << ijet + 1 << ":" << endl;
      cout << "---------------------------------------------------------------------------------------------" << endl;
    
      cout << std::setw(18) << "name"
      << std::setw(10) << "pt"
      << std::setw(9) << "m"
      << std::setw(8) << "y"
      << std::setw(8) << "phi"
      << std::setw(8) << "constit"
      << std::setw(8) << "delta_R"
      << std::setw(8) << "sym"
      << std::setw(8) << "mu"
      << std::setw(8) << "mxdropz"  // max_dropped_z from verbose logging
    
    << endl;
    
    PseudoJet original_jet = jets[ijet];
    
    // set precision for display
    cout << setprecision(4) << fixed;
    
    cout << std::setw(18) << "Original Jet"
      << std::setw(10) << original_jet.pt()
      << std::setw(9) << original_jet.m()
      << std::setw(8) << original_jet.rap()
      << std::setw(8) << original_jet.phi()
      << std::setw(8) << original_jet.constituents().size()
    << endl;
    
    for (unsigned i_sd = 0; i_sd < sd_vec.size(); i_sd++) {
      // the current Soft Drop
      const SoftDrop & sd = (*sd_vec[i_sd]).soft_drop();
      PseudoJet sd_jet = sd(jets[ijet]);
      
      cout << std::setw(18) << sd_vec[i_sd]->name();
      if (sd_jet != 0.0) {
        cout << std::setw(10) << sd_jet.pt()
          << std::setw(9) << sd_jet.m()
          << std::setw(8) << sd_jet.rap()
          << std::setw(8) << sd_jet.phi()
          << std::setw(8) << sd_jet.constituents().size()
          << std::setw(8) << sd_jet.structure_of<contrib::SoftDrop>().delta_R()
          << std::setw(8) << sd_jet.structure_of<contrib::SoftDrop>().symmetry()
          << std::setw(8) << sd_jet.structure_of<contrib::SoftDrop>().mu()
          // the next line is part of the verbose information that is off by default
          << std::setw(8) << sd_jet.structure_of<contrib::SoftDrop>().max_dropped_symmetry();
      } else {
        cout << " ---- untagged jet ----";
      }
      cout << endl;

    }
  }

  // clean up
  for (unsigned i_sd = 0; i_sd < sd_vec.size(); i_sd++) {
    delete sd_vec[i_sd];
  }

  return 0;
}

//----------------------------------------------------------------------
/// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}

//----------------------------------------------------------------------
/// overloaded jet info output
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  if (jet == 0) {
    ostr << " 0 ";
  } else {
    ostr << " pt = " << jet.pt()
         << " m = " << jet.m()
         << " y = " << jet.rap()
         << " phi = " << jet.phi();
  }
  return ostr;
}
