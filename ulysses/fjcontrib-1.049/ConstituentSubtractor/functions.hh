// $Id: functions.hh 1240 2020-02-23 13:51:05Z peter.berta $
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
#include <vector>
#include <iomanip>

#include "fastjet/PseudoJet.hh"

using namespace fastjet;
using namespace std;

// \class JetWidth
// computes the jet width used as one shape in the example
class JetWidth : public fastjet::FunctionOfPseudoJet<double>{
public:

  // action of the function
  double result(const fastjet::PseudoJet &jet) const{
    if (!jet.has_constituents()){
      return -0.1;
    }    
    double width = 1e-6;
    double ptSum = 0;
    
    if (jet.constituents().size() < 2) return width;

    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    for (unsigned int iconst=0; iconst<constituents.size(); iconst++){
      fastjet::PseudoJet cons=constituents.at(iconst);
      double dR = sqrt(cons.squared_distance(jet));
      //cout << iconst << "  " << cons.squared_distance(jet) << "  " << cons.pt() << "  " << cons.m() << endl;
      double pt = cons.pt();
      width += dR * pt;
      ptSum += pt;
    }

    return width/ptSum;
  }

};



// read in input particles                                                                                                                                                              
void read_event(vector<PseudoJet> &hard_event, vector<PseudoJet> &full_event, vector<PseudoJet> *hard_event_charged=0, vector<PseudoJet> *pileup_event_charged=0){
  string line;
  int  nsub  = 0; // counter to keep track of which sub-event we're reading                                                                                                             
  bool format_ptRapPhi=false;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there are extra "pollution"                                                                                                               
    // characters (e.g. line-feed).                                                                                                                     

    if (line.substr(0,9) == "#PTRAPPHI") format_ptRapPhi=true;
                                
    if (line.substr(0,4) == "#END") {break;}
    if (line.substr(0,9) == "#SUBSTART") {
      // if more sub events follow, make copy of first one (the hard one) here                                                                                                          
      if (nsub == 1) hard_event = full_event;
      nsub += 1;
    }
    if (line.substr(0,1) == "#") {continue;}
    PseudoJet particle(0,0,1,1);
    double charge,pid;

    if (format_ptRapPhi){
      double pt,y,phi;
      linestream >> pt >> y >> phi >> pid >> charge;
      particle.reset_PtYPhiM(pt,y,phi);
    }
    else{
      double px,py,pz,E;
      linestream >> px >> py >> pz >> E >> pid >> charge;
      particle.reset(px,py,pz,E);
    }

    // push event onto back of full_event vector                                                                                                                                      

    if ( nsub <= 1 ) {
      if (hard_event_charged && fabs(charge)>0.99) hard_event_charged->push_back(particle);
    }
    else {
      if (pileup_event_charged && fabs(charge)>0.99) pileup_event_charged->push_back(particle);
    }

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
