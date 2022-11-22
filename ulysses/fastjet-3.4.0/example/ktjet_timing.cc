//----------------------------------------------------------------------
// ktjet_timing.cc: Program similar to fastjet_timing.cc, to be used
//                  for timing and result comparisons. For usage, see
//                  explanations given at the start of fastjet_timing.cc
// 
// Note: - some options present in fastjet_timing are missing
//       - one or two behave differently, notably -write
//
#include<iostream>
#include<sstream>
#include<valarray>
#include<vector>
#include<cstddef> // for size_t
#include "CmdLine.hh"
#include "fastjet/internal/numconsts.hh"

#include <cstdio> // recent g++ compiler are a bit more picky

/** Need to include these KtJet Headers */
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"
using namespace std;
using namespace KtJet;

inline double pow2(const double x) {return x*x;}

/// a program to test and time the kt algorithm as implemented in ktjet
int main (int argc, char ** argv) {

  CmdLine cmdline(argc,argv);
  //bool clever = !cmdline.present("-dumb");
  int  repeat = cmdline.int_val("-repeat",1);
  int  combine = cmdline.int_val("-combine",1);
  bool write   = cmdline.present("-write");
  double ktR   = cmdline.double_val("-r",1.0);
  double inclkt = cmdline.double_val("-incl",-1.0);
  int    excln  = cmdline.int_val   ("-excln",-1);
  double excld  = cmdline.double_val("-excld",-1.0);
  int  nev     = cmdline.int_val("-nev",1);
  bool   massless = cmdline.present("-massless");
  bool   get_all_dij   = cmdline.present("-get-all-dij");


  for (int iev = 0; iev < nev; iev++) {
  vector<KtJet::KtLorentzVector> jets;
  string line;
  int  ndone = 0;
  while (getline(cin, line)) {
      //cout << line<<endl;
    istringstream linestream(line);
    if (line == "#END") {
      ndone += 1;
      if (ndone == combine) {break;}
    }
    if (line.substr(0,1) == "#") {continue;}
    valarray<double> fourvec(4);
    linestream >> fourvec[0] >> fourvec[1] >> fourvec[2] >> fourvec[3];
    if (massless) {
      linestream >> fourvec[0] >> fourvec[1] >> fourvec[2];
      fourvec[3] = sqrt(pow2(fourvec[0])+pow2(fourvec[1])+pow2(fourvec[2]));}
    else {
      linestream >> fourvec[0] >> fourvec[1] >> fourvec[2] >> fourvec[3];
    }
    KtJet::KtLorentzVector p(fourvec[0],fourvec[1],fourvec[2],fourvec[3]);
    jets.push_back(p);
  }
  
  // set KtEvent flags
  int type, angle, recom;
  ostringstream info;
  if (cmdline.present("-eekt")) {
    type  = 1; // e+e-
    angle = 1; // angular
    recom = 1; // E
    info << "Algorithm: KtJet e+e- kt algorithm" ;
  } else {
    type  = 4; // PP
    angle = 2; // delta R
    recom = 1; // E
    info << "Algorithm: KtJet (long.inv.) with R = " << ktR ;
  }
  //double rparameter = 1.0;

  for (int i = 0; i < repeat ; i++) {
    // Construct the KtEvent object 
    KtJet::KtEvent ev(jets,type,angle,recom,ktR);

    if (i!=0) {continue;}
    int nparticles = jets.size();
    cout << "Number of particles = "<< nparticles << endl;
    cout << info.str() << endl;

    // Print out the number of final state jets
    //std::cout << "Number of final state jets: " << ev.getNJets() << std::endl;
    if (write) {
      /** Retrieve the final state jets from KtEvent sorted by Pt*/
      std::vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
      
      /** Print out jets 4-momentum and Pt */
      std::vector<KtJet::KtLorentzVector>::const_iterator itr = jets.begin();
      for( ; itr != jets.end() ; ++itr) {
	std::cout << "Jets Pt2: " << pow2((*itr).perp()) << std::endl; 
      }
    }

    if (inclkt >= 0.0) {
      // Retrieve the final state jets from KtEvent sorted by Pt
      std::vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
    
      // Print out index, rap, phi, pt 
      for (size_t j = 0; j < jets.size(); j++) {
	if (jets[j].perp() < inclkt) {break;}
	double phi = jets[j].phi();
	if (phi < 0.0) {phi += fastjet::twopi;}
	printf("%5u %15.8f %15.8f %15.8f\n",j,jets[j].rapidity(),phi,jets[j].perp());
      }
    }

    if (excln > 0) {
      ev.findJetsN(excln);
      vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
      cout << "Printing "<<excln<<" exclusive jets\n";
      for (size_t j = 0; j < jets.size(); j++) {
        double phi = jets[j].phi();
        if (phi < 0) phi += fastjet::twopi;
	printf("%5u %15.8f %15.8f %15.8f\n",j,
	       jets[j].rapidity(),phi,jets[j].perp());
      }
    }

    if (excld > 0.0) {
      ev.findJetsD(excld);
      vector<KtJet::KtLorentzVector> jets = ev.getJetsPt();
      cout << "Printing exclusive jets for d = "<<excld<<"\n";
      for (size_t j = 0; j < jets.size(); j++) {
        double phi = jets[j].phi();
        if (phi < 0) phi += fastjet::twopi;
	printf("%5u %15.8f %15.8f %15.8f\n",j,
	       jets[j].rapidity(),phi,jets[j].perp());
      }
    }

    if (get_all_dij) {
      for (int i = nparticles-1; i > 0; i--) {
        printf("d for n = %4d -> %4d is %14.5e\n", i+1, i, ev.getDMerge(i));
      }
    }

  }

}
}
