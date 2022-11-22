//----------------------------------------------------------------------
// ktjet example program that should do the same thing as the
// fastjet_example program (as of 3 Feb 2006)
//
// NB: the ../Makefile may need to be modified to set proper
//     paths for access to the CLHEP and KtJet libraries.
//----------------------------------------------------------------------
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 

/** Need to include these KtJet Headers */
#include "KtJet/KtEvent.h"
#include "KtJet/KtLorentzVector.h"

using namespace std;
using namespace KtJet;

// a declaration of a function that pretty prints a list of jets
void print_jets (const vector<KtLorentzVector> &);

/// an example program showing how the fastjet_example program would
/// be translated for use with ktjet.
int main (int argc, char ** argv) {
  
  vector<KtLorentzVector> input_particles;
  
  // read in input particles
  double px, py , pz, E;
  while (cin >> px >> py >> pz >> E) {
    // create a KtLorentzVector with these components and put it onto
    // back of the input_particles vector
    input_particles.push_back(KtLorentzVector(px,py,pz,E)); 
  }

  // run the inclusive jet clustering in PP mode using the covariant
  // E-scheme for recobination (type=4, angle=2, recom=1, rparameter=1.0)
  double Rparam = 1.0;
  KtEvent clust_seq(input_particles,4,2,1,Rparam);

  // extract the inclusive jets with pt > 5 GeV, sorted by pt
  double ptmin = 5.0;
  vector<KtLorentzVector> temporary_jets = clust_seq.getJetsPt();
  vector<KtLorentzVector> inclusive_jets;
  for (unsigned int i = 0; i < temporary_jets.size(); i++) {
    if (temporary_jets[i].perp() >= ptmin) {
      inclusive_jets.push_back(temporary_jets[i]);}
    else {break;}
  }

  // print them out
  cout << "Printing inclusive jets with pt > "<< ptmin<<" GeV\n";
  cout << "---------------------------------------\n";
  print_jets(inclusive_jets);
  cout << endl;

  // Extract the exclusive jets with dcut = 25 GeV^2.
  double dcut = 25.0; 
  // Note that KtJet's definition of dij differs from Ellis&Soper (and
  // fastjet) in the case where Rparam /= 1.0 (though in this case one
  // should perhaps not be using the exclusive kt algorithm in any case).
  clust_seq.findJetsD(dcut * Rparam*Rparam);
  vector<KtLorentzVector> exclusive_jets = clust_seq.getJetsPt();

  // print them out
  cout << "Printing exclusive jets with dcut = "<< dcut<<" GeV^2\n";
  cout << "--------------------------------------------\n";
  print_jets(exclusive_jets);


}


//----------------------------------------------------------------------
// a function that pretty prints a list of jets
void print_jets (const vector<KtLorentzVector> & jets) {

  // label the columns
  printf("%5s %15s %15s %15s %15s\n","jet #", "rapidity", 
	 "phi", "pt", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < jets.size(); i++) {
    int n_constituents = jets[i].getConstituents().size();
    double phi = jets[i].phi();
    if (phi < 0.0) {phi += 6.283185307179586476925286766559005768394;}
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, jets[i].rapidity(), phi,
	   jets[i].perp(), n_constituents);
  }

}
