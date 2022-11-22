///
/// Use this file from root. First compile it and then run showjets(...)
/// on a file produced by ClusterSequence::print_jets_for_root(...)
///
/// root> .L jet-plots.C+
/// root> showjets(filename [, label])
///
#include<vector>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cstdlib>
#include "TLatex.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH2.h"
#include "THStack.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TColor.h"

using namespace std;


class JetHist {
private:
  vector<TH2D *> _jets;
  TH2D * _background;
  string _comment;
public:
  static double default_etamax;
  static int    default_nbins;

  JetHist(const string & filename, double etamax=default_etamax, int nbins=default_nbins);
  ~JetHist();
  string comment() {return _comment;}
  THStack stack;
  TH2D * jet(int i) {return i>= 0 ? _jets[i] : _background;}

};

double JetHist::default_etamax = 6.0;
int    JetHist::default_nbins  = 40; // y: 2*nbins; phi: nbins

// get jet "histograms" from filename which is expected to be made of repeated
// blocks as follows:
//       jet# eta phi pt ...
//        ipart eta phi pt
//        ipart eta phi pt
//        ...
//       #END
JetHist::JetHist (const string & filename, double etamax, int nbins) {
  ifstream file(filename.c_str());
  string line;
  //double etamax=6;
  //double etamax=5;
  double phimax = 2*3.14159265;
  //int    nbins=40;

  // construct a histogram for the background to the jets
  ostringstream bname;
  bname << filename <<"-background";
  _background = new TH2D(bname.str().c_str(),bname.str().c_str(),
  			 2*nbins,-etamax,etamax,nbins,0.0,phimax);
  //_background = new TH2D(bname.str().c_str(),bname.str().c_str(),
  //			 2*nbins,-etamax,etamax,2,0.0,phimax);
  _background->SetFillColor(kWhite);
  // these were supposed to have labelled the axes, but it doesn't work.
  _background->GetXaxis()->SetTitle("#eta");
  _background->GetYaxis()->SetTitle("#phi");
  _background->GetZaxis()->SetTitle("p_{#perp}");
  stack.Add(_background);

  while (getline(file,line)) {
    if (line.substr(0,1) != " ") { // all interesting lines start with space?
      // extract a comment if there is one
      if (line.substr(0,2) == "# ") {
	_comment = line.substr(2,line.length()-2);
      }
      continue;
    } 
    ostringstream name;
    name << filename<<"-jet-"<< _jets.size();
    TH2D * hist = new TH2D(name.str().c_str(),name.str().c_str(),
			   2*nbins,-etamax,etamax,nbins,0.0,phimax);
    int    i;
    double eta, phi, pt;
    //cout << filename <<": jet "<<_jets.size()<<endl;
    bool have_line = true;
    while (have_line || getline(file,line)) {
      have_line = false;
      if (line.substr(0,4) == "#END") {break;}
      istringstream sline(line);
      sline >> i >> eta >> phi >> pt;
      //cout << i << " "<<eta<<" "<<phi<<" "<<pt<<endl;
      hist->Fill(eta,phi,pt); // fill at phi,eta with weight pt
      
      // workaround for bug in stacks: fill all lower elements of the stack
      // with a fake amount -- this, miraculously will lead to correct coloring
      // of the top of the stack!
      //for (unsigned int j = 0; j < _jets.size(); j++) {
      //  _jets[j]->Fill(phi,eta,1e-7); 
      //}
    }
    // give it a colour (whatever that means...)
    //hist->SetFillColor(_jets.size());
    int njet = _jets.size();
    //hist->SetFillColor(njet+2);
    hist->SetFillColor(njet%50+2); // %50 seems tomake to diff to many-jet case
    //if (njet == 0) hist->SetFillColor(kRed);
    //else if (njet == 1) hist->SetFillColor(kBlue);
    //else  hist->SetFillColor(kGreen);

    // add it to the list of jets (so we can delete it later...)
    _jets.push_back(hist);

    // put it onto the stack
    stack.Add(hist);
    //if (njet == 2) break;
    //break;
  }
}

// clean up --------------
JetHist::~JetHist () {
  for (unsigned int i = 0; i < _jets.size(); i++) {
    delete _jets[i];
  }
  delete _background;
}

//----------------------------------------------------------------------
/// set up a reasonable bunch of colours
void set_default_colours(TCanvas * lego) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(1);
  gStyle->SetFrameFillColor(0);

  Int_t cancolor = 0;
  lego->SetFillColor(cancolor);

  int ngrey = 3;
  for (int ir = 0; ir < ngrey; ir++) {
    for (int ig = 0; ig < ngrey; ig++) {
      for (int ib = 0; ib < ngrey; ib++) {
        int icol = 7+ir + ngrey *ig + ngrey*ngrey * ib;
        TColor * color=(TColor*)(gROOT->GetListOfColors()->At(icol));
        if (icol == 7) {
          // avoid white -- put grey instead
          color->SetRGB(0.5,0.5,0.5);
        } else {
          color->SetRGB(1-ir*1.0/ngrey,1-ig*1.0/ngrey,1-ib*1.0/ngrey);
        }
      }
    }
  }
}

//----------------------------------------------------------------------
/// show the jets contained in filename (as produced by
/// ClusterSequence::print_jets_for_root()), with an optional label
TCanvas * showjets (const char * filename, const char * label = 0) {

  // display the various 2-d drawing options
  gROOT->Reset();

  // set up canvas
  TCanvas * lego = new TCanvas("lego","lego options",400,50,800,600);
  lego->SetTheta(30.0);
  lego->SetPhi(20.0);

  // orientation used for plots in subtraction paper
  //lego->SetTheta(62.15);
  //lego->SetPhi(9.15);

  ////vector<double> col 

  set_default_colours(lego);

  TPaveLabel pl;
   
  JetHist * jets = new JetHist(filename);
  jets->stack.Draw("lego1"); // cyl does not work with 5.16
  if (label != 0) {
    Float_t x1=0.63, y1=0.875, x2=0.95, y2=0.925;
    pl.DrawPaveLabel(x1,y1,x2,y2,label,"brNDC");
  } else if (jets->comment() != "") {
    Float_t x1=0.15, y1=0.875, x2=0.95, y2=0.925;
    pl.DrawPaveLabel(x1,y1,x2,y2,jets->comment().c_str(),"brNDC");
  }

  // normal histogram labels not working, so draw them by hand
  TLatex l;
  l.SetTextAlign(22);
  l.SetTextSize(0.05);
  //l.DrawLatex(0.0,0.85,"anti-k_{t}, R=1");

  l.SetTextSize(0.04);
  l.DrawLatex(0.20,-0.98,"y");
  l.SetTextAlign(32);
  l.DrawLatex(-0.7,0.8,"p_{t} [GeV]");
  l.DrawLatex(-0.6,-0.78,"#phi");

  // do not delete jets -- otherwise you lose everything!;

  return lego;
  ///
  lego->Update();
}
