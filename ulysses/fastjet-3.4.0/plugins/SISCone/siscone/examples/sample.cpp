///////////////////////////////////////////////////////////////////////////////
// File: sample.cpp                                                          //
// Description: example program for the Csiscone class (see documentation)   //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision::                                                              $//
// $Date::                                                                  $//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "siscone/momentum.h"
#include "siscone/siscone.h"

#define R     0.7
#define f     0.5
#define f_alt 0.75

using namespace std;
using namespace siscone;

int main(){
  vector<Cmomentum> particles;    // list of particles
  Csiscone siscone;               // main object for the cone algorithm
  int i;                          // loop index
  int N;                          // number of particles
  double px,py,pz,E;              // particles 4-momentum
  char fline[512];                // line to read from a file

  // read particles
  FILE *flux;
  flux = fopen("events/single-event.dat", "r");
  if (flux==NULL){
    cerr << "cannot read event" << endl;
    return 1;
  }

  N=0;
  while (fgets(fline, 512, flux)!=NULL){
    if (fline[0]!='#'){ // skip lines beginning with '#'
      if (sscanf(fline, "%le%le%le%le", &px, &py, &pz, &E)==4){
        particles.push_back(Cmomentum(px, py, pz, E));
        N++;
      } else {
        cout << "error in reading event file Giving up." << endl;
        fclose(flux);
        return 2;
      }
    }
  }
  fclose(flux);

  // compute jets
  // first compute with multiple passes (default)
  i=siscone.compute_jets(particles, R, f);
  cout << "  " << i << " jets found in multi-pass run" << endl;

  // then, recompute it with a different f
  i=siscone.recompute_jets(f_alt);
  cout << "  " << i << " jets found with alternative f" << endl;

  // one pass
  i=siscone.compute_jets(particles, R, f, 1);
  cout << "  " << i << " jets found in single-pass run" << endl;

  // show jets
  vector<Cjet>::iterator it_j;
  int i1;
  fprintf(stdout, "#             pT        eta      phi       px         py         pz         E    \n");
  for (it_j = siscone.jets.begin(), i1=0 ; 
       it_j != siscone.jets.end() ; it_j++, i1++){
    fprintf(stdout, "Jet %3d: %10.3f %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	    i1, it_j->v.perp(), it_j->v.eta, it_j->v.phi, it_j->v.px, it_j->v.py,  it_j->v.pz,  it_j->v.E);
  }

  return 0;
}
