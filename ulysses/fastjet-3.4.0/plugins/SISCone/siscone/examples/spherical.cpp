///////////////////////////////////////////////////////////////////////////////
// File: spherical.cpp                                                       //
// Description: example program for the CSphsiscone class (spherical SISCone)//
// This file is part of the SISCone project.                                 //
// WARNING: this is not the main SISCone trunk but                           //
//          an adaptation to spherical coordinates                           //
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
// $Revision:: 227                                                          $//
// $Date:: 2008-06-12 20:00:44 -0400 (Thu, 12 Jun 2008)                     $//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "siscone/spherical/momentum.h"
#include "siscone/spherical/siscone.h"

#define R     0.7
#define f     0.5
#define f_alt 0.75

using namespace std;
using namespace siscone_spherical;

int main(){
  vector<CSphmomentum> particles; // list of particles
  CSphsiscone siscone;            // main object for the cone algorithm
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
        particles.push_back(CSphmomentum(px, py, pz, E));
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
  vector<CSphjet>::iterator it_j;
  int i1;
  fprintf(stdout, "#            theta      phi       px         py         pz         E    \n");
  for (it_j = siscone.jets.begin(), i1=0 ; 
       it_j != siscone.jets.end() ; it_j++, i1++){
    fprintf(stdout, "Jet %3d: %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f\n",
	    i1, it_j->v._theta, it_j->v._phi, it_j->v.px, it_j->v.py,  it_j->v.pz,  it_j->v.E);
  }

  return 0;
}
