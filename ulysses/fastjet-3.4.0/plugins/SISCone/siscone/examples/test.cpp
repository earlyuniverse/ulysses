///////////////////////////////////////////////////////////////////////////////
// File: test.cpp                                                            //
// Description: example program that implements tests with random particles  //
//              and output various informations                              //
//                                                                           //
// Note: for a usage example of SISCone, we advise looking at main.cpp       //
//       or http://projects.hepforge.org/siscone/usage.html                  //
//                                                                           //
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
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <math.h>

#include "siscone/momentum.h"
#include "siscone/siscone.h"

#define N_default  500
#define R          0.7
#define F          0.75

using namespace std;
using namespace siscone;

int main(int argc, char* argv[]){
  vector<Cmomentum> particles;
  Cmomentum *v;
  double phi=0, eta=0, pt=1;
  unsigned int N;

  unsigned int i;
  FILE *flux;

  if (argc==1){
    //cout << "using default number of particles" << endl;
    N = N_default;
  } else {
    sscanf(argv[1], "%u", &N);
    //cout << "using " << N << " particles" << endl;
  }    

  // Initialise random number generator
  timeval timestamp;
  gettimeofday(&timestamp, NULL);
  srand(timestamp.tv_usec);

  // build particle list
  cout << "build particle list" << endl;
  flux = fopen("particles.dat", "w+");
  for (i=0;i<N;i++){
    // uniform eta between -5 and 5
    eta = -5.0+10.0*rand()/(RAND_MAX+1.0);

    // uniform azimuth
    phi = 2.0*M_PI*rand()/(RAND_MAX+1.0);

    // logarithmically uniform pt (between 1e-3 and 100 GeV)
    pt = exp(log(0.001)+log(1e5)*rand()/(RAND_MAX+1.0));

    particles.push_back(Cmomentum(pt*cos(phi), pt*sin(phi), pt*sinh(eta), pt*cosh(eta)));

    fprintf(flux, "%e\t%e\t%e\n",   particles[i].eta, particles[i].phi,particles[i].perp());
  }
  fclose(flux);

  cout << "SISCone: initialise engine" << endl;
  Csiscone siscone;

  // cluster the event
  cout << "cluster the event" << endl;
  siscone.compute_jets(particles, R, F);

#ifdef DEBUG_STABLE_CONES 
  cout << "hash_candidates=" << siscone.nb_hash_cones_total << " in " << siscone.nb_hash_occupied_total << " cells" << endl;
#endif
  // save list of stable cones
  cout << "save stable cone results:" << endl;
  unsigned int pass;
  flux = fopen("protocones.dat", "w+");
  for (pass=0;pass<siscone.protocones_list.size();pass++){
    cout << "    pass " << pass << " found " << siscone.protocones_list[pass].size()
  	 << " stable cones" << endl;
    fprintf(flux, "# pass %d: %u stable cones\n", pass, 
	    (unsigned int) siscone.protocones_list[pass].size());
    for (i=0;i<siscone.protocones_list[pass].size();i++){
      v = &(siscone.protocones_list[pass][i]);
      fprintf(flux, "%e\t%e\t%e\n", v->eta, v->phi, v->perp());
    }
  }
  fclose(flux);
  
  cout << "bye..." << endl;

  return 0;
}
