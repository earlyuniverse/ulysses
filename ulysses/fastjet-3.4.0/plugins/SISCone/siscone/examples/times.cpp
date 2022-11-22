///////////////////////////////////////////////////////////////////////////////
// File: times.cpp                                                           //
// Description: example program that computes execution times                //
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

#define Nruns 32
#define R    0.7
#define f    0.5

using namespace std;
using namespace siscone;

timeval time_start, time_end;

// compute time spent between time_start and time_end
int time_spent(){
  timeval time_diff;
  
  // compute different with initial time
  time_diff.tv_sec = time_end.tv_sec-time_start.tv_sec;
  if (time_end.tv_usec > time_start.tv_usec){
    time_diff.tv_usec = time_end.tv_usec-time_start.tv_usec;
  } else {
    time_diff.tv_sec--;
    time_diff.tv_usec = (1000000+time_end.tv_usec)-time_start.tv_usec;
  }
  
  return 1000000*time_diff.tv_sec+time_diff.tv_usec;
}



int main(){
  vector<Cmomentum> particles;
  Csiscone siscone;
  double eta,phi;

  // number of events and particles
  int i, N;
  int n_ev, part_inc;

  // time statistics variables
  int time_siscone;

  // save files
  FILE *flux;

  // initialise random number generator
  cout << "initialise random number generator" << endl;
  timeval timestamp;

  gettimeofday(&timestamp, NULL);
  srand(timestamp.tv_usec);

  flux = fopen("times.dat", "w+");

  N = 1;
  part_inc = 1;
  do{
    fprintf(stdout, "\r%5d particles\n", N);
    time_siscone=0;

    for (n_ev=0;n_ev<Nruns;n_ev++){
      // build particle list
      particles.clear();
      for (i=0;i<N;i++){
	eta = -3.0+6.0*rand()/(RAND_MAX+1.0);
	phi = 2.0*M_PI*rand()/(RAND_MAX+1.0);
	particles.push_back(Cmomentum(cos(phi), sin(phi), tanh(eta), 1.0));
      }
      
      // run siscone
      gettimeofday(&time_start, NULL);
      siscone.compute_jets(particles, R, f);
      gettimeofday(&time_end, NULL);
      time_siscone+=time_spent();
    }

    fprintf(flux, "%d\t%e\n", N, time_siscone/(1.0*Nruns));

    N+=part_inc;
    if (N==(part_inc<<3))
      part_inc <<= 1;
    //  } while (N<=1024);
  } while (N<=1024);
  
  fclose(flux);
  fprintf(stdout, "\n");

  cout << "bye..." << endl;

  return 0;
}
