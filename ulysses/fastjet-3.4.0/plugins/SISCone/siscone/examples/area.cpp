///////////////////////////////////////////////////////////////////////////////
// File: main.cpp                                                            //
// Description: main program that runs siscone from the command line         //
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
// $Revision:: 171                                                          $//
// $Date:: 2007-06-19 10:26:05 -0400 (Tue, 19 Jun 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include "siscone/momentum.h"
#include "siscone/siscone.h"
#include "siscone/area.h"
#include "options.h"

using namespace std;
using namespace siscone;

int main(int argc, char *argv[]){
  vector<Cmomentum> particles;
  Carea siscone_with_area;
  int i,N;
  double px,py,pz,E;
  Coptions opts;
  char fline[512];

  if (opts.parse_options(argc, argv))
    exit(1);

  // deal with help message
  if (opts.help_flag){
    opts.print_help();
    exit(0);
  }

  // deal with version flag
  if (opts.version_flag){
    opts.print_version();
    exit(0);
  }

  // various files used to read input data and store results
  FILE *flux;
  FILE *fpart;

  // read particles
  if (opts.verbose_flag) cout << "reading particles" << endl;
  flux = fopen(opts.ev_name, "r");
  if (flux==NULL){
    cerr << "cannot read event" << endl;
    return 1;
  }

  N=0;
  fpart = fopen("particles.dat", "w+");
  while ((opts.N_stop!=0) && (fgets(fline, 512, flux)!=NULL)){
    if (fline[0]!='#'){ // skip lines beginning with '#'
      if (sscanf(fline, "%le%le%le%le", &px, &py, &pz, &E)==4){    
	particles.push_back(Cmomentum(px, py, pz, E));
	fprintf(fpart, "%e\t%e\n",   particles[N].eta, particles[N].phi);
	N++;
	opts.N_stop--;
      } else {
	cout << "error in reading event file Giving up." << endl;
	fclose(flux);
	fclose(fpart);
	exit(2);
      }
    }
  }
  fclose(flux);
  fclose(fpart);
  if (opts.verbose_flag) 
    cout << "  working with " << N << " particles" << endl;

  // compute jets
  if (opts.verbose_flag) cout << "computing jet contents" << endl;
  i=siscone_with_area.compute_areas(particles, opts.R, opts.f, opts.npass, opts.SM_var);
  if (opts.verbose_flag){
    unsigned int pass;
    for (pass=0;pass<siscone_with_area.protocones_list.size();pass++)
      cout << "    pass " << pass << " found " << siscone_with_area.protocones_list[pass].size()
	   << " stable cones" << endl;
    cout << "  Final result: " << i << " jets found" << endl;
  }

  // save jets
  if (opts.verbose_flag) 
    cout << "saving result" << endl;
  flux = fopen("jets_with_area.dat", "w+");
  vector<Cjet_area>::iterator ja;
  for (ja=siscone_with_area.jet_areas.begin();ja!=siscone_with_area.jet_areas.end();ja++){
    fprintf(flux, "%e\t%e\t%e\t%e\t%e\n",
            ja->v.perp(), ja->v.eta, ja->v.phi,
            ja->active_area, ja->passive_area);
  }

  fclose(flux);

  if (opts.verbose_flag) 
    cout << "bye..." << endl;

  return 0;
}
