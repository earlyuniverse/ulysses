//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2021, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#include "fastjet/PxConePlugin.hh"

#include "fastjet/ClusterSequence.hh"
#include <sstream>

// pxcone stuff
#include "pxcone.h"


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

thread_safety_helpers::FirstTimeTrue PxConePlugin::_first_time;

string PxConePlugin::description () const {
  ostringstream desc;
  
  desc << "PxCone jet algorithm with " 
       << "cone_radius = "        << cone_radius        () << ", "
       << "min_jet_energy = "     << min_jet_energy     () << ", "
       << "overlap_threshold  = " << overlap_threshold  () << ", "
       << "E_scheme_jets  = "     << E_scheme_jets      () << ", "
       << "mode (1=e+e-, 2=hh) = " << _mode
       << " (NB: non-standard version of PxCone, containing small bug fixes by Gavin Salam)";

  return desc.str();
}


void PxConePlugin::run_clustering(ClusterSequence & clust_seq) const {
  // print a banner if we run this for the first time
  //_print_banner(clust_seq.fastjet_banner_stream());
 
  // only have hh mode
  //int mode = 2;

  int    ntrak = clust_seq.jets().size(), itkdm = 4;
  double *ptrak = new double[ntrak*4+1];
  for (int i = 0; i < ntrak; i++) {
    ptrak[4*i+0] = clust_seq.jets()[i].px();
    ptrak[4*i+1] = clust_seq.jets()[i].py();
    ptrak[4*i+2] = clust_seq.jets()[i].pz();
    ptrak[4*i+3] = clust_seq.jets()[i].E();
  }  

  // max number of allowed jets
  int mxjet = ntrak;
  int njet;
  double *pjet  = new double[mxjet*5+1];
  int    *ipass = new int[ntrak+1];
  int    *ijmul = new int[mxjet+1];
  int ierr;

  // run pxcone
  pxcone(
    _mode  ,    // 1=>e+e-, 2=>hadron-hadron
    ntrak  ,    // Number of particles
    itkdm  ,    // First dimension of PTRAK array: 
    ptrak  ,    // Array of particle 4-momenta (Px,Py,Pz,E)
    cone_radius()  ,    // Cone size (half angle) in radians
    min_jet_energy() ,    // Minimum Jet energy (GeV)
    overlap_threshold()  ,    // Maximum fraction of overlap energy in a jet
    mxjet  ,    // Maximum possible number of jets
    njet   ,    // Number of jets found
    pjet ,       // 5-vectors of jets
    ipass,      // Particle k belongs to jet number IPASS(k)-1
                // IPASS = -1 if not assosciated to a jet
    ijmul,      // Jet i contains IJMUL[i] particles
    ierr        // = 0 if all is OK ;   = -1 otherwise
    );

  if (ierr != 0) throw Error("An error occurred while running PXCONE");

  // now transfer information back 
  valarray<int> last_index_created(njet);

  vector<vector<int> > jet_particle_content(njet);

  // get a list of particles in each jet
  for (int itrak = 0; itrak < ntrak; itrak++) {
    int jet_i = ipass[itrak] - 1;
    if (jet_i >= 0) jet_particle_content[jet_i].push_back(itrak);
  }

  // now transfer the jets back into our own structure -- we will
  // mimic the cone code with a sequential recombination sequence in
  // which the jets are built up by adding one particle at a time
  for(int ipxjet = njet-1; ipxjet >= 0; ipxjet--) {
    const vector<int> & jet_trak_list = jet_particle_content[ipxjet];
    int jet_k = jet_trak_list[0];
  
    for (unsigned ilist = 1; ilist < jet_trak_list.size(); ilist++) {
      int jet_i = jet_k;
      // retrieve our misappropriated index for the jet
      int jet_j = jet_trak_list[ilist];
      // do a fake recombination step with dij=0
      double dij = 0.0;
      //clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
      if (ilist != jet_trak_list.size()-1 || E_scheme_jets()) {
        // our E-scheme recombination in cases where it doesn't matter
        clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, jet_k);
      } else {
        // put in pxcone's momentum for the last recombination so that the
        // final inclusive jet corresponds exactly to PXCONE's
        clust_seq.plugin_record_ij_recombination(jet_i, jet_j, dij, 
                                PseudoJet(pjet[5*ipxjet+0],pjet[5*ipxjet+1],
                                          pjet[5*ipxjet+2],pjet[5*ipxjet+3]),
                                                 jet_k);
      }
    }
  
    // NB: put a sensible looking d_iB just to be nice...
    double d_iB = clust_seq.jets()[jet_k].perp2();
    clust_seq.plugin_record_iB_recombination(jet_k, d_iB);
  }


  //// following code is for testing only
  //cout << endl;
  //for (int ijet = 0; ijet < njet; ijet++) {
  //  PseudoJet jet(pjet[ijet][0],pjet[ijet][1],pjet[ijet][2],pjet[ijet][3]);
  //  cout << jet.perp() << " " << jet.rap() << endl;
  //}
  //cout << "-----------------------------------------------------\n";
  //vector<PseudoJet> ourjets(clust_seq.inclusive_jets());
  //for (vector<PseudoJet>::const_iterator ourjet = ourjets.begin();
  //     ourjet != ourjets.end(); ourjet++) {
  //  cout << ourjet->perp() << " " << ourjet->rap() << endl;
  //}
  ////cout << endl;

  delete[] ptrak;
  delete[] ipass;
  delete[] ijmul;
  delete[] pjet;
}

// print a banner for reference to the 3rd-party code
void PxConePlugin::_print_banner(ostream *ostr) const{
  if (! _first_time()) return;

  // make sure the user has not set the banner stream to NULL
  if (!ostr) return;  

  (*ostr) << "#-------------------------------------------------------------------------" << endl;
  (*ostr) << "# You are running the PxCone plugin for FastJet                           " << endl;
  (*ostr) << "# Original code by the Luis Del Pozo, David Ward and Michael H. Seymour   " << endl;
  (*ostr) << "# If you use this plugin, please cite                                     " << endl;
  (*ostr) << "#   M. H. Seymour and C. Tevlin, JHEP 0611 (2006) 052 [hep-ph/0609100].   " << endl;
  (*ostr) << "# in addition to the usual FastJet reference.                             " << endl;
  (*ostr) << "#-------------------------------------------------------------------------" << endl;

  // make sure we really have the output done.
  ostr->flush();
}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
