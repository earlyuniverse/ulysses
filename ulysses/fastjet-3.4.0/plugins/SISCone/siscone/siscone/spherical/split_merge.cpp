///////////////////////////////////////////////////////////////////////////////
// File: split_merge.cpp                                                     //
// Description: source file for splitting/merging (contains the CJet class)  //
// This file is part of the SISCone project.                                 //
// WARNING: this is not the main SISCone trunk but                           //
//          an adaptation to spherical coordinates                           //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006-2008 Gavin Salam and Gregory Soyez                     //
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

#include <siscone/siscone_error.h>
#include "split_merge.h"
#include "momentum.h"
#include <limits>   // for max
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cassert>
#include <cmath>

namespace siscone_spherical{

using namespace std;

/********************************************************
 * class CSphjet implementation                         *
 * real Jet information.                                *
 * This class contains information for one single jet.  *
 * That is, first, its momentum carrying information    *
 * about its centre and pT, and second, its particle    *
 * contents                                             *
 ********************************************************/
// default ctor
//--------------
CSphjet::CSphjet(){
  n = 0;
  v = CSphmomentum();
  E_tilde = 0.0;
  sm_var2 = 0.0;
  pass = CJET_INEXISTENT_PASS; // initialised to a value that should
                               // notappear in the end (after clustering)
}

// default dtor
//--------------
CSphjet::~CSphjet(){

}

// ordering of jets in E (e.g. used in final jets ordering)
//----------------------------------------------------------
bool jets_E_less(const CSphjet &j1, const CSphjet &j2){
  return j1.v.E > j2.v.E;
}


/********************************************************
 * CSphsplit_merge_ptcomparison implementation          *
 * This deals with the ordering of the jets candidates  *
 ********************************************************/

// odering of two jets
// The variable of the ordering is pt or mt 
// depending on 'split_merge_scale' choice
//
// with EPSILON_SPLITMERGE defined, this attempts to identify
// delicate cases where two jets have identical momenta except for
// softish particles -- the difference of pt's may not be correctly
// identified normally and so lead to problems for the fate of the
// softish particle.
//
// NB: there is a potential issue in momentum-conserving events,
// whereby the harder of two jets becomes ill-defined when a soft
// particle is emitted --- this may have a knock-on effect on
// subsequent split-merge steps in cases with sufficiently large R
// (but we don't know what the limit is...)
//------------------------------------------------------------------
bool CSphsplit_merge_ptcomparison::operator ()(const CSphjet &jet1, const CSphjet &jet2) const{
  double q1, q2;

  // compute the value for comparison for both jets
  // This depends on the choice of variable (mt is the default)
  q1 = jet1.sm_var2;
  q2 = jet2.sm_var2;

  bool res = q1 > q2;

  // if we enable the refined version of the comparison (see defines.h),
  // we compute the difference more precisely when the two jets are very
  // close in the ordering variable.
#ifdef EPSILON_SPLITMERGE
  if ( (fabs(q1-q2) < EPSILON_SPLITMERGE*max(q1,q2)) &&
       (jet1.v.ref != jet2.v.ref) ) {
#ifdef DEBUG_SPLIT_MERGE
    cout << "Using high-precision ordering tests" << endl;
#endif
    // get the momentum of the difference
    CSphmomentum difference;
    double E_tilde_difference;
    get_difference(jet1,jet2,&difference,&E_tilde_difference);
    
    // use the following relation: pt1^2 - pt2^2 = (pt1+pt2)*(pt1-pt2)
    double qdiff;
    CSphmomentum sum = jet1.v ;
    sum +=  jet2.v;
    double E_tilde_sum = jet1.E_tilde + jet2.E_tilde;
    
    // depending on the choice of ordering variable, set the result
    switch (split_merge_scale){
    case SM_Etilde:  
      qdiff = E_tilde_sum*E_tilde_difference;
      break;
    case SM_E:  
      qdiff = sum.E*difference.E;
      break;
    default:
      throw siscone::Csiscone_error("Unsupported split-merge scale choice: "
				    + SM_scale_name());
    }
    res = qdiff > 0;
  }
#endif  // EPSILON_SPLITMERGE

  return res;
}


/// return a name for the sm scale choice 
/// NB: used internally and also by fastjet
std::string split_merge_scale_name(Esplit_merge_scale sms) {
  switch(sms) {
  case SM_E:
    return "E (IR unsafe for pairs of identical decayed heavy particles)";
  case SM_Etilde:
    return "Etilde (sum of E.[1+sin^2(theta_{i,jet})])";
  default:
    return "[SM scale without a name]";
  }
}


// get the difference between 2 jets
//  - j1         first jet
//  - j2         second jet
//  - v          jet1-jet2
//  - pt_tilde   jet1-jet2 pt_tilde
// return true if overlapping, false if disjoint
//-----------------------------------------------
void CSphsplit_merge_ptcomparison::get_difference(const CSphjet &j1, const CSphjet &j2, 
						  CSphmomentum *v, double *E_tilde) const {
  int i1,i2;

  // initialise
  i1=i2=0;
  *v = CSphmomentum();
  *E_tilde = 0.0;

  CSph3vector jet1_axis = j1.v;
  //jet1_axis /= j1.v._norm;
  jet1_axis /= j1.v.E;
  CSph3vector jet2_axis = j2.v;
  //jet2_axis /= j2.v._norm;
  jet2_axis /= j2.v.E;

  // compute overlap
  // at the same time, we store union in indices
  // note tat for Etilde, we'll add the direct energy contributino at the end
  do{
    if (j1.contents[i1]==j2.contents[i2]) {
      const CSphmomentum & p = (*particles)[j1.contents[i1]];      
      (*E_tilde) += p.E*((norm2_cross_product3(p,jet1_axis)-norm2_cross_product3(p,jet2_axis))/(*particles_norm2)[j1.contents[i1]]);
      i1++;
      i2++;
    } else if (j1.contents[i1]<j2.contents[i2]){
      const CSphmomentum & p = (*particles)[j1.contents[i1]];      
      (*v) += p;
      (*E_tilde) += p.E*norm2_cross_product3(p,jet1_axis)/(*particles_norm2)[j1.contents[i1]];
      i1++;
    } else if (j1.contents[i1]>j2.contents[i2]){
      const CSphmomentum &p = (*particles)[j2.contents[i2]];      
      (*v) -= p;
      (*E_tilde) -= p.E*norm2_cross_product3(p,jet2_axis)/(*particles_norm2)[j2.contents[i2]];
      i2++;
    } else {
      throw siscone::Csiscone_error("get_non_overlap reached part it should never have seen...");
    }
  } while ((i1<j1.n) && (i2<j2.n));

  // deal with particles at the end of the list...
  while (i1 < j1.n) {
    const CSphmomentum &p = (*particles)[j1.contents[i1]];      
    (*v) += p;
    (*E_tilde) += p.E*norm2_cross_product3(p,jet1_axis)/(*particles_norm2)[j1.contents[i1]];
    i1++;
  }
  while (i2 < j2.n) {
    const CSphmomentum &p = (*particles)[j2.contents[i2]];      
    (*v) -= p;
    (*E_tilde) -= p.E*norm2_cross_product3(p,jet2_axis)/(*particles_norm2)[j2.contents[i2]];
    i2++;
  }

  // add the direct energy contribution to Etilde
  (*E_tilde) += v->E;
}


/********************************************************
 * class CSphsplit_merge implementation                 *
 * Class used to split and merge jets.                  *
 ********************************************************/
// default ctor
//--------------
CSphsplit_merge::CSphsplit_merge(){
  merge_identical_protocones = false;
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
#ifdef MERGE_IDENTICAL_PROTOCONES_DEFAULT_TRUE
  merge_identical_protocones = true;
#endif
#endif
  _user_scale = NULL;
  indices = NULL;

  // ensure that ptcomparison points to our set of particles (though params not correct)
  ptcomparison.particles = &particles;
  ptcomparison.particles_norm2 = &particles_norm2;
  candidates.reset(new multiset<CSphjet,CSphsplit_merge_ptcomparison>(ptcomparison));

  // no hardest cut (col-unsafe)
  SM_var2_hardest_cut_off = -numeric_limits<double>::max();

  // no energy cutoff for the particles to put in p_uncol_hard
  stable_cone_soft_E2_cutoff = -1.0;

  // no pt-weighted splitting
  use_E_weighted_splitting = false;
}


// default dtor
//--------------
CSphsplit_merge::~CSphsplit_merge(){
  full_clear();
}


// initialisation function
//  - _particles  list of particles
//  - protocones  list of protocones (initial jet candidates)
//  - R2          cone radius (squared)
//  - Emin        minimal energy allowed for jets
//-------------------------------------------------------------
int CSphsplit_merge::init(vector<CSphmomentum> & /*_particles*/, vector<CSphmomentum> *protocones, double R2, double Emin){
  // browse protocones
  return add_protocones(protocones, R2, Emin);
}


// initialisation function for particle list
//  - _particles  list of particles
//-------------------------------------------------------------
int CSphsplit_merge::init_particles(vector<CSphmomentum> &_particles){
  full_clear();

  // compute the list of particles
  // here, we do not need to throw away particles 
  // with infinite rapidity (colinear with the beam)
  particles = _particles;
  n = particles.size();

  // store the particle norm^2
  particles_norm2.resize(n);
  for (int i=0;i<n;i++){
    particles_norm2[i] = particles[i].norm2();
  }

  // ensure that ptcomparison points to our set of particles (though params not correct)
  ptcomparison.particles = &particles;
  ptcomparison.particles_norm2 = &particles_norm2;

  // set up the list of particles left.
  init_pleft();

  indices = new int[n];

  return 0;
}


// build initial list of remaining particles
//------------------------------------------
int CSphsplit_merge::init_pleft(){
  // at this level, we only rule out particles with 
  // infinite rapidity
  // in the parent particle list, index contain the run 
  // at which particles are puts in jets:
  //  - -1 means infinity rapidity
  //  -  0 means not included
  //  -  i mean included at run i
  int i,j;

  // copy particles removing the ones with infinite rapidity
  j=0;
  p_remain.clear();
  for (i=0;i<n;i++){
    // set ref for checkxor
    particles[i].ref.randomize();

    //REMOVED: check if rapidity is not infinite or ill-defined
    //if (fabs(particles[i].pz) < (particles[i].E)){
      p_remain.push_back(particles[i]);
      // set up parent index for tracability
      p_remain[j].parent_index = i;
      // left particles are marked with a 1
      // IMPORTANT NOTE: the meaning of index in p_remain is
      //   somehow different as in the initial particle list.
      //   here, within one pass, we use the index to track whether
      //   a particle is included in the current pass (index set to 0
      //   in add_protocones) or still remain (index still 1)
      p_remain[j].index = 1;

      j++;
      // set up parent-particle index
      particles[i].index = 0;
    //} else {
    //  particles[i].index = -1;
    //}
  }
  n_left = p_remain.size();
  n_pass = 0;

  merge_collinear_and_remove_soft();

  return 0;
}


// partial clearance
// we want to keep   particle list and indices
// for future usage, so do not clear it !
// this is done in full_clear
//----------------------------------------
int CSphsplit_merge::partial_clear(){
  // release jets

  // set up the auto_ptr for the multiset with the _current_ state of
  // ptcomparison (which may have changed since we constructed the
  // class)
  candidates.reset(new multiset<CSphjet,CSphsplit_merge_ptcomparison>(ptcomparison));

  // start off with huge number
  most_ambiguous_split = numeric_limits<double>::max();

  jets.clear();
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if (merge_identical_protocones)
    cand_refs.clear();
#endif

  p_remain.clear();

  return 0;
}


// full clearance
//----------------
int CSphsplit_merge::full_clear(){
  partial_clear();

  // clear previously allocated memory
  if (indices != NULL){
    delete[] indices;
  }
  particles.clear();

  return 0;
}


// build the list 'p_uncol_hard' from p_remain by clustering collinear particles
// note that thins in only used for stable-cone detection 
// so the parent_index field is unnecessary
//-------------------------------------------------------------------------
int CSphsplit_merge::merge_collinear_and_remove_soft(){
  int i,j;
  vector<CSphmomentum> p_sorted;
  bool collinear;
  double dphi;

  p_uncol_hard.clear();

  // we first sort the particles according to their theta angle
  for (i=0;i<n_left;i++)
    p_sorted.push_back(p_remain[i]);
  sort(p_sorted.begin(), p_sorted.end(), momentum_theta_less);

  // then we cluster particles looping over the particles in the following way
  // if (a particle i has same eta-phi a one after (j))
  // then add momentum i to j
  // else add i to the p_uncol_hard list
  i = 0;
  while (i<n_left){
    // check if the particle passes the stable_cone_soft_E2_cutoff
    if (p_sorted[i].E*p_sorted[i].E<stable_cone_soft_E2_cutoff) {
      i++;
      continue;
    }

    // check if there is(are) particle(s) with the 'same' theta
    collinear = false;
    j=i+1;
    while ((j<n_left) && (fabs(p_sorted[j]._theta-p_sorted[i]._theta)<EPSILON_COLLINEAR) && (!collinear)){
      dphi = fabs(p_sorted[j]._phi-p_sorted[i]._phi);
      if (dphi>M_PI) dphi = twopi-dphi;
      if (dphi<EPSILON_COLLINEAR){
	// i is collinear with j; add the momentum (i) to the momentum (j) 
#ifdef DEBUG_SPLIT_MERGE
	cout << "# collinear merging at point " << p_sorted[i]._theta << ", " << p_sorted[j]._phi << endl;
#endif
	p_sorted[j] += p_sorted[i];
	//p_sorted[j].build_thetaphi();
	p_sorted[j].build_norm();
	// set collinearity test to true
	collinear = true;
      }
      j++;
    }
    // if no collinearity detected, add the particle to our list
    if (!collinear)
      p_uncol_hard.push_back(p_sorted[i]);
    i++;
  }

  return 0;
}


// add a list of protocones
//  - protocones  list of protocones (initial jet candidates)
//  - R2          cone radius (squared)
//  - Emin        minimal energy allowed for jets
//-------------------------------------------------------------
int CSphsplit_merge::add_protocones(vector<CSphmomentum> *protocones, double R2, double Emin){
  int i;
  CSphmomentum *c;
  CSphmomentum *v;
  double tan2R;
  CSphjet jet;

  if (protocones->size()==0)
    return 1;

  E_min = Emin;
  double R = sqrt(R2);
  tan2R = tan(R);
  tan2R *= tan2R;

#ifdef DEBUG_SPLIT_MERGE
    cout << "particle list: ";
    for (int i2=0;i2<n_left;i2++)
      cout << p_remain[i2].parent_index << " " 
	   << p_remain[i2].px << " "  << p_remain[i2].py << " "
	   << p_remain[i2].pz << " "  << p_remain[i2].E  << endl;
    cout << endl;
#endif

  // browse protocones
  // for each of them, build the list of particles in them
  for (vector<CSphmomentum>::iterator p_it = protocones->begin();p_it != protocones->end();p_it++){
    // initialise variables
    c = &(*p_it);

    // browse particles to create cone contents
    // note that jet is always initialised with default values at this level
    jet.v = CSphmomentum();
    jet.contents.clear();
    for (i=0;i<n_left;i++){
      v = &(p_remain[i]);
      if (is_closer(v, c, tan2R)){
	jet.contents.push_back(v->parent_index);
	jet.v+= *v;
	v->index=0;
      }
    }
    jet.n=jet.contents.size();

    // compute Etilde for that jet.
    // we can't do that before as it requires knowledge of the jet axis
    // which has just been computed.
    compute_Etilde(jet);

    // set the momentum in protocones 
    // (it was only known through its spatial coordinates up to now)
    *c = jet.v;
    c->build_thetaphi();

    // set the jet range
    jet.range=CSphtheta_phi_range(c->_theta,c->_phi,R);

#ifdef DEBUG_SPLIT_MERGE
    cout << "adding protojet: ";

    unsigned int phirange=jet.range.phi_range;
    for (unsigned int i2=0;i2<32;i2++) fprintf(stdout, "%d", (phirange&(1<<i2)) >> i2 );
    fprintf(stdout, "\t");
    unsigned int thetarange=jet.range.theta_range;
    for (unsigned int i2=0;i2<32;i2++) fprintf(stdout, "%d", (thetarange&(1<<i2)) >> i2);
    fprintf(stdout, "\t");

    for (int i2=0;i2<jet.n;i2++)
      cout << jet.contents[i2] << " ";
    cout << endl;
#endif

    // add it to the list of jets
    insert(jet);
  }
  
  // update list of included particles
  n_pass++;

#ifdef DEBUG_SPLIT_MERGE
  cout << "remaining particles: "; 
#endif
  int j=0;
  for (i=0;i<n_left;i++){
    if (p_remain[i].index){
      // copy particle
      p_remain[j]=p_remain[i];
      p_remain[j].parent_index = p_remain[i].parent_index;
      p_remain[j].index=1;
      // set run in initial list
      particles[p_remain[j].parent_index].index = n_pass;
#ifdef DEBUG_SPLIT_MERGE
      cout << p_remain[j].parent_index << " ";
#endif
      j++;
    }
  }
#ifdef DEBUG_SPLIT_MERGE
  cout << endl;
#endif
  n_left = j;
  p_remain.resize(j);

  merge_collinear_and_remove_soft();

  return 0;
}

/*
 * remove the hardest protocone and declare it a a jet 
 *  - protocones  list of protocones (initial jet candidates)
 *  - R2          cone radius (squared)
//  - Emin        minimal energy allowed for jets
 * return 0 on success, 1 on error
 *
 * The list of remaining particles (and the uncollinear-hard ones)
 * is updated.
 */
int CSphsplit_merge::add_hardest_protocone_to_jets(std::vector<CSphmomentum> *protocones, double R2, double Emin){

  int i;
  CSphmomentum *c;
  CSphmomentum *v;
  double R, tan2R;
  CSphjet jet, jet_candidate;
  bool found_jet = false;

  if (protocones->size()==0)
    return 1;

  E_min = Emin;
  R = sqrt(R2);
  tan2R = tan(R);
  tan2R *= tan2R;

  // browse protocones
  // for each of them, build the list of particles in them
  for (vector<CSphmomentum>::iterator p_it = protocones->begin();p_it != protocones->end();p_it++){
    // initialise variables
    c = &(*p_it);

    // browse particles to create cone contents
    // note that jet is always initialised with default values at this level
    jet_candidate.v = CSphmomentum();
    jet_candidate.contents.clear();
    for (i=0;i<n_left;i++){
      v = &(p_remain[i]);
      if (is_closer(v, c, tan2R)){
	jet_candidate.contents.push_back(v->parent_index);
	jet_candidate.v+= *v;
	v->index=0;
      }
    }
    jet_candidate.n=jet_candidate.contents.size();

    // compute Etilde for that jet.
    // we can't do that before as it requires knowledge of the jet axis
    // which has just been computed.
    compute_Etilde(jet_candidate);

    // set the momentum in protocones 
    // (it was only known through its spatial coordinates up to now)
    *c = jet_candidate.v;
    c->build_thetaphi();

    // set the jet range
    jet_candidate.range=CSphtheta_phi_range(c->_theta,c->_phi,R);

    // check that the protojet has large enough pt
    if (jet_candidate.v.E<E_min)
      continue;

    // assign the split-merge (or progressive-removal) squared scale variable
    if (_user_scale) {
      // sm_var2 is the signed square of the user scale returned
      // for the jet candidate
      jet_candidate.sm_var2 = (*_user_scale)(jet_candidate);
      jet_candidate.sm_var2 *= abs(jet_candidate.sm_var2);
    } else {
      jet_candidate.sm_var2 = get_sm_var2(jet_candidate.v, jet_candidate.E_tilde);
    }

    // now check if it is possibly the hardest
    if ((! found_jet) ||
	(_user_scale ? _user_scale->is_larger(jet_candidate, jet)
	             : ptcomparison(jet_candidate, jet))){
      jet = jet_candidate;
      found_jet = true;
    }
  }

  // make sure at least one of the jets has passed the selection
  if (!found_jet) return 1;  

  // add the jet to the list of jets
  jets.push_back(jet);
  jets[jets.size()-1].v.build_thetaphi();
  jets[jets.size()-1].v.build_norm();

#ifdef DEBUG_SPLIT_MERGE
  cout << "PR-Jet " << jets.size() << " [size " << jet.contents.size() << "]:";
#endif
    
  // update the list of what particles are left
  int p_remain_index = 0;
  int contents_index = 0;
  //sort(next_jet.contents.begin(),next_jet.contents.end());
  for (int index=0;index<n_left;index++){
    if ((contents_index<(int) jet.contents.size()) &&
	(p_remain[index].parent_index == jet.contents[contents_index])){
      // this particle belongs to the newly found jet
      // set pass in initial list
      particles[p_remain[index].parent_index].index = n_pass;
#ifdef DEBUG_SPLIT_MERGE
      cout << " " << jet.contents[contents_index];
#endif
      contents_index++;
    } else {
      // this particle still has to be clustered
      p_remain[p_remain_index] = p_remain[index];
      p_remain[p_remain_index].parent_index = p_remain[index].parent_index;
      p_remain[p_remain_index].index=1;
      p_remain_index++;
    }
  }
  p_remain.resize(n_left-jet.contents.size());
  n_left = p_remain.size();
  jets[jets.size()-1].pass = particles[jet.contents[0]].index;

  // update list of included particles
  n_pass++;

#ifdef DEBUG_SPLIT_MERGE
  cout << endl;
#endif

  // male sure the list of uncol_hard particles (used for the next
  // stable cone finding) is updated [could probably be more
  // efficient]
  merge_collinear_and_remove_soft();
  
  return 0;
}

/*
 * really do the splitting and merging
 * At the end, the vector jets is filled with the jets found.
 * the 'contents' field of each jets contains the indices
 * of the particles included in that jet. 
 *  - overlap_tshold    threshold for splitting/merging transition
 *  - Emin              minimal energy allowed for jets
 * return the number of jets is returned
 ******************************************************************/
int CSphsplit_merge::perform(double overlap_tshold, double Emin){
  // iterators for the 2 jets
  cjet_iterator j1;
  cjet_iterator j2;

  E_min = Emin;

  if (candidates->size()==0)
    return 0;

  if (overlap_tshold>=1.0 || overlap_tshold <= 0) {
    ostringstream message;
    message << "Illegal value for overlap_tshold, f = " << overlap_tshold;
    message << "  (legal values are 0<f<1)";
    throw siscone::Csiscone_error(message.str());
  }

  // overlap (the contents of this variable depends on the choice for
  // the split--merge variable.)
  // Note that the square of the ovelap is used
  double overlap2;

  // avoid to compute tshold*tshold at each overlap
  double overlap_tshold2 = overlap_tshold*overlap_tshold;

  do{
    if (candidates->size()>0){
      // browse for the first jet
      j1 = candidates->begin();
      
      // if hardest jet does not pass threshold then nothing else will
      // either so one stops the split merge.
      //if (j1->sm_var2<SM_var2_hardest_cut_off) {break;}
      if (j1->sm_var2<SM_var2_hardest_cut_off) {break;}

      // browse for the second jet
      j2 = j1;
      j2++;
      int j2_relindex = 1; // used only in ifdef, but costs little so keep it outside

      while (j2 != candidates->end()){
#ifdef DEBUG_SPLIT_MERGE
        if (j2_relindex==1) show();
        cout << "check overlap between cdt 1 and cdt " << j2_relindex+1 << " with overlap " << endl;
#endif
	// check overlapping
	if (get_overlap(*j1, *j2, &overlap2)){
	  // check if overlapping energy passes threshold
	  // Note that this depends on the ordering variable
#ifdef DEBUG_SPLIT_MERGE
          cout << "overlap between cdt 1 and cdt " << j2_relindex+1 << " with overlap " 
               << sqrt(overlap2)/j2->v.E << endl<<endl;
#endif
	  // We use the energy for the overlap computation
	  if (overlap2<overlap_tshold2*sqr(j2->v.E)){
#ifdef DEBUG_SPLIT_MERGE
            cout << "  --> split" << endl<<endl;
#endif
	    // split jets
	    split(j1, j2);
	    
	    // update iterators
	    j2 = j1 = candidates->begin();
            j2_relindex = 0;
	  } else {
#ifdef DEBUG_SPLIT_MERGE
            cout << "  --> merge" << endl<<endl;
#endif
	    // merge jets
	    merge(j1, j2);
	    
	    // update iterators
	    j2 = j1 = candidates->begin();
            j2_relindex = 0;
	  }
	}
        // watch out: split/merge might have caused new jets with E <
        // Emin to disappear, so the total number of jets may
        // have changed by more than expected and j2 might already by
        // the end of the candidates list...
        j2_relindex++;
	if (j2 != candidates->end()) j2++;
      } // end of loop on the second jet
      
      if (j1 != candidates->end()) {
        // all "second jet" passed without overlapping
        // (otherwise we won't leave the j2 loop)
        // => store jet 1 as real jet
        jets.push_back(*j1);
        jets[jets.size()-1].v.build_thetaphi();
        jets[jets.size()-1].v.build_norm();
        // a bug where the contents has non-zero size has been cropping
        // up in many contexts -- so catch it!
        assert(j1->contents.size() > 0);
        jets[jets.size()-1].pass = particles[j1->contents[0]].index;
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
        cand_refs.erase(j1->v.ref);
#endif
        candidates->erase(j1);
      }
    }
  } while (candidates->size()>0);

  // sort jets by Energy
  sort(jets.begin(), jets.end(), jets_E_less);
#ifdef DEBUG_SPLIT_MERGE
      show();
#endif

  return jets.size();
}



// save the event on disk
//  - flux   stream used to save jet contents
//--------------------------------------------
int CSphsplit_merge::save_contents(FILE *flux){
  jet_iterator it_j;
  CSphjet *j1;
  int i1, i2;

  fprintf(flux, "# %d jets found\n", (int) jets.size());
  fprintf(flux, "# columns are: px, py, pz, E and number of particles for each jet\n");
  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j1 = &(*it_j);
    fprintf(flux, "%e\t%e\t%e\t%e\t%d\n", 
	    j1->v.px, j1->v.py, j1->v.pz, j1->v.E, j1->n);
  }
  
  fprintf(flux, "# jet contents\n");
  fprintf(flux, "# columns are: px, py, pz, E, particle index and jet number\n");
  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j1 = &(*it_j);
    for (i2=0;i2<j1->n;i2++)
      fprintf(flux, "%e\t%e\t%e\t%e\t%d\t%d\n", 
      	      particles[j1->contents[i2]].px, particles[j1->contents[i2]].py,
      	      particles[j1->contents[i2]].pz, particles[j1->contents[i2]].E,
      	      j1->contents[i2], i1);
  }
  
  return 0;
}


// show current jets/candidate status
//------------------------------------
int CSphsplit_merge::show(){
  jet_iterator it_j;
  cjet_iterator it_c;
  CSphjet *j;
  const CSphjet *c;
  int i1, i2;

  for (it_j = jets.begin(), i1=0 ; it_j != jets.end() ; it_j++, i1++){
    j = &(*it_j);
    fprintf(stdout, "jet %2d: %e\t%e\t%e\t%e\t", i1+1,
	    j->v.px, j->v.py, j->v.pz, j->v.E);

    unsigned int phirange=j->range.phi_range;
    for (i2=0;i2<32;i2++) fprintf(stdout, "%d", (phirange&(1<<i2)) >> i2 );
    fprintf(stdout, "\t");
    unsigned int thetarange=j->range.theta_range;
    for (i2=0;i2<32;i2++) fprintf(stdout, "%d", (thetarange&(1<<i2)) >> i2);
    fprintf(stdout, "\t");
    
    for (i2=0;i2<j->n;i2++)
      fprintf(stdout, "%d ", j->contents[i2]);
    fprintf(stdout, "\n");
  }
  
  for (it_c = candidates->begin(), i1=0 ; it_c != candidates->end() ; it_c++, i1++){
    c = &(*it_c);
    fprintf(stdout, "cdt %2d: %e\t%e\t%e\t%e\t%e\t", i1+1,
	    c->v.px, c->v.py, c->v.pz, c->v.E, sqrt(c->sm_var2));

    unsigned int phirange=c->range.phi_range;
    for (i2=0;i2<32;i2++) fprintf(stdout, "%d", (phirange&(1<<i2)) >> i2 );
    fprintf(stdout, "\t");
    unsigned int thetarange=c->range.theta_range;
    for (i2=0;i2<32;i2++) fprintf(stdout, "%d", (thetarange&(1<<i2)) >> i2);
    fprintf(stdout, "\t");
    
    for (i2=0;i2<c->n;i2++)
      fprintf(stdout, "%d ", c->contents[i2]);
    fprintf(stdout, "\n");
  }
  
  fprintf(stdout, "\n");
  return 0;
}


// get the overlap between 2 jets
//  - j1        first jet
//  - j2        second jet
//  - overlap2  returned overlap^2 (determined by the choice of SM variable)
// return true if overlapping, false if disjoint
//---------------------------------------------------------------------
bool CSphsplit_merge::get_overlap(const CSphjet &j1, const CSphjet &j2, double *overlap2){
  // check if ranges overlap
  if (!is_range_overlap(j1.range,j2.range))
    return false;

  int i1,i2;
  bool is_overlap;

  // initialise
  i1=i2=idx_size=0;
  is_overlap = false;
  CSphmomentum v;

  // compute overlap
  // at the same time, we store union in indices
  do{
    if (j1.contents[i1]<j2.contents[i2]){
      indices[idx_size] = j1.contents[i1];
      i1++;
    } else if (j1.contents[i1]>j2.contents[i2]){
      indices[idx_size] = j2.contents[i2];
      i2++;
    } else { // (j1.contents[i1]==j2.contents[i2])
      v += particles[j2.contents[i2]];
      indices[idx_size] = j2.contents[i2];
      i1++;
      i2++;
      is_overlap = true;
    }
    idx_size++;
  } while ((i1<j1.n) && (i2<j2.n));

  // finish computing union
  // (only needed if overlap !)
  if (is_overlap){
    while (i1<j1.n){
      indices[idx_size] = j1.contents[i1];
      i1++;
      idx_size++;
    }
    while (i2<j2.n){
      indices[idx_size] = j2.contents[i2];
      i2++;
      idx_size++;
    }
  }

  // assign the overlapping var as return variable
  (*overlap2) = sqr(v.E); //get_sm_var2(v, E_tilde);

  return is_overlap;
}



// split the two given jet.
// during this procedure, the jets j1 & j2 are replaced
// by 2 new jets. Common particles are associted to the 
// closest initial jet.
//  - it_j1  iterator of the first jet in 'candidates'
//  - it_j2  iterator of the second jet in 'candidates'
//  - j1     first jet (CSphjet instance)
//  - j2     second jet (CSphjet instance)
// return true on success, false on error
////////////////////////////////////////////////////////////////
bool CSphsplit_merge::split(cjet_iterator &it_j1, cjet_iterator &it_j2){
  int i1, i2;
  CSphjet jet1, jet2;
  double E1_weight, E2_weight;
  CSphmomentum tmp;
  CSphmomentum *v;

  // shorthand to avoid having "->" all over the place
  const CSphjet & j1 = * it_j1;
  const CSphjet & j2 = * it_j2;

  i1=i2=0;
  jet2.v = jet1.v = CSphmomentum();

  // compute centroids
  // When use_E_weighted_splitting is activated, the
  // "geometrical" distance is weighted by the inverse
  // of the E of the protojet
  // This is stored in E{1,2}_weight
  E1_weight = (use_E_weighted_splitting) ? 1.0/j1.v.E/j1.v.E : 1.0;
  E2_weight = (use_E_weighted_splitting) ? 1.0/j2.v.E/j2.v.E : 1.0;

  // compute jet splitting
  do{
    if (j1.contents[i1]<j2.contents[i2]){
      // particle i1 belong only to jet 1
      v = &(particles[j1.contents[i1]]);
      jet1.contents.push_back(j1.contents[i1]);
      jet1.v += *v;
      //jet1.pt_tilde += pt[j1.contents[i1]];
      i1++;
      jet1.range.add_particle(v->_theta,v->_phi);
    } else if (j1.contents[i1]>j2.contents[i2]){
      // particle i2 belong only to jet 2
      v = &(particles[j2.contents[i2]]);
      jet2.contents.push_back(j2.contents[i2]);
      jet2.v += *v;
      //jet2.pt_tilde += pt[j2.contents[i2]];
      i2++;
      jet2.range.add_particle(v->_theta,v->_phi);
    } else { // (j1.contents[i1]==j2.contents[i2])
      // common particle, decide which is the closest centre
      v = &(particles[j1.contents[i1]]);

      //TODO: improve this brutal use of atan2 and sqrt !!!!

      //? what when == ? 
      // When use_E_weighted_splitting is activated, the
      // "geometrical" distance is weighted by the inverse
      // of the E of the protojet
      double d1 = get_distance(&(j1.v), v)*E1_weight;
      double d2 = get_distance(&(j2.v), v)*E2_weight;
      // do bookkeeping on most ambiguous split
      if (fabs(d1-d2) < most_ambiguous_split) 
        most_ambiguous_split = fabs(d1-d2);

      if (d1<d2){
	// particle i1 belong only to jet 1
	jet1.contents.push_back(j1.contents[i1]);
	jet1.v += *v;
	//jet1.pt_tilde += pt[j1.contents[i1]];
	jet1.range.add_particle(v->_theta,v->_phi);
      } else {
	// particle i2 belong only to jet 2
	jet2.contents.push_back(j2.contents[i2]);
	jet2.v += *v;
	//jet2.pt_tilde += pt[j2.contents[i2]];
	jet2.range.add_particle(v->_theta,v->_phi);
      }      

      i1++;
      i2++;
    }
  } while ((i1<j1.n) && (i2<j2.n));
  
  while (i1<j1.n){
    v = &(particles[j1.contents[i1]]);
    jet1.contents.push_back(j1.contents[i1]);
    jet1.v += *v;
    //jet1.pt_tilde += pt[j1.contents[i1]];
    i1++;
    jet1.range.add_particle(v->_theta,v->_phi);
  }
  while (i2<j2.n){
    v = &(particles[j2.contents[i2]]);
    jet2.contents.push_back(j2.contents[i2]);
    jet2.v += *v;
    //jet2.pt_tilde += pt[j2.contents[i2]];
    i2++;
    jet2.range.add_particle(v->_theta,v->_phi);
  }

  // finalise jets
  jet1.n = jet1.contents.size();
  jet2.n = jet2.contents.size();

  // now the jet axis is known, we can compute Etilde
  compute_Etilde(jet1);
  compute_Etilde(jet2);

  // remove previous jets
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  cand_refs.erase(j1.v.ref);
  cand_refs.erase(j2.v.ref);
#endif
  candidates->erase(it_j1);
  candidates->erase(it_j2);

  // reinsert new ones
  insert(jet1);
  insert(jet2);

  return true;
}

// merge the two given jet.
// during this procedure, the jets j1 & j2 are replaced
// by 1 single jets containing both of them.
//  - it_j1  iterator of the first jet in 'candidates'
//  - it_j2  iterator of the second jet in 'candidates'
// return true on success, false on error
////////////////////////////////////////////////////////////////
bool CSphsplit_merge::merge(cjet_iterator &it_j1, cjet_iterator &it_j2){
  CSphjet jet;
  int i;

  // build new jet
  // note: particles within j1 & j2 have already been stored in indices
  for (i=0;i<idx_size;i++){
    jet.contents.push_back(indices[i]);
    jet.v += particles[indices[i]];
    //jet.pt_tilde += pt[indices[i]];
  }
  jet.n = jet.contents.size();

  compute_Etilde(jet);

  // deal with ranges
  jet.range = range_union(it_j1->range, it_j2->range);

  // remove old candidates
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if (merge_identical_protocones){
    cand_refs.erase(it_j1->v.ref);
    cand_refs.erase(it_j2->v.ref);
  }
#endif
  candidates->erase(it_j1);
  candidates->erase(it_j2);

  // reinsert new candidate
  insert(jet);

  return true;
}

/**
 * Check whether or not a jet has to be inserted in the 
 * list of protojets. If it has, set its sm_variable and
 * insert it to the list of protojets.
 */
bool CSphsplit_merge::insert(CSphjet &jet){

  // eventually check that no other candidate are present with the
  // same cone contents. We recall that this automatic merging of
  // identical protocones can lead to infrared-unsafe situations.
#ifdef ALLOW_MERGE_IDENTICAL_PROTOCONES
  if ((merge_identical_protocones) && (!cand_refs.insert(jet.v.ref).second))
    return false;
#endif

  // check that the protojet has large enough energy
  if (jet.v.E<E_min)
    return false;

  // assign SM variable
  jet.sm_var2 = get_sm_var2(jet.v, jet.E_tilde);

  // insert the jet.
  candidates->insert(jet);

  return true;
}

/**
 * given a 4-momentum and its associated pT, return the 
 * variable that has to be used for SM
 * \param v          4 momentum of the protojet
 * \param pt_tilde   pt_tilde of the protojet
 */
double CSphsplit_merge::get_sm_var2(CSphmomentum &v, double &E_tilde){
  switch(ptcomparison.split_merge_scale) {
  case SM_E:       return v.E*v.E;
  case SM_Etilde:  return E_tilde*E_tilde;
  default:
    throw siscone::Csiscone_error("Unsupported split-merge scale choice: "
				  + ptcomparison.SM_scale_name());
  }

  //return 0.0;
}



/// compute Etilde for a given jet
void CSphsplit_merge::compute_Etilde(CSphjet &jet){
  jet.v.build_norm();
  jet.E_tilde=0.0;
  CSph3vector jet_axis = jet.v;
  //if (jet.v._norm==0){
  //  jet_axis = CSph3vector(0.0,0.0,0.0);
  //} else {
  jet_axis/=jet.v.E;
    //}
  //cout << "~~~ Axis: " << jet.v.px << " " << jet.v.py << " " << jet.v.pz << " " << jet.v._norm << endl;
  //cout << "~~~ Axis: " << jet_axis.px << " " << jet_axis.py << " " << jet_axis.pz << endl;
  for (vector<int>::iterator cont_it=jet.contents.begin(); cont_it!=jet.contents.end(); cont_it++){
    const CSphmomentum &p = particles[*cont_it];
    jet.E_tilde+=p.E*(1.0+norm2_cross_product3(p,jet_axis)/particles_norm2[*cont_it]);
  }
}

}
