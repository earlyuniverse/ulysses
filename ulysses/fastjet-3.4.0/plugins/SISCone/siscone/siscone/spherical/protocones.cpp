///////////////////////////////////////////////////////////////////////////////
// File: protocones.cpp                                                      //
// Description: source file for stable cones determination (Cstable_cones)   //
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

/*******************************************************
 * Introductory note:                                  *
 * Since this file has many member functions, we have  *
 * structured them in categories:		       *
 * INITIALISATION FUNCTIONS                            *
 *  - ctor()                                           *
 *  - ctor(particle_list)                              *
 *  - dtor()                                           *
 *  - init(particle_list)                              *
 * ALGORITHM MAIN ENTRY                                *
 *  - get_stable_cone(radius)                          *
 * ALGORITHM MAIN STEPS                                *
 *  - init_cone()                                      *
 *  - test_cone()                                      *
 *  - update_cone()                                    *
 *  - proceed_with_stability()                         *
 * ALGORITHM MAIN STEPS FOR COCIRCULAR SITUATIONS      *
 *  - cocircular_pt_less(v1, v2)                       *
 *  - prepare_cocircular_list()                        *
 *  - test_cone_cocircular()                           *
 *  - test_stability(candidate, border_list)           *
 *  - updat_cone_cocircular()                          *
 * RECOMPUTATION OF CONE CONTENTS                      *
 *  - compute_cone_contents()                          *
 *  - recompute_cone_contents()                        *
 *  - recompute_cone_contents_if_needed()              *
 * VARIOUS TOOLS                                       *
 *  - circle_intersect()                               *
 *  - is_inside()                                      *
 *  - abs_dangle()                                     *
 *******************************************************/

#include <siscone/defines.h>
#include <siscone/siscone_error.h>
#include <siscone/circulator.h>
#include "protocones.h"
#include <math.h>
#include <iostream>
#include <algorithm>

namespace siscone_spherical{

using namespace std;

/**********************************************************************
 * CSphstable_cones implementation                                    *
 * Computes the list of stable comes from a particle list.            *
 * This class does the first fundamental task of te cone algorithm:   *
 * it is used to compute the list of stable cones given a list        *
 * of particles.                                                      *
 **********************************************************************/

////////////////////////////////////////////////////////
// INITIALISATION FUNCTIONS                           //
//  - ctor()                                          //
//  - ctor(particle_list)                             //
//  - dtor()                                          //
//  - init(particle_list)                             //
////////////////////////////////////////////////////////

// default ctor
//--------------
CSphstable_cones::CSphstable_cones(){
  nb_tot = 0;
  hc = NULL;
}

// ctor with initialisation
//--------------------------
CSphstable_cones::CSphstable_cones(vector<CSphmomentum> &_particle_list)
  : CSphvicinity(_particle_list){

  nb_tot = 0;
  hc = NULL;
}

// default dtor
//--------------
CSphstable_cones::~CSphstable_cones(){
  if (hc!=NULL) delete hc;
}

/*
 * initialisation
 *  - _particle_list  list of particles
 *  - _n              number of particles
 *********************************************************************/
void CSphstable_cones::init(vector<CSphmomentum> &_particle_list){
  // check already allocated mem
  if (hc!=NULL){
    delete hc;
  }
  if (protocones.size()!=0)
    protocones.clear();

  multiple_centre_done.clear();

  // initialisation
  set_particle_list(_particle_list);
}


////////////////////////////////////////////////////////
// ALGORITHM MAIN ENTRY                               //
//  - get_stable_cone(radius)                         //
////////////////////////////////////////////////////////

/*
 * compute stable cones.
 * This function really does the job i.e. computes
 * the list of stable cones (in a seedless way)
 *  - _radius:  radius of the cones
 * The number of stable cones found is returned
 *********************************************************************/
int CSphstable_cones::get_stable_cones(double _radius){
  int p_idx;

  // check if everything is correctly initialised
  if (n_part==0){
    return 0;
  }

  R  = _radius;
  R2 = R*R;
  tan2R = tan(R);
  tan2R *= tan2R;

  // allow hash for cones candidates
  hc = new sph_hash_cones(n_part, R);

  // browse all particles
  for (p_idx=0;p_idx<n_part;p_idx++){
    // step 0: compute the child list CL.
    //         Note that this automatically sets the parent P
    build(&plist[p_idx], 2.0*R);

    // special case: 
    //   if the vicinity is empty, the parent particle is a 
    //   stable cone by itself. Add it to protocones list.
    if (vicinity_size==0){
      protocones.push_back(*parent);
      continue;
    }

#ifdef DEBUG_STABLE_CONES
    cout << endl << endl;
    cout << "plot 'particles.dat' u 2:1 pt 1 ps 3" << endl;
    cout << "set label 1 'x' at " << parent->_phi << ", " << parent->_theta << endl;
#endif

    // step 1: initialise with the first cone candidate
    init_cone();

    do{
      // step 2: test cone stability for that pair (P,C)
      test_cone();

      // step 3: go to the next cone child candidate C
    } while (!update_cone());
  }

  return proceed_with_stability();
}


////////////////////////////////////////////////////////
// ALGORITHM MAIN STEPS                               //
//  - init_cone()                                     //
//  - test_cone()                                     //
//  - update_cone()                                   //
//  - proceed_with_stability()                        //
////////////////////////////////////////////////////////

/*
 * initialise the cone.
 * We take the first particle in the angular ordering to compute 
 * this one
 * return 0 on success, 1 on error
 *********************************************************************/
int CSphstable_cones::init_cone(){
  // The previous version of the algorithm was starting the 
  // loop around vicinity elements with the "most isolated" child.
  // given the nodist method to calculate the cone contents, we no
  // longer need to worry about which cone comes first...
  first_cone=0;

  // now make sure we have lists of the cocircular particles
  prepare_cocircular_lists();

  //TODO? deal with a configuration with only degeneracies ? 
  // The only possibility seems a regular hexagon with a parent point
  // in the centre. And this situation is by itself unclear.
  // Hence, we do nothing here !

  // init set child C
  centre = vicinity[first_cone];
  child = centre->v;
  centre_idx = first_cone;

  // build the initial cone (nodist: avoids calculating distances --
  // just deduces contents by circulating around all in/out operations)
  // this function also sets the list of included particles
  compute_cone_contents();
  
  return 0;
}


/*
 * test cones.
 * We check if the cone(s) built with the present parent and child 
 * are stable
 * return 0 on success 1 on error
 *********************************************************************/
int CSphstable_cones::test_cone(){
  siscone::Creference weighted_cone_ref;
  
  // depending on the side we are taking the child particle,
  // we test different configuration.
  // Each time, two configurations are tested in such a way that
  // all 4 possible cases (parent or child in or out the cone)
  // are tested when taking the pair of particle parent+child
  // and child+parent.

  // here are the tests entering the first series:
  //  1. check if the cone is already inserted
  //  2. check cone stability for the parent and child particles
  
  //UPDATED(see below): if (centre->side){
  //UPDATED(see below):   // test when both particles are not in the cone
  //UPDATED(see below):   // or when both are in.
  //UPDATED(see below):   // Note: for the totally exclusive case, test emptyness before
  //UPDATED(see below):   cone_candidate = cone;
  //UPDATED(see below):   if (cone.ref.not_empty()){
  //UPDATED(see below):     hc->insert(&cone_candidate, parent, child, false, false);
  //UPDATED(see below):   }
  //UPDATED(see below): 
  //UPDATED(see below):   cone_candidate = cone;
  //UPDATED(see below):   cone_candidate+= *parent + *child;
  //UPDATED(see below):   hc->insert(&cone_candidate, parent, child, true, true);
  //UPDATED(see below): } else {
  //UPDATED(see below):   // test when 1! of the particles is in the cone
  //UPDATED(see below):   cone_candidate = cone + *parent;
  //UPDATED(see below):   hc->insert(&cone_candidate, parent, child, true, false);
  //UPDATED(see below): 
  //UPDATED(see below):   cone_candidate = cone + *child;
  //UPDATED(see below):   hc->insert(&cone_candidate, parent, child, false, true);
  //UPDATED(see below): }
  //UPDATED(see below): 
  //UPDATED(see below): nb_tot+=2;

  // instead of testing 2 inclusion/exclusion states for every pair, we test the 4 of them
  // when the parent has an energy bigger than the child
  if (parent->E >= child->E){
    // test when both particles are not in the cone
    // Note: for the totally exclusive case, test emptiness before
    cone_candidate = cone;
    if (cone.ref.not_empty()){
      hc->insert(&cone_candidate, parent, child, false, false);
    }

    // test when 1! of the particles is in the cone
    cone_candidate += *parent;
    hc->insert(&cone_candidate, parent, child, true, false);

    cone_candidate  = cone;
    cone_candidate += *child;
    hc->insert(&cone_candidate, parent, child, false, true);

    // test when both are in.
    cone_candidate += *parent;
    hc->insert(&cone_candidate, parent, child, true, true);
    
    nb_tot += 4;
  }

  
  return 0;
}


/*
 * update the cone
 * go to the next child for that parent and update 'cone' appropriately
 * return 0 if update candidate found, 1 otherwise
 ***********************************************************************/
int CSphstable_cones::update_cone(){
#ifdef DEBUG_STABLE_CONES
  cout << "call 'circles_plot.gp' '" << centre->centre.px << "' '" 
       << centre->centre.py << "' '" << centre->centre.pz << "'" << endl
       << "pause -1 '(" << centre->angle << " " << (centre->side ? '+' : '-') << ")";
#endif

  // get the next child and centre
  centre_idx++;
  if (centre_idx==vicinity_size)
    centre_idx=0;
  if (centre_idx==first_cone)
    return 1;

  // update the cone w.r.t. the old child
  // only required if the old child is entering inside in which
  // case we need to add it. We also know that the child is 
  // inside iff its side is -.
  if (!centre->side){
#ifdef DEBUG_STABLE_CONES
    cout << " old_enter";
#endif
    // update cone
    cone += (*child);

    // update info on particles inside
    centre->is_inside->cone = true;

    // update stability check quantities
    dpt += fabs(child->px)+fabs(child->py)+fabs(child->pz);
  }

  // update centre and child to correspond to the new position
  centre = vicinity[centre_idx];
  child = centre->v;

  // check cocircularity
  // note that if cocirculaity is detected (i.e. if we receive 1
  // in the next test), we need to recall 'update_cone' directly
  // since tests and remaining part of te update has been performed
  //if (cocircular_check())
  if (cocircular_check()){
#ifdef DEBUG_STABLE_CONES
    cout << " Co-circular case detected" << endl;
#endif
    return update_cone();
  }

  // update the cone w.r.t. the new child
  // only required if the new child was already inside in which
  // case we need to remove it. We also know that the child is 
  // inside iff its side is +.
  if ((centre->side) && (cone.ref.not_empty())){
#ifdef DEBUG_STABLE_CONES
    cout << " new exit";
#endif

    // update cone
    cone -= (*child);

    // update info on particles inside
    centre->is_inside->cone = false;

    // update stability check quantities
    dpt += fabs(child->px)+fabs(child->py)+fabs(child->pz); //child->perp2();
  }

  // check that the addition and subtraction of vectors does
  // not lead to too much rounding error
  // for that, we compute the sum of pt modifications and of |pt|
  // since last recomputation and once the ratio overpasses a threshold
  // we recompute vicinity.
  if ((dpt>PT_TSHOLD*(fabs(cone.px)+fabs(cone.py)+fabs(cone.pz))) && (cone.ref.not_empty())){
    recompute_cone_contents();
  }
  if (cone.ref.is_empty()){
    cone = CSphmomentum();
    dpt=0.0;
  }

#ifdef DEBUG_STABLE_CONES
  cout << "'" << endl; 
#endif

  return 0; 
}


/*
 * compute stability of all enumerated candidates.
 * For all candidate cones which are stable w.r.t. their border particles,
 * pass the last test: stability with quadtree intersection
 ************************************************************************/
int CSphstable_cones::proceed_with_stability(){
  int i;
  sph_hash_element *elm;

  for (i=0;i<=hc->mask;i++){
    // test ith cell of the hash array
    elm = hc->hash_array[i];

    // browse elements therein
    while (elm!=NULL){
      // test stability
      if (elm->is_stable){
	// stability is not ensured by all pairs of "edges" already browsed
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
	//  => testing stability with quadtree intersection
	if (quadtree->circle_intersect(elm->eta, elm->phi, R2)==elm->ref)
#else
	//  => testing stability with the particle-list intersection
	if (circle_intersect(elm->centre)==elm->centre.ref)
#endif
	  protocones.push_back(CSphmomentum(elm->centre,1.0));
      }
      
      // jump to the next one
      elm = elm->next;
    }
  }
  
  // free hash
  // we do that at this level because hash eats rather a lot of memory
  // we want to free it before running the split/merge algorithm
#ifdef DEBUG_STABLE_CONES
  nb_hash_cones = hc->n_cones;
  nb_hash_occupied = hc->n_occupied_cells;
#endif

  delete hc;
  hc=NULL;

  return protocones.size();
}


////////////////////////////////////////////////////////
// ALGORITHM MAIN STEPS FOR COCIRCULAR SITUATIONS     //
//  - cocircular_pt_less(v1, v2)                      //
//  - prepare_cocircular_list()                       //
//  - test_cone_cocircular()                          //
//  - test_stability(candidate, border_vect)          //
//  - updat_cone_cocircular()                         //
////////////////////////////////////////////////////////

/// pt-ordering of momenta used for the cocircular case
//NEVER USED
//bool cocircular_pt_less(CSphmomentum *v1, CSphmomentum *v2){
//  return v1->perp2() < v2->perp2();
//}

/*
 * run through the vicinity of the current parent and for each child
 * establish which other members are cocircular... Note that the list
 * associated with each child contains references to vicinity
 * elements: thus two vicinity elements each associated with one given
 * particle may appear in a list -- this needs to be watched out for
 * later on... 
 **********************************************************************/
void CSphstable_cones::prepare_cocircular_lists() {
  siscone::circulator<vector<CSphvicinity_elm*>::iterator > here(vicinity.begin(), 
								 vicinity.begin(), 
								 vicinity.end());

  siscone::circulator<vector<CSphvicinity_elm*>::iterator > search(here);

  do {
    CSphvicinity_elm* here_pntr = *here();
    search.set_position(here);

    // search forwards for things that should have "here" included in
    // their cocircularity list
    while (true) {
      ++search;
      if ( siscone::abs_dphi((*search())->angle, here_pntr->angle) < 
	   here_pntr->cocircular_range 
           && search() != here()) {
        (*search())->cocircular.push_back(here_pntr);
      } else {
        break;
      }
    }

    // search backwards
    search.set_position(here);
    while (true) {
      --search;
      if ( siscone::abs_dphi((*search())->angle, here_pntr->angle) < 
	   here_pntr->cocircular_range 
           && search() != here()) {
        (*search())->cocircular.push_back(here_pntr);
      } else {
        break;
      }
    }

    ++here;
  } while (here() != vicinity.begin());
}

/* 
 * Testing cocircular configurations in p^3 time,
 * rather than 2^p time; we will test all contiguous subsets of points
 * on the border --- note that this is till probably overkill, since
 * in principle we only have to test situations where up to a
 * half-circle is filled (but going to a full circle is simpler)
 ******************************************************************/
void CSphstable_cones::test_cone_cocircular(CSphmomentum & borderless_cone,
					    list<CSphmomentum *> & border_list) {
  // in spherical coordinates, we don't have a universal x-y axis system
  // to measure the angles. So we first determine one minimising
  // the uncertainties
  CSph3vector angl_dir1, angl_dir2;
  centre->centre.get_angular_directions(angl_dir1, angl_dir2);
  angl_dir1/=angl_dir1._norm;
  angl_dir2/=angl_dir2._norm;

  // now we have te reference axis, create the CSphborder_store structure
  vector<CSphborder_store> border_vect;
  border_vect.reserve(border_list.size());
  for (list<CSphmomentum *>::iterator it = border_list.begin();
       it != border_list.end(); it++) {
    border_vect.push_back(CSphborder_store(*it, centre->centre, angl_dir1, angl_dir2));
  }

  // get them into order of angle
  sort(border_vect.begin(), border_vect.end());

  // set up some circulators, since these will help us go around the
  // circle easily
  siscone::circulator<vector<CSphborder_store>::iterator > 
    start(border_vect.begin(), border_vect.begin(),border_vect.end());
  siscone::circulator<vector<CSphborder_store>::iterator > mid(start), end(start);
  
  // test the borderless cone
  CSphmomentum candidate = borderless_cone;
  //candidate.build_etaphi();
  if (candidate.ref.not_empty())
    test_stability(candidate, border_vect);

  do {
    // reset status wrt inclusion in the cone
    mid = start;
    do {
      mid()->is_in = false;
    } while (++mid != start);

    // now run over all inclusion possibilities with this starting point
    candidate = borderless_cone;
    while (++mid != start) { 
      // will begin with start+1 and go up to start-1
      mid()->is_in = true;
      candidate += *(mid()->mom);
      test_stability(candidate, border_vect);
    }

  } while (++start != end);

  // mid corresponds to momentum that we need to include to get the
  // full cone
  mid()->is_in = true;
  candidate += *(mid()->mom);
  test_stability(candidate, border_vect);
}


/**
 * carry out the computations needed for the stability check of the
 * candidate, using the border_vect to indicate which particles
 * should / should not be in the stable cone; if the cone is stable
 * insert it into the hash.
 **********************************************************************/
void CSphstable_cones::test_stability(CSphmomentum & candidate, const vector<CSphborder_store> & border_vect) {
  
  // this almost certainly has not been done...
  //candidate.build_etaphi();

  bool stable = true;
  for (unsigned i = 0; i < border_vect.size(); i++) {
    if (is_closer(&candidate, border_vect[i].mom,tan2R) ^ (border_vect[i].is_in)) {
      stable = false;
      break; // it's unstable so there's no point continuing
    }
  }

  if (stable) hc->insert(&candidate);
}

/*
 * check if we are in a situation of cocircularity.
 * if it is the case, update and test in the corresponding way
 * return 'false' if no cocircularity detected, 'true' otherwise
 * Note that if cocircularity is detected, we need to 
 * recall 'update' from 'update' !!!
 ***************************************************************/
bool CSphstable_cones::cocircular_check(){
  // check if many configurations have the same centre.
  // if this is the case, branch on the algorithm for this
  // special case.
  // Note that those situation, being considered separately in 
  // test_cone_multiple, must only be considered here if all
  // angles are on the same side (this avoid multiple counting)

  if (centre->cocircular.empty()) return false;

  // first get cone into status required at end...
  if ((centre->side) && (cone.ref.not_empty())){
    // update cone
    cone -= (*child);

    // update info on particles inside
    centre->is_inside->cone = false;

    // update stability check quantities
    dpt += fabs(child->px)+fabs(child->py)+fabs(child->pz); //child->perp2();
  }


  // now establish the list of unique children in the list
  // first make sure parent and child are in!

  list<siscone::Cvicinity_inclusion *> removed_from_cone;
  list<siscone::Cvicinity_inclusion *> put_in_border;
  list<CSphmomentum *> border_list;
  
  CSphmomentum cone_removal;
  CSphmomentum border = *parent;
  border_list.push_back(parent);

  // make sure child appears in the border region
  centre->cocircular.push_back(centre);

  // now establish the full contents of the cone minus the cocircular
  // region and of the cocircular region itself
  for(list<CSphvicinity_elm *>::iterator it = centre->cocircular.begin();
      it != centre->cocircular.end(); it++) {

    if ((*it)->is_inside->cone) {
      cone_removal           += *((*it)->v);
      (*it)->is_inside->cone  = false;
      removed_from_cone.push_back((*it)->is_inside);
    }

    // if a point appears twice (i.e. with + and - sign) in the list of 
    // points on the border, we take care not to include it twice.
    // Note that this situation may appear when a point is at a distance
    // close to 2R from the parent
    if (!(*it)->is_inside->cocirc) {
      border += *((*it)->v);
      (*it)->is_inside->cocirc  = true;
      put_in_border.push_back((*it)->is_inside);
      border_list.push_back((*it)->v);
      //cout << "  adding particle " << (*it)->v->_theta << ", " << (*it)->v->_phi << " to the border list" << endl; 
    }
  }


  // figure out whether this pairing has been observed before
  CSphmomentum borderless_cone = cone;
  borderless_cone -= cone_removal;
  bool consider = true;
  for (unsigned int i=0;i<multiple_centre_done.size();i++){
    if ((multiple_centre_done[i].first ==borderless_cone.ref) &&
	(multiple_centre_done[i].second==border.ref))
      consider = false;
  }

  // now prepare the hard work
  if (consider) {
    // record the fact that we've now seen this combination
    multiple_centre_done.push_back(pair<siscone::Creference,siscone::Creference>(borderless_cone.ref, 
										 border.ref));

    // first figure out whether our cone momentum is good
    double local_dpt = fabs(cone_removal.px) + fabs(cone_removal.py);
    double total_dpt = dpt + local_dpt;

    recompute_cone_contents_if_needed(borderless_cone, total_dpt);
    if (total_dpt == 0) {
      // a recomputation has taken place -- so take advantage of this
      // and update the member cone momentum
      cone = borderless_cone + cone_removal;
      dpt  = local_dpt;
    }

    test_cone_cocircular(borderless_cone, border_list);
  }


  // relabel things that were in the cone but got removed
  for(list<siscone::Cvicinity_inclusion *>::iterator is_in = removed_from_cone.begin();
      is_in != removed_from_cone.end(); is_in++) {
    (*is_in)->cone = true;
  }

  // relabel things that got put into the border 
  for(list<siscone::Cvicinity_inclusion *>::iterator is_in = put_in_border.begin();
      is_in != put_in_border.end(); is_in++) {
    (*is_in)->cocirc = false;
  }

  // we're done with everything -- return true to signal to user that we've
  // been through the co-circularity rigmarole
  return true;
}


////////////////////////////////////////////////////////
// RECOMPUTATION OF CONE CONTENTS                     //
//  - compute_cone_contents()                         //
//  - recompute_cone_contents()                       //
//  - recompute_cone_contents_if_needed()             //
////////////////////////////////////////////////////////

/**
 * compute the cone contents by going once around the full set of
 * circles and tracking the entry/exit status each time 
 * given parent, child and centre compute the momentum
 * of the particle inside the cone
 * This sets up the inclusion information, which can then be directly 
 * used to calculate the cone momentum.
 **********************************************************************/
void CSphstable_cones::compute_cone_contents() {
  siscone::circulator<vector<CSphvicinity_elm*>::iterator > 
    start(vicinity.begin()+first_cone, vicinity.begin(), vicinity.end());

  siscone::circulator<vector<CSphvicinity_elm*>::iterator > here(start);

  // note that in the following algorithm, the cone contents never includes
  // the child. Indeed, if it has positive sign, then it will be set as
  // outside at the last step in the loop. If it has negative sign, then the 
  // loop will at some point go to the corresponding situation with positive
  // sign and set the inclusion status to 0.

  do {
    // as we leave this position a particle enters if its side is
    // negative (i.e. the centre is the one at -ve angle wrt to the
    // parent-child line
    if (!(*here())->side) ((*here())->is_inside->cone) = 1;
    
    // move on to the next position
    ++here;
    
    // as we arrive at this position a particle leaves if its side is positive
    if ((*here())->side) ((*here())->is_inside->cone) = 0;
  } while (here != start);

  // once we've reached the start the 'is_inside' information should be
  // 100% complete, so we can use it to calculate the cone contents
  // and then exit
  recompute_cone_contents();
  return;

}


/*
 * compute the cone momentum from particle list.
 * in this version, we use the 'pincluded' information
 * from the CSphvicinity class
 */
void CSphstable_cones::recompute_cone_contents(){
  unsigned int i;

  // set momentum to 0
  cone = CSphmomentum();

  // Important note: we can browse only the particles
  // in vicinity since all particles in the cone are
  // withing a distance 2R w.r.t. parent hence in vicinity.
  // Among those, we only add the particles for which 'is_inside' is true !
  // This methos rather than a direct comparison avoids rounding errors
  for (i=0;i<vicinity_size;i++){
    // to avoid double-counting, only use particles with + angle
    if ((vicinity[i]->side) && (vicinity[i]->is_inside->cone))
      cone += *vicinity[i]->v;
  }
  
  // set check variables back to 0
  dpt = 0.0;
}


/*
 * if we have gone beyond the acceptable threshold of change, compute
 * the cone momentum from particle list.  in this version, we use the
 * 'pincluded' information from the CSphvicinity class, but we don't
 * change the member cone, only the locally supplied one
 */
void CSphstable_cones::recompute_cone_contents_if_needed(CSphmomentum & this_cone, 
							 double & this_dpt){
  
  if (this_dpt > PT_TSHOLD*(fabs(this_cone.px)+fabs(this_cone.py))) {
    if (cone.ref.is_empty()) {
      this_cone = CSphmomentum();
    } else {
      // set momentum to 0
      this_cone = CSphmomentum();
      
      // Important note: we can browse only the particles
      // in vicinity since all particles in the this_cone are
      // withing a distance 2R w.r.t. parent hence in vicinity.
      // Among those, we only add the particles for which 'is_inside' is true !
      // This methos rather than a direct comparison avoids rounding errors
      for (unsigned int i=0;i<vicinity_size;i++){
        // to avoid double-counting, only use particles with + angle
        if ((vicinity[i]->side) && (vicinity[i]->is_inside->cone))
          this_cone += *vicinity[i]->v;
      }
      
    }
    // set check variables back to 0
    this_dpt = 0.0;
  }

}


////////////////////////////////////////////////////////
// VARIOUS TOOLS                                      //
//  - circle_intersect()                              //
//  - is_inside()                                     //
//  - abs_dangle()                                    //
////////////////////////////////////////////////////////


/*
 * circle intersection.
 * computes the intersection with a circle of given centre and radius.
 * The output takes the form of a checkxor of the intersection's particles
 *  - cx    circle centre x coordinate
 *  - cy    circle centre y coordinate
 * return the checkxor for the intersection
 ******************************************************************/
siscone::Creference CSphstable_cones::circle_intersect(CSph3vector &cone_centre){
  siscone::Creference intersection;
  int i;

  for (i=0;i<n_part;i++){
    // really check if the distance is less than R
    if (is_closer(&cone_centre, &(plist[i]), tan2R))
      intersection+=plist[i].ref;
  }
  
  return intersection;
}

}
