///////////////////////////////////////////////////////////////////////////////
// File: vicinity.cpp                                                        //
// Description: source file for particle vicinity (Cvicinity class)          //
// This file is part of the SISCone project.                                 //
// WARNING: this is not the main SISCone trunk but                           //
//          an adaptation to spherical coordinates                           //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006-2008 Gavin Salam and Gregory Soyez                          //
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

#include "vicinity.h"
#include <math.h>
#include <algorithm>
#include <iostream>

namespace siscone_spherical{

using namespace std;

/*************************************************************
 * CSphvicinity_elm implementation                           *
 * element in the vicinity of a parent.                      *
 * class used to manage one points in the vicinity           *
 * of a parent point.                                        *
 *************************************************************/

// ordering pointers to CSphvicinity_elm
//---------------------------------------
bool ve_less(CSphvicinity_elm *ve1, CSphvicinity_elm *ve2){
  return ve1->angle < ve2->angle;
}


/*************************************************************
 * CSphvicinity implementation                               *
 * list of element in the vicinity of a parent.              *
 * class used to manage the points which are in the vicinity *
 * of a parent point. The construction of the list can be    *
 * made from a list of points or from a quadtree.            *
 *************************************************************/

// default constructor
//---------------------
CSphvicinity::CSphvicinity(){
  n_part = 0;

  ve_list = NULL;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  quadtree = NULL;
#endif

  parent = NULL;
  VR2 = VR = 0.0;

}

// constructor with initialisation
//---------------------------------
CSphvicinity::CSphvicinity(vector<CSphmomentum> &_particle_list){
  parent = NULL;
  ve_list = NULL;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  quadtree = NULL;
#endif
  cosVR = VR2 = tan2R = VR = 0.0;

  set_particle_list(_particle_list);
}

// default destructor
//--------------------
CSphvicinity::~CSphvicinity(){
  if (ve_list!=NULL)
    delete[] ve_list;

#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  if (quadtree!=NULL)
    delete quadtree;
#endif
}

/*
 * set the particle_list
 *  - particle_list   list of particles (type CSphmomentum)
 *  - n               number of particles in the list
 ************************************************************/ 
void CSphvicinity::set_particle_list(vector<CSphmomentum> &_particle_list){
  int i,j;
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  double eta_max=0.0;
#endif
  
  // if the particle list is not empty, destroy it !
  if (ve_list!=NULL){
    delete[] ve_list;
  }
  vicinity.clear();
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  if (quadtree!=NULL)
    delete quadtree;
#endif

  // allocate memory array for particles
  // Note: - we compute max for |eta|
  //       - we allocate indices to particles
  n_part = 0;
  plist.clear();
  pincluded.clear();
  for (i=0;i<(int) _particle_list.size();i++){
    // if a particle is colinear with the beam (infinite rapidity)
    // we do not take it into account
    //if (fabs(_particle_list[i].pz)!=_particle_list[i].E){
      plist.push_back(_particle_list[i]);
      pincluded.push_back(siscone::Cvicinity_inclusion()); // zero inclusion status

      // the parent_index is handled in the split_merge because 
      // of our multiple-pass procedure.
      // Hence, it is not required here any longer.
      // plist[n_part].parent_index = i;
      plist[n_part].index = n_part;

      // make sure the reference is randomly created
      plist[n_part].ref.randomize();

#ifdef USE_QUADTREE_FOR_STABILITY_TEST
      if (fabs(plist[n_part].eta)>eta_max) eta_max=fabs(plist[n_part].eta);
#endif
      n_part++;
      //}
  }

  // allocate quadtree and vicinity_elm list
  // note: we set phi in [-pi:pi] as it is the natural range for atan2!
  ve_list = new CSphvicinity_elm[2*n_part];
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
  eta_max+=0.1;
  quadtree = new siscone::Cquadtree(0.0, 0.0, eta_max, M_PI);
#endif

  // append particle to the vicinity_elm list
  j = 0;
  for (i=0;i<n_part;i++){
#ifdef USE_QUADTREE_FOR_STABILITY_TEST
    quadtree->add(&plist[i]);
#endif
    ve_list[j].v = ve_list[j+1].v = &plist[i];
    ve_list[j].is_inside = ve_list[j+1].is_inside = &(pincluded[i]);
    j+=2;
  }

}


/*
 * build the vicinity list from a list of points.
 *  - _parent   reference particle
 *  - _VR       vicinity radius
 ************************************************************/
void CSphvicinity::build(CSphmomentum *_parent, double _VR){
  int i;

  // set parent and radius
  parent = _parent;

  VR  = _VR;
  VR2 = VR*VR;
  cosVR = cos(VR);
  R2  = 0.25*VR2;
  R   = 0.5*VR;
  double tmp = tan(R);
  tan2R = tmp*tmp;

  D2_R = 2.0*(1-cos(R));
  //tmp = sqrt(D2_R);
  inv_R_EPS_COCIRC  = 1.0 / R / EPSILON_COCIRCULAR;
  inv_R_2EPS_COCIRC = 0.5 / R / EPSILON_COCIRCULAR;

  // clear vicinity
  vicinity.clear();

  // init parent variables
  // we cpte the direction of the centre and two orthogonal ones 
  // to measure the angles. Those are taken orthogonal to the
  // axis of smallest components (of the centre) to increase precision
  parent_centre = (*parent)/parent->_norm;
  parent_centre.get_angular_directions(angular_dir1, angular_dir2);
  angular_dir1 /= angular_dir1._norm;
  angular_dir2 /= angular_dir2._norm;
  
  // really browse the particle list
  for (i=0;i<n_part;i++){
    append_to_vicinity(&plist[i]);
  }

  // sort the vicinity
  sort(vicinity.begin(), vicinity.end(), ve_less);

  vicinity_size = vicinity.size();
}


/// strictly increasing function of the angle 
//TODO//
inline double sort_angle(double s, double c){
  if (s==0) return (c>0) ? 0.0 : 2.0;
  double t=c/s;
  return (s>0) ? 1-t/(1+fabs(t)) : 3-t/(1+fabs(t));
}


/*
 * append a particle to the 'vicinity' list after
 * having computed the angular-ordering quantities
 *  - v   vector to test
 **********************************************************/
void CSphvicinity::append_to_vicinity(CSphmomentum *v){
  // skip the particle itself)
  if (v==parent)
    return;

  int i=2*(v->index);

  // compute the distance of the i-th particle with the parent
  double dot = dot_product3(parent_centre,*v);
  CSph3vector vnormal = *v;
  vnormal/=v->_norm;
  dot/=v->_norm;

  // really check if the distance is less than VR
  if (dot>cosVR){
    CSph3vector cross = cross_product3(parent_centre,vnormal);

    // for the centres
    CSph3vector median = (parent_centre+vnormal);
    double amplT = sqrt((tan2R*(1+dot)+(dot-1))*(1+dot));
    CSph3vector transverse = amplT*cross/cross._norm;

    // first angle (+)
    ve_list[i].centre = median + transverse;
    ve_list[i].centre.build_norm();
    ve_list[i].centre/=ve_list[i].centre._norm;    
    CSph3vector diff = ve_list[i].centre - parent_centre;
    //ve_list[i].angle = atan2(dot_product3(angular_dir2, diff),dot_product3(angular_dir1, diff));
    ve_list[i].angle = sort_angle(dot_product3(angular_dir2, diff),dot_product3(angular_dir1, diff));
    ve_list[i].side = true;
    ve_list[i].cocircular.clear();
    vicinity.push_back(&(ve_list[i]));

    // second angle (-)    
    ve_list[i+1].centre = median - transverse;
    ve_list[i+1].centre.build_norm();
    ve_list[i+1].centre/=ve_list[i+1].centre._norm;    
    diff = ve_list[i+1].centre - parent_centre;
    ve_list[i+1].angle = sort_angle(dot_product3(angular_dir2, diff),dot_product3(angular_dir1, diff));
    ve_list[i+1].side = false;
    ve_list[i+1].cocircular.clear();
    vicinity.push_back(&(ve_list[i+1]));

    // now work out the cocircularity range for the two points (range
    // of angle within which the points stay within a distance
    // EPSILON_COCIRCULAR of circule
    // P = parent; C = child; O = Origin (center of circle)
    CSph3vector OP = parent_centre - ve_list[i+1].centre;
    CSph3vector OC = vnormal - ve_list[i+1].centre;
    
    // two sources of error are (GPS CCN29-19) epsilon/(R sin theta)
    // and sqrt(2*epsilon/(R (1-cos theta))) and the way things work
    // out, it is the _smaller_ of the two that is relevant [NB have
    // changed definition of theta here relative to that used in
    // CCN29] [NB2: write things so as to avoid zero denominators and
    // to minimize the multiplications, divisions and above all sqrts
    // -- that means that c & s are defined including a factor of VR2]
    double inv_err1 = cross_product3(OP,OC)._norm * inv_R_EPS_COCIRC;
    double inv_err2_sq = (D2_R-dot_product3(OP,OC)) * inv_R_2EPS_COCIRC;
    ve_list[i].cocircular_range = siscone::pow2(inv_err1) > inv_err2_sq ? 
      1.0/inv_err1 : 
      sqrt(1.0/inv_err2_sq);
    ve_list[i+1].cocircular_range = ve_list[i].cocircular_range;
  }
}

}
