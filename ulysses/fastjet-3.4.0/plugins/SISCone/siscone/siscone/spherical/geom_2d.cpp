///////////////////////////////////////////////////////////////////////////////
// File: geom_2d.cpp                                                         //
// Description: source file for two-dimensional geometry tools               //
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

#include "geom_2d.h"
#include <algorithm>

namespace siscone_spherical{

#define PHI_RANGE_MASK 0xFFFFFFFF

/*********************************************************
 * class CSphtheta_phi_range implementation              *
 * class for holding a covering range in eta-phi         *
 *                                                       *
 * This class deals with ranges in the eta-phi plane. It *
 * implements methods to test if two ranges overlap and  *
 * to take the union of two overlapping intervals.       *
 *********************************************************/

using namespace std;

// static member default init
//----------------------------
double CSphtheta_phi_range::theta_min = 0.0;
double CSphtheta_phi_range::theta_max = M_PI;

// default ctor
//--------------
CSphtheta_phi_range::CSphtheta_phi_range(){
  theta_range = 0;
  phi_range = 0;
}

// ctor with initialisation
// we initialise with a centre (in eta,phi) and a radius
//  - c_theta  theta coordinate of the centre
//  - c_phi    phi coordinate of the centre
//  - R        radius
//-------------------------------------------------------
CSphtheta_phi_range::CSphtheta_phi_range(double c_theta, double c_phi, double R){
  // determination of the eta range
  //-------------------------------
  double xmin = max(c_theta-R,theta_min+0.00001);
  double xmax = min(c_theta+R,theta_max-0.00001);

  unsigned int cell_min = get_theta_cell(xmin);
  unsigned int cell_max = get_theta_cell(xmax);

  // warning: if cell_max==2^31, 2*cell_max==0 hence, 
  // even if the next formula is formally (2*cell_max-cell_min),
  // expressing it as (cell_max-cell_min)+cell_max is safe.
  theta_range = (cell_max-cell_min)+cell_max;

  // determination of the phi range
  // !! taking care of periodicity !!
  // !! and the theta dependence   !!
  //---------------------------------
  double ymin,ymax;
  double extra = asin(R/M_PI);
  // if the theta range comes too close to the endpoints (theta=0 or
  // theta=pi), then keep the full phi range
  if (xmin<=theta_min+extra){
    ymin = -M_PI+0.00001;
    ymax =  M_PI-0.00001;
  } else if (xmax>=theta_max-extra){
    ymin = -M_PI+0.00001;
    ymax =  M_PI-0.00001;
  } else {
    extra = max(1.0/sin(xmin), 1.0/sin(xmax));
    ymin = (c_phi-R)*extra;
    while (ymin<-M_PI) ymin+=twopi;
    while (ymin> M_PI) ymin-=twopi;
    ymax = (c_phi-R)*extra;
    while (ymax<-M_PI) ymax+=twopi;
    while (ymax> M_PI) ymax-=twopi;
  }
  cell_min = get_phi_cell(ymin);
  cell_max = get_phi_cell(ymax);

  // Also, if the interval goes through pi, inversion is needed
  if (ymax>ymin)
    phi_range = (cell_max-cell_min)+cell_max;
  else {
    phi_range = (cell_min==cell_max) 
      ? PHI_RANGE_MASK
      : ((PHI_RANGE_MASK^(cell_min-cell_max)) + cell_max);
  }
}

// assignment of range
//  - r   range to assign to current one
//---------------------------------------
CSphtheta_phi_range& CSphtheta_phi_range::operator = (const CSphtheta_phi_range &r){
  theta_range = r.theta_range;
  phi_range = r.phi_range;

  return *this;
}

// add a particle to the range
//  - eta  eta coordinate of the particle
//  - phi  phi coordinate of the particle
// \return 0 on success, 1 on error
//----------------------------------------
int CSphtheta_phi_range::add_particle(const double theta, const double phi){
  // get the theta cell
  unsigned int theta_cell = get_theta_cell(theta);
  
  // deal with the eta coordinate
  theta_range |= theta_cell;

  // deal with the phi coordinate
  //
  // watch out: if the theta_cell includes theta==0 or theta==pi,
  // incude the full phi range
  if ((theta_cell == 0x1) || (theta_cell == 0x80000000))
    phi_range = 0xffffffff;
  else
    phi_range |= get_phi_cell(phi);

  return 0;
}


// test overlap
//  - r1  first range
//  - r2  second range
// return true if overlap, false otherwise.
//------------------------------------------
bool is_range_overlap(const CSphtheta_phi_range &r1, const CSphtheta_phi_range &r2){
  // check overlap in eta AND phi
  return ((r1.theta_range & r2.theta_range) && (r1.phi_range & r2.phi_range));
}

// compute union
// Note: we assume that the two intervals overlap
//  - r1  first range
//  - r2  second range
// \return union of the two ranges
//------------------------------------------
const CSphtheta_phi_range range_union (const CSphtheta_phi_range &r1, const CSphtheta_phi_range &r2){
  CSphtheta_phi_range tmp;

  // compute union in eta
  tmp.theta_range = r1.theta_range | r2.theta_range;

  // compute union in phi
  tmp.phi_range = r1.phi_range | r2.phi_range;

  return tmp;
}

}
