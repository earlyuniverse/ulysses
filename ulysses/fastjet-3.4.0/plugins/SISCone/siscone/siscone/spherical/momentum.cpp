///////////////////////////////////////////////////////////////////////////////
// File: momentum.cpp                                                        //
// Description: source file for 4-momentum class Cmomentum                   //
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

#include "momentum.h"
#include <math.h>
#include <stdlib.h>

namespace siscone_spherical{

/*************************************************************************
 * class CSph3vector                                                     *
 * This class contains the information for particle or group of          *
 * particles management.                                                 *
 *************************************************************************/
 
// default ctor
//--------------
CSph3vector::CSph3vector(){
  _theta = _phi = _norm = 0.0;
  px = py = pz = 0.0;
  ref = siscone::Creference();
}

// ctor with initialisation
//--------------------------
CSph3vector::CSph3vector(double _px, double _py, double _pz){
  px = _px;
  py = _py;
  pz = _pz;

  // compute the norm
  build_norm();

  ref = siscone::Creference();
}

// default dtor
//--------------
CSph3vector::~CSph3vector(){

}


// assignment of vectors
//-----------------------
CSph3vector& CSph3vector::operator = (const CSph3vector &v){
  px = v.px;
  py = v.py;
  pz = v.pz;

  _norm  = v._norm;
  _theta = v._theta;
  _phi   = v._phi;

  ref = v.ref;
  return *this;
}

// addition of vectors
//------------------------------------------------
const CSph3vector CSph3vector::operator + (const CSph3vector &v){
  CSph3vector tmp = *this;
  return tmp+=v;
}

// subtraction of vectors
//------------------------------------------------
const CSph3vector CSph3vector::operator - (const CSph3vector &v){
  CSph3vector tmp = *this;
  return tmp-=v;
}

// division by constant
//------------------------------------------------
const CSph3vector CSph3vector::operator / (const double &r){
  CSph3vector tmp = *this;
  return tmp/=r;
}

// incrementation
//------------------------------------------------
CSph3vector& CSph3vector::operator += (const CSph3vector &v){
  px+=v.px;
  py+=v.py;
  pz+=v.pz;

  return *this;
}

// decrementation
//------------------------------------------------
CSph3vector& CSph3vector::operator -= (const CSph3vector &v){
  px-=v.px;
  py-=v.py;
  pz-=v.pz;

  return *this;
}

// multiplication by a constant
//------------------------------------------------
CSph3vector& CSph3vector::operator *= (const double &r){
  px*=r;
  py*=r;
  pz*=r;

  return *this;
}

// division by a constant
//------------------------------------------------
CSph3vector& CSph3vector::operator /= (const double &r){
  px/=r;
  py/=r;
  pz/=r;

  _norm/=r;

  return *this;
}

// build norm from 3-momentum info
void CSph3vector::build_norm(){
  _norm = norm();
}

// build norm from 3-momentum info
void CSph3vector::build_thetaphi(){
  _theta = theta();
  _phi = phi();
}


// for this direction, compute the two reference directions
// used to measure angles
void CSph3vector::get_angular_directions(CSph3vector &angular_dir1, CSph3vector &angular_dir2){
  if (px < py){
    if (pz < px){
      // z smallest
      angular_dir1 = CSph3vector(-py, px, 0.0);
    } else {
      // x smallest
      angular_dir1 = CSph3vector(0.0, -pz, py);
    }
  } else {
    if (pz < py){
      // z smallest
      angular_dir1 = CSph3vector(-py, px, 0.0);
    } else {
      // y smallest
      angular_dir1 = CSph3vector(-pz, 0.0, px);
    }
  }
  angular_dir2 = cross_product3(*this, angular_dir1);
  // We'll simply take x & y so the reflection symmetry is not broken
  //angular_dir1 = CSph3vector(0.0, -pz, py);
  //angular_dir2 = CSph3vector(-pz, 0.0, -px);
}

/*************************************************************************
 * class CSphmomentum                                                    *
 * This class contains the information for particle or group of          *
 * particles management.                                                 *
 * It includes all Lorentz properties as well as tools for summing them. *
 *************************************************************************/
 
// default ctor
//--------------
CSphmomentum::CSphmomentum(){
  E=0.0;
  index = -1;
}

// ctor with initialisation
//--------------------------
CSphmomentum::CSphmomentum(double _px, double _py, double _pz, double _E)
  : CSph3vector(_px, _py, _pz) {
  E  = _E;

  // compute the angles
  build_thetaphi();
}

// ctor with initialisation
//--------------------------
CSphmomentum::CSphmomentum(CSph3vector &_v, double _E)
  : CSph3vector(_v.px, _v.py, _v.pz) {
  E  = _E;
}

// default dtor
//--------------
CSphmomentum::~CSphmomentum(){

}

// assignment of vectors
//-----------------------
CSphmomentum& CSphmomentum::operator = (const CSphmomentum &v){
  px = v.px;
  py = v.py;
  pz = v.pz;
  E  = v.E;

  _norm  = v._norm;
  _theta = v._theta;
  _phi   = v._phi;

  ref = v.ref;
  return *this;
}

// addition of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
const CSphmomentum CSphmomentum::operator + (const CSphmomentum &v){
  CSphmomentum tmp = *this;
  return tmp+=v;
}

// incrementation of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
CSphmomentum& CSphmomentum::operator += (const CSphmomentum &v){
  px+=v.px;
  py+=v.py;
  pz+=v.pz;
  E +=v.E;

  ref+=v.ref;

  return *this;
}

// decrementation of vectors
// !!! WARNING !!! no updating of eta and phi !!!
//------------------------------------------------
CSphmomentum& CSphmomentum::operator -= (const CSphmomentum &v){
  px-=v.px;
  py-=v.py;
  pz-=v.pz;
  E -=v.E;

  ref-=v.ref;
  return *this;
}


// ordering of two vectors
// the default ordering is w.r.t. their references
//-------------------------------------------------
bool operator < (const CSphmomentum &v1, const CSphmomentum &v2){
  return v1.ref < v2.ref;
}

// ordering of vectors in eta (e.g. used in collinear tests)
//-----------------------------------------------------------
bool momentum_theta_less(const CSphmomentum &v1, const CSphmomentum &v2){
  return v1._theta < v2._theta;
}

// ordering of vectors in pt
//---------------------------
bool momentum_pt_less(const CSphmomentum &v1, const CSphmomentum &v2){
  return v1.perp2() < v2.perp2();
}

}

