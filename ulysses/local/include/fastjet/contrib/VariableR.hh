//  VariableR Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2009-2016
//  David Krohn, Gregory Soyez, Jesse Thaler, and Lian-Tao Wang
//
//  $Id: VariableR.hh 903 2016-03-09 16:31:58Z jthaler $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_VARIABLER_HH__
#define __FASTJET_CONTRIB_VARIABLER_HH__

#include <fastjet/internal/base.hh>

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "VariableRPlugin.hh"
#include <map>
#include <queue>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
    
   ////////
   //
   //  For backwards compatibility with v1.0, keep old constructors
   //  These just call the VR constructor with certain parameter choice
   //
   ////////
   
   class AKTVR : public VariableRPlugin {
   public:
      AKTVR (double rho, double max_r) : VariableRPlugin(rho,0.0,max_r,AKTLIKE) {}
   };

   class CAVR : public VariableRPlugin {
   public:
      CAVR (double rho, double max_r) : VariableRPlugin(rho,0.0,max_r,CALIKE) {}
   };

   class KTVR : public VariableRPlugin {
   public:
      KTVR (double rho, double max_r) : VariableRPlugin(rho,0.0,max_r,KTLIKE) {}
   };
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VARIABLER_HH__
