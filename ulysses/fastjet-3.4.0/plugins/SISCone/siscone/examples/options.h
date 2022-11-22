///////////////////////////////////////////////////////////////////////////////
// File: options.h                                                           //
// Description: management of the cmdline options of the main program        //
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

#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include "siscone/siscone.h"

/**
 * \class Coptions

 * options for the 'cone' sample
 */
class Coptions{
 public:
  /// default ctor
  Coptions();

  /// default dtor
  ~Coptions();

  /// parse oprions
  /// \param argc  number of arguments from the command line
  /// \param argv  arguments from the command line
  /// \return 1 on error, 0 on success
  int parse_options(int argc, char **argv);

  /// print the help message
  int print_help();

  /// print program version
  int print_version();

  // flags
  int help_flag;     ///< do we need to print the help message
  int version_flag;  ///< do we need to print the version description
  int verbose_flag;  ///< do we need to print the help message

  // options
  int N_stop;        ///< maximum number of particle
  double R;          ///< cone radius
  double f;          ///< split/merge threshold
  double ptmin;      ///< minimal pT for jet candidates
  char *ev_name;     ///< event to read
  int npass;         ///< number of passes (0 for \infty)

  /// variable for split-merge
  siscone::Esplit_merge_scale SM_var;
};

#endif
