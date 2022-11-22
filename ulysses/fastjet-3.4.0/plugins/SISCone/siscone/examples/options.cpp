///////////////////////////////////////////////////////////////////////////////
// File: options.cpp                                                         //
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

#include "options.h"
#include <string.h>
#include <getopt.h>
#include <iostream>

using namespace std;
using namespace siscone;

#define N_DEFAULT          -1
#define R_DEFAULT          0.7
#define THRESHOLD_DEFAULT  0.5
#define PTMIN_DEFAULT      0.0
#define NPASS_DEFAULT      0
#define DEFAULT_EVENT      "events/single-event.dat"
#define SM_DEFAULT         SM_pttilde

/*******************************************
 * Coptions implementation                 *
 * options for the 'cone' sample           *
 *******************************************/

// default ctor
//--------------
Coptions::Coptions(){
  // set default flags values
  help_flag=0;
  version_flag=0;
  verbose_flag=1;

  // set default options values
  N_stop = N_DEFAULT;
  R = R_DEFAULT;
  f = THRESHOLD_DEFAULT;
  npass = NPASS_DEFAULT;
  ev_name = NULL;
  SM_var = SM_DEFAULT;
}


// default dtor
//--------------
Coptions::~Coptions(){
  if (ev_name!=NULL)
    delete[] ev_name;
}


// parse oprions
//  - argc  number of arguments from the command line
//  - argv  arguments from the command line
// return 1 on error, 0 on success
//---------------------------------
int Coptions::parse_options(int argc, char **argv){
  int opt_param;
  int option_index;
  bool stop=false;

  // browse the command-line options{
  static struct option siscone_options[]={
    // options that set a flag
    {"verbose",  no_argument, &verbose_flag, 1},
    {"quiet",    no_argument, &verbose_flag, 0},
    {"help",     no_argument, &help_flag   , 1},
    {"version",  no_argument, &version_flag, 1},
    // options setting parameters
    {"number",   required_argument, NULL, 'N'},
    {"radius",   required_argument, NULL, 'R'},
    {"fraction", required_argument, NULL, 'f'},
    {"ptmin",    required_argument, NULL, 'p'},
    {"npass",    required_argument, NULL, 'n'},
    {"event",    required_argument, NULL, 'e'},
    {"sm",       required_argument, NULL, 's'},
    {0,0,0,0}
  };

  
  do{
    // getopt_long stores the option index here.
    option_index=0;

    // retreive options
    opt_param = getopt_long(argc, argv, "hvqN:R:f:p:n:e:s:", 
			    siscone_options, &option_index);

    // Detect the end of the options.
    if (opt_param == -1)
      stop=true;

    // branch according to 'opt_param'
    switch (opt_param){
    case 'h': help_flag = 1;      break;  // help
    case 'v': verbose_flag = 1;   break;  // verbose
    case 'q': verbose_flag = 0;   break;  // quiet
    case 'N':  // max number of paprticles
      sscanf(optarg, "%d", &N_stop);
      if (N_stop<=0){
	cout << "Warning: the specified number of particles must be positive. Using default one" << endl;
	N_stop = N_DEFAULT;
      }
      break;
    case 'R':
      sscanf(optarg, "%lf", &R);
      if (R<=0){
	cout << "Warning: the specified cone radius must be positive. Using default one" << endl;
	R = R_DEFAULT;
      }      
      break;
    case 'f':
      sscanf(optarg, "%lf", &f);
      if ((f<0) || (f>1)){
	cout << "Warning: the specified split/merge threshold must be in [0,1]. Using default one" << endl;
	f = THRESHOLD_DEFAULT;
      }            
      break;
    case 'p':
      sscanf(optarg, "%lf", &ptmin);
      if (ptmin<0){
	cout << "Warning: the specified minimal pT must be non-negative. Using default one" << endl;
	ptmin = PTMIN_DEFAULT;
      }            
      break;
    case 'n':  // max number of paprticles
      sscanf(optarg, "%d", &npass);
      if (npass<0){
	cout << "Warning: the specified number of passes must be non negative. Using default one" << endl;
	npass = NPASS_DEFAULT;
      }
      break;
    case 'e':
      if (ev_name==NULL){
	ev_name = new char[strlen(optarg)+1];
	strcpy(ev_name, optarg);
      }
      break;
    case 's':
      char tmp[512];
      strcpy(tmp, optarg);
      if (strcmp(tmp, "pttilde")==0){
	SM_var = SM_pttilde;
      } else if (strcmp(tmp, "mt")==0){
	SM_var = SM_mt;
      } else if (strcmp(tmp, "pt")==0){
	SM_var = SM_pt;
      } else if (strcmp(tmp, "Et")==0){
	SM_var = SM_Et;
      } else {
	cout << "Warning: the specified varible for split--merge is not valid (should be pttilde, pt, mt or Et). Using pttilde as the default one." << endl;
	SM_var = SM_pttilde;
      }
      break;
    case 0:
    case -1:
      break;
    case '?':
      fprintf(stderr, "Giving up.\n");
      return 1;
      break;
    default:
      if (!help_flag){
	fprintf(stderr, "unrecognized option %c. Giving up.\n", opt_param);
	return 1;
      }
    }
  } while (!stop);
  
  if (ev_name==NULL){
    ev_name = new char[strlen(DEFAULT_EVENT)+1];
    strcpy(ev_name, DEFAULT_EVENT);
  }
  
  return 0;
}


// print the help message
//------------------------
int Coptions::print_help(){
  cout << siscone_package_name() << " " << siscone_version() << endl;
  cout << "Usage: " << siscone_package_name() << " <args>" << endl;
  cout << endl;
  cout << "Here is an exhaustive list of the arguments:" << endl;
  cout << "Parameters control (with default values):" << endl;
  cout << "  -n <val>, --number=<val>  : set the maximum number of particles allowed (all)" << endl;
  cout << "  -R <val>, --radius=<val>  : set the radius (" << R_DEFAULT << ")" << endl;
  cout << "  -f <val>, --fraction=<val>: set the overlap parameter (" << THRESHOLD_DEFAULT << ")" << endl;
  cout << "  -p <val>, --ptmin=<val>   : set the minimal pT for protojets (" << PTMIN_DEFAULT << ")" << endl;
  cout << "  -n <val>, --npass=<val>   : set the maximal number of passes (0 for no limit) (" << NPASS_DEFAULT << ")" << endl;
  cout << "  -e <val>, --event=<val>   : set the event filename (" << DEFAULT_EVENT << ")" << endl;
  cout << "  -s <val>, --sm=<val>      : variable for split--merge: pttilde, mt, pt or Et (pttilde)" << endl;
  cout << endl;
  cout << "Output flags" << endl;
  cout << "  --version    : show version information" << endl;
  cout << "  -h, --help   : show this message" << endl;
  cout << "  -v, --verbose: be verbose (on by default)" << endl;
  cout << "  -q, --quiet  : be quiet" << endl;
  cout << endl;

  return 0;
}


// print program version
//-----------------------
int Coptions::print_version(){
  cout << siscone_package_name() << " " << siscone_version() << endl;
  cout << "Copyright (C) 2006." << endl;
  cout << siscone_package_name() << " comes with NO WARRANTY," << endl;
  cout << "to the extent permitted by law." << endl;
  cout << "You may redistribute copies of " << siscone_package_name() << endl;
  cout << "under the terms of the GNU General Public License." << endl;
  cout << "For more information about these matters," << endl;
  cout << "see the files named COPYING." << endl;
  cout << "Please send bugs or comments to AUTHORS" << endl;

  return 0;
}
