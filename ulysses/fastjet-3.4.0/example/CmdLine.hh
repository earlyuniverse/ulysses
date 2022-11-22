///////////////////////////////////////////////////////////////////////////////
// File: CmdLine.hh                                                          //
// Part of the CmdLine library
//                                                                           //
// Copyright (c) 2007 Gavin Salam                                            //
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


#ifndef __CMDLINE__
#define __CMDLINE__

#include<string>
#include<sstream>
#include<map>
#include<vector>
using namespace std;

/// \if internal_doc
/// @ingroup internal
/// \class CmdLine
/// Class designed to deal with command-line arguments in a fashion similar
/// to what was done in f90 iolib.
///
/// Note that functionality might be slightly different? 
/// Currently do not implement access to arguments by index
/// though data structure would in principle allow this quite easily.
///
/// GPS 03/01/05
/// [NB: wonder if some of this might be more efficiently written 
/// with templates for different type that can be read in...]
///
/// Other question: dealing with list of options is rather common
/// occurrence -- command-line arguments, but also card files; maybe one
/// could somehow use base/derived classes to share common functionality? 
/// \endif
class CmdLine {
  mutable map<string,int> __options;
  vector<string> __arguments;
  mutable map<string,bool> __options_used;
  //string __progname;
  string __command_line;

 public :
  CmdLine() {};
  /// initialise a CmdLine from a C-style array of command-line arguments
  CmdLine(const int argc, char** argv);
  /// initialise a CmdLine from a C++ vector of arguments 
  CmdLine(const vector<string> & args);

  /// true if the option is present
  bool    present(const string & opt) const;
  /// true if the option is present and corresponds to a value
  bool    present_and_set(const string & opt) const;

  /// return a reference to the vector of command-line arguments (0 is
  /// command).
  inline const vector<string> & arguments() const {return __arguments;}

  /// returns the value of the argument converted to type T
  template<class T> T value(const string & opt) const;
  template<class T> T value(const string & opt, const T & defval) const;


  /// return the integer value corresponding to the given option
  int     int_val(const string & opt) const;
  /// return the integer value corresponding to the given option or default if option is absent
  int     int_val(const string & opt, const int & defval) const;

  /// return the double value corresponding to the given option
  double  double_val(const string & opt) const;
  /// return the double value corresponding to the given option or default if option is absent
  double  double_val(const string & opt, const double & defval) const;

  /// return the string value corresponding to the given option
  string  string_val(const string & opt) const;
  /// return the string value corresponding to the given option or default if option is absent
  string  string_val(const string & opt, const string & defval) const;

  /// return the full command line
  string  command_line() const;

  /// return true if all options have been asked for at some point or other
  bool all_options_used() const;

 private:
  /// builds the internal structures needed to keep track of arguments and options
  void init();

  /// report failure of conversion
  void _report_conversion_failure(const string & opt, 
                                  const string & optstring) const;


};



/// returns the value of the argument converted to type T
template<class T> T CmdLine::value(const string & opt) const {
  T result;
  string optstring = string_val(opt);
  istringstream optstream(optstring);
  optstream >> result;
  if (optstream.fail()) _report_conversion_failure(opt, optstring);
  return result;
}

/// for the string case, just copy the string...
template<> inline string CmdLine::value<string>(const string & opt) const {
  return string_val(opt);}



template<class T> T CmdLine::value(const string & opt, const T & defval) const {
  if (this->present_and_set(opt)) {return value<T>(opt);} 
  else {return defval;}
}

#endif
