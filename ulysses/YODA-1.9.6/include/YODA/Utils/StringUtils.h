// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_STRINGUTILS_H
#define YODA_STRINGUTILS_H

#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>
//#include <iostream>

namespace YODA {
  namespace Utils {


    /// Exception to be thrown by lexical_cast below
    struct bad_lexical_cast : public std::runtime_error {
      bad_lexical_cast(const std::string& what) : std::runtime_error(what) {}
    };

    /// Convert between any types via stringstream
    template<typename T, typename U>
    T lexical_cast(const U& in) {
      try {
        std::stringstream ss;
        ss << in;
        T out;
        ss >> out;
        return out;
      } catch (const std::exception& e) {
        throw bad_lexical_cast(e.what());
      }
    }


    /// Generic convenient conversion to string
    template <typename T>
    inline std::string toStr(const T& x) {
      // std::ostringstream ss; ss << x;
      // return ss.str();
      return lexical_cast<std::string>(x);
    }


    /// Convert a string to lower-case
    inline std::string toLower(const std::string& s) {
      std::string out = s;
      std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) tolower);
      return out;
    }


    /// Convert a string to upper-case
    inline std::string toUpper(const std::string& s) {
      std::string out = s;
      std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) toupper);
      return out;
    }


    /// Replace XML reserved characters with &XXXX; encodings
    inline std::string encodeForXML(const std::string& in) {
      std::string out = in;
      typedef std::pair<std::string, std::string> CharsToEntities;
      std::vector<CharsToEntities> cs2es;
      cs2es.push_back(std::make_pair("&", "&amp;"));
      cs2es.push_back(std::make_pair("<", "&lt;"));
      cs2es.push_back(std::make_pair(">", "&gt;"));
      for (std::vector<CharsToEntities>::const_iterator c2e = cs2es.begin(); c2e != cs2es.end(); ++c2e) {
        std::string::size_type pos = -1;
        while ( ( pos = out.find(c2e->first, pos + 1) ) != std::string::npos ) {
          out.replace(pos, 1, c2e->second);
        }
      }
      return out;
    }


    /// Does the string @a s contain the given substring?
    inline bool contains(const std::string& s, const std::string& sub) {
      return s.find(sub) != std::string::npos;
    }

    /// Does the string @a s start with the given substring?
    inline bool startswith(const std::string& s, const std::string& sub) {
      return s.find(sub) == 0;
    }

    /// Does the string @a s end with the given substring?
    inline bool endswith(const std::string& s, const std::string& sub) {
      const size_t pos = s.find(sub);
      return pos != std::string::npos && pos == (s.size()-sub.size());
    }


    /// In-place trim from start
    inline std::string& iltrim(std::string& s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char c){ return !std::isspace(c); }));

      return s;
    }

    /// In-place trim from end
    inline std::string& irtrim(std::string& s) {
      s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c){ return !std::isspace(c); }).base(), s.end());
      return s;
    }

    /// In-place trim from both ends
    inline std::string& itrim(std::string& s) {
      return iltrim(irtrim(s));
    }


    /// Trim from start
    inline std::string ltrim(const std::string& s) { //< @note Could just pass by value, but const& sends a symbolic message
      std::string s2 = s;
      return iltrim(s2);
    }

    /// Trim from end
    inline std::string rtrim(const std::string& s) { //< @note Could just pass by value, but const& sends a symbolic message
      std::string s2 = s;
      return irtrim(s2);
    }

    /// Trim from both ends
    inline std::string trim(const std::string& s) { //< @note Could just pass by value, but const& sends a symbolic message
      std::string s2 = s;
      return itrim(s2);
    }


  }


  // Annoyingly, this doesn't work without significant SFINAE etc., due to clashes with existing string & stream operators
  // /// string + any appending operator for convenience in the YODA namespace
  // template <typename T>
  // std::string& operator + (std::string& s, const T& x) {
  //   s += Utils::toStr(x);
  //   return s;
  // }
  // template <typename T>
  // std::string operator + (const T& x, std::string& s) {
  //   return Utils::toStr(x) + s;
  // }
  // template <typename T>
  // std::string operator + (const char* s, const T& x) {
  //   return s + Utils::toStr(x);
  // }
  // template <typename T>
  // std::string operator + (const T& x, const char* s) {
  //   return Utils::toStr(x) + x;
  // }


}

#endif
