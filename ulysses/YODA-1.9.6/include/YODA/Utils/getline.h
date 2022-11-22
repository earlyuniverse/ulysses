// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2017 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_GETLINE_H
#define YODA_GETLINE_H

#include <iostream>
#include <string>

namespace YODA {
  namespace Utils {


    // Portable version of getline taken from
    // http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
    inline std::istream& getline(std::istream& is, std::string& t) {
      t.clear();

      // The characters in the stream are read one-by-one using a std::streambuf.
      // That is faster than reading them one-by-one using the std::istream.
      // Code that uses streambuf this way must be guarded by a sentry object.
      // The sentry object performs various tasks,
      // such as thread synchronization and updating the stream state.
      std::istream::sentry se(is, true);
      std::streambuf* sb = is.rdbuf();

      for (;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
          return is;
        case '\r':
          if (sb->sgetc() == '\n')
            sb->sbumpc();
          return is;
        case EOF:
          // Also handle the case when the last line has no line ending
          if (t.empty())
            is.setstate(std::ios::eofbit);
          return is;
        default:
          t += (char)c;
        }
      }
    }


  }
}

#endif
