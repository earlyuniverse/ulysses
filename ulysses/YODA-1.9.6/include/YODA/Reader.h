// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_READER_H
#define YODA_READER_H

#include "YODA/AnalysisObject.h"
#include "YODA/Utils/Traits.h"
#include "YODA/Index.h"
#include <string>
#include <fstream>
#include <vector>
#include <type_traits>
#include <iostream>
#include <sstream>

namespace YODA {


  /// Pure virtual base class for various output writers.
  class Reader {
  public:

    /// Virtual destructor
    virtual ~Reader() {}


    /// @name Reading multiple analysis objects,
    /// @{

    /// @brief Read in a collection of objects @a objs from output stream @a stream.
    ///
    /// This version fills (actually, appends to) a variable supplied container
    /// Note: SFINAE is used to check for a void push_back(const AnalysisObject*) method
    ///
    /// @todo Extend SFINAE Pushable cf. Writer to allow adding to containers of smart ptr type
    template<typename CONT>
    typename std::enable_if<YODA::Pushable<CONT,AnalysisObject*>::value>::type
    read(std::istream& stream, CONT& aos) {
      // if CONT==std::vector<AnalysisObject*>, the compiler should select
      // the virtual method below, since it prefers non-templated methods in the lookup
      // otherwise we would enter a endless recursion. Check in case of problems.
      std::vector<AnalysisObject*> v_aos;
      read(stream, v_aos);
      for (const AnalysisObject* ao : v_aos) aos.push_back(ao);
    }

    /// @brief Read in a collection of objects @a objs from output stream @a stream.
    ///
    /// This version fills (actually, appends to) a supplied vector, avoiding copying,
    /// and is hence CPU efficient.
    ///
    virtual void read(std::istream& stream, std::vector<AnalysisObject*>& aos) = 0;

    /// @brief Read in a collection of objects from output stream @a stream.
    ///
    /// This version returns a vector by value, involving copying, and is hence less
    /// CPU efficient than the alternative version where a vector is filled by reference.
    std::vector<AnalysisObject*> read(std::istream& stream) {
      std::vector<AnalysisObject*> rtn;
      read(stream, rtn);
      return rtn;
    }


    /// @brief Read in a collection of objects @a objs from file @a filename.
    ///
    ///
    /// This version fills (actually, appends to) a variable supplied container
    /// Note: SFINAE is used to check for a void push_back(const AnalysisObject*) method
    ///
    /// @todo Extend SFINAE Pushable cf. Writer to allow adding to containers of smart ptr type
    template<typename CONT>
    typename std::enable_if<YODA::Pushable<CONT,AnalysisObject*>::value>::type
    read(const std::string& filename, CONT& aos) {
      // if CONT==std::vector<AnalysisObject*>, the compiler should select
      // the virtual method below, since it prefers non-templated methods in the lookup
      // otherwise we would enter a endless recursion. Check in case of problems.
      std::vector<AnalysisObject*> v_aos;
      read(filename, v_aos);
      for (const AnalysisObject* ao : v_aos) aos.push_back(ao);
    }

    /// @brief Read in a collection of objects @a objs from file @a filename.
    ///
    /// This version fills (actually, appends to) a supplied vector, avoiding copying,
    /// and is hence CPU efficient.
    ///
    void read(const std::string& filename, std::vector<AnalysisObject*>& aos) {
      if (filename != "-") {
        try {
          std::ifstream instream;
          instream.open(filename.c_str());
          if (instream.fail())
            throw ReadError("Reading from filename " + filename + " failed");
          read(instream, aos);
          instream.close();
        } catch (std::ifstream::failure& e) {
          throw ReadError("Reading from filename " + filename + " failed: " + e.what());
        }
      } else {
        try {
          read(std::cin, aos);
        } catch (std::runtime_error& e) {
          throw ReadError("Reading from stdin failed: " + std::string(e.what()));
        }
      }
    }

    /// @brief Read in a collection of objects from output stream @a stream.
    ///
    /// This version returns a vector by value, involving copying, and is hence less
    /// CPU efficient than the alternative version where a vector is filled by reference.
    std::vector<AnalysisObject*> read(const std::string& filename) {
      std::vector<AnalysisObject*> rtn;
      read(filename, rtn);
      return rtn;
    }

    /// @}


    /// @brief Make file index
    ///
    /// Makes an index of a file's contents.
    ///
    /// @param[in] filename Path to the file to index.
    /// @return @sa Index
    Index mkIndex(const std::string& filename) {
      if (filename != "-") {
        try {
          std::ifstream instream;
          instream.open(filename.c_str());
          Index idx(mkIndex(instream));
          instream.close();
          return idx;
        } catch (std::ifstream::failure& e) {
          throw WriteError("Reading of filename " + filename +
                           " failed: " + e.what());
        }
      } else {
        try {
          return mkIndex(std::cin);
        } catch (std::runtime_error& e) {
          throw ReadError("Reading from stdin failed: " + std::string(e.what()));
        }
      }
    }

    /// @brief Make stream index
    ///
    /// Makes an index of an input stream's contents.
    ///
    /// @param[in] stream Input stream to index.
    /// @return @sa Index
    virtual Index mkIndex(std::istream& stream) = 0;
  };


  /// Factory function to make a writer object by format name or a filename
  Reader& mkReader(const std::string& format_name);


}

#endif
