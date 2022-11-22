// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2020 The YODA collaboration (see AUTHORS for details)
//
#ifndef YODA_Fillable_h
#define YODA_Fillable_h

#include "YODA/AnalysisObject.h"
#include "YODA/Exceptions.h"

namespace YODA {


  /// A base class for all fillable objects
  class Fillable {
  public:

    /// @name Constructors
    //@{

    /// Virtual destructor for inheritance
    virtual ~Fillable() = default;


    /// Fill-dimension of this data object
    virtual size_t fillDim() const = 0;


    /// @name Modifiers
    //@{

    /// Reset
    virtual void reset() = 0;

    /// Rescale as if all fill weights had been different by factor @a scalefactor
    virtual void scaleW(double scalefactor) = 0;
    // setAnnotation("ScaledBy", annotation<double>("ScaledBy", 1.0) * scalefactor);


    /// @name Whole histo data
    //@{

    /// Get the number of fills
    virtual double numEntries(bool includeoverflows=true) const = 0;

    /// Get the effective number of fills
    virtual double effNumEntries(bool includeoverflows=true) const = 0;

    /// Get sum of weights in histo
    virtual double sumW(bool includeoverflows=true) const = 0;

    /// Get sum of squared weights in histo
    virtual double sumW2(bool includeoverflows=true) const = 0;

    //@}

  };


}

#endif
