// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Dbn0D.h"
#include <cmath>

namespace YODA {


  double Dbn0D::errW() const {
    return sqrt(sumW2());
  }

  double Dbn0D::relErrW() const {
    if (effNumEntries() == 0 || sumW() == 0) {
      throw LowStatsError("Requested relative error of a distribution with no net fill weights");
    }
    return errW()/sumW();
  }


  Dbn0D& Dbn0D::add(const Dbn0D& d) {
    _numEntries += d._numEntries;
    _sumW     += d._sumW;
    _sumW2    += d._sumW2;
    return *this;
  }

  Dbn0D& Dbn0D::subtract(const Dbn0D& d) {
    _numEntries += d._numEntries; //< @todo Hmm, add or subtract?!?
    _sumW     -= d._sumW;
    _sumW2    += d._sumW2;
    return *this;
  }


}
