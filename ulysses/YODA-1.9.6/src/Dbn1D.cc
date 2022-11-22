// -*- C++ -*-
//
// This file is part of YODA -- Yet more Objects for Data Analysis
// Copyright (C) 2008-2021 The YODA collaboration (see AUTHORS for details)
//
#include "YODA/Dbn1D.h"

namespace YODA {


  double Dbn1D::xMean() const {
    if (effNumEntries() == 0 || sumW() == 0) {
      throw LowStatsError("Requested mean of a distribution with no net fill weights");
    }
    // This is ok, even for negative sum(w)
    return sumWX()/sumW();
  }


  double Dbn1D::xVariance() const {
    // Weighted variance defined as
    // sig2 = ( sum(wx**2) * sum(w) - sum(wx)**2 ) / ( sum(w)**2 - sum(w**2) )
    // see http://en.wikipedia.org/wiki/Weighted_mean
    if (effNumEntries() == 0) {
      throw LowStatsError("Requested variance of a distribution with no net fill weights");
    } else if (effNumEntries() <= 1.0) { //fuzzyLessEquals(effNumEntries(), 1.0)) {
      throw LowStatsError("Requested variance of a distribution with <= 1 effective entry");
    }
    const double num = sumWX2()*sumW() - sqr(sumWX());
    const double den = sqr(sumW()) - sumW2();
    if (den==0.) {
      throw WeightError("Undefined weighted variance");
    }
    /// @todo Isn't this sensitive to the overall scale of the weights?
    /// Shouldn't it check if den is bigger then num by a set number of
    /// orders of magnitude and vice versa?
    // if (fabs(num) < 1e-10 && fabs(den) < 1e-10) {
    //   throw WeightError("Numerically unstable weights in width calculation");
    // }
    // The weighted variance as defined above
    const double var = num/den;
    /// We take the modulus of the weighted variance since the expression above can be negative with weighted means
    /// @todo Is this the correct approach? There is no information online other than "weights are non-negative"...
    return fabs(var);
  }


  double Dbn1D::xStdErr() const {
    // Handle zero/negative sum weight
    if (effNumEntries() == 0) {
      throw LowStatsError("Requested std error of a distribution with no net fill weights");
    }
    /// @todo Unbiased should check that Neff > 1 and divide by N-1?
    return std::sqrt(xVariance() / effNumEntries());
  }


  double Dbn1D::xRMS() const {
    // Weighted RMS defined as
    // rms = sqrt(sum{w x^2} / sum{w})
    if (effNumEntries() == 0) {
      throw LowStatsError("Requested RMS of a distribution with no net fill weights");
    }
    const double meansq = sumWX2() / sumW();
    return std::sqrt(meansq);
  }


  Dbn1D& Dbn1D::add(const Dbn1D& d) {
    _dbnW     += d._dbnW;
    _sumWX    += d._sumWX;
    _sumWX2   += d._sumWX2;
    return *this;
  }


  Dbn1D& Dbn1D::subtract(const Dbn1D& d) {
    _dbnW     -= d._dbnW;
    _sumWX    -= d._sumWX;
    _sumWX2   -= d._sumWX2;
    return *this;
  }


}
