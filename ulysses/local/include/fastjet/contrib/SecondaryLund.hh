// $Id: SecondaryLund.hh 1289 2021-11-09 11:53:53Z scyboz $
//
// Copyright (c) 2018-, Frederic A. Dreyer, Keith Hamilton, Alexander Karlberg,
// Gavin P. Salam, Ludovic Scyboz, Gregory Soyez, Rob Verheyen
//
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

#ifndef __FASTJET_CONTRIB_SECONDARYLUND_HH__
#define __FASTJET_CONTRIB_SECONDARYLUND_HH__

#include "LundGenerator.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{
  
//----------------------------------------------------------------------
/// \class SecondaryLund
/// Base class for definitions of the leading emission
class SecondaryLund {
public:
  /// SecondaryLund constructor
  SecondaryLund() {}

  /// destructor
  virtual ~SecondaryLund() {}

  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<LundDeclustering> & declusts) const = 0;

  int operator()(const std::vector<LundDeclustering> & declusts) const {
    return result(declusts);
  }

  /// description of the class
  virtual std::string description() const;
};
  
//----------------------------------------------------------------------
/// \class SecondaryLund_mMDT
/// Contains a definition for the leading emission using mMDTZ
class SecondaryLund_mMDT : public SecondaryLund {
public:
  /// SecondaryLund_mMDT constructor
  SecondaryLund_mMDT(double zcut = 0.025) : zcut_(zcut) {}

  /// destructor
  virtual ~SecondaryLund_mMDT() {}
  
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<LundDeclustering> & declusts) const;
  
  /// description of the class
  virtual std::string description() const;  
private:
  /// zcut parameter
  double zcut_;
};
  
//----------------------------------------------------------------------
/// \class SecondaryLund_dotmMDT
/// Contains a definition for the leading emission using dotmMDT
class SecondaryLund_dotmMDT : public SecondaryLund {
public:
  /// SecondaryLund_dotmMDT constructor
  SecondaryLund_dotmMDT(double zcut = 0.025) : zcut_(zcut) {}

  /// destructor
  virtual ~SecondaryLund_dotmMDT() {}
  
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<LundDeclustering> & declusts) const;
  
  /// description of the class
  virtual std::string description() const;  

private:
  /// zcut parameter
  double zcut_;
};
  
//----------------------------------------------------------------------
/// \class SecondaryLund_Mass
/// Contains a definition for the leading emission using mass
class SecondaryLund_Mass : public SecondaryLund {
public:
  /// SecondaryLund_Mass constructor (default mass reference is W mass)
  SecondaryLund_Mass(double ref_mass = 80.4) : mref2_(ref_mass*ref_mass) {}

  /// destructor
  virtual ~SecondaryLund_Mass() {}
  
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<LundDeclustering> & declusts) const;
  
  /// description of the class
  virtual std::string description() const;  
private:
  double mref2_;
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SECONDARYLUND_HH__

