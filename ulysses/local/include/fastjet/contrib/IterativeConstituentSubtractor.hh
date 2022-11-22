//
// ConstituentSubtractor package
// Questions/comments: peter.berta@cern.ch
//
// Copyright (c) 2018-, Peter Berta
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

#ifndef __FASTJET_CONTRIB_ITERATIVECONSTITUENTSUBTRACTOR_HH__
#define __FASTJET_CONTRIB_ITERATIVECONSTITUENTSUBTRACTOR_HH__


#include "ConstituentSubtractor.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class IterativeConstituentSubtractor
/// A class to perform subtraction of background, e.g. pileup, from a set of input particles from whole event. The output is a vector of corrected PseudoJets. More info on the algorithm can be found in
/// https://arxiv.org/abs/1905.03470
///
/// For usage, see
/// example_iterative.cc
///

  class IterativeConstituentSubtractor : public fastjet::contrib::ConstituentSubtractor{
  public:

    ///
    /// default ctor which should be used before the event loop
    IterativeConstituentSubtractor();

    ///
    /// initialization (construction of ghosts and some checks if the provided parameters make sense). This should be used before the event loop.
    virtual void initialize();


    ///
    /// default dtor
    virtual ~IterativeConstituentSubtractor(){}
    
    ///
    /// a description of what this class does
    virtual std::string description() const;
    
    
    
    ///
    /// do iterative subtraction. The particles with |eta|>max_eta are discarded at the beginning, i.e. they are not used, nor returned. The ghosts are added automatically inside this function up to max_eta.
    virtual std::vector<fastjet::PseudoJet> subtract_event(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const *hard_proxies=0);

    ///
    /// this should be not used
    virtual std::vector<fastjet::PseudoJet> subtract_event(std::vector<fastjet::PseudoJet> const &particles, double max_eta);
    
    ///
    /// function to set the _ghost_removal member. If set to true, then the proxies with non-zero remaining pt are discarded for the next iteration
    void set_ghost_removal(bool ghost_removal);

    ///
    /// Set maximal distance and alpha parameters for the iterative CS
    void set_parameters(std::vector<double> const &max_distances, std::vector<double> const &alphas);

    void set_nearby_hard_parameters(std::vector<double> const &nearby_hard_radii, std::vector<double> const &nearby_hard_factors);


  protected:
    std::vector<double> _max_distances;
    std::vector<double> _alphas;
    std::vector<double> _nearby_hard_radii;
    std::vector<double> _nearby_hard_factors;
    bool _use_nearby_hard_iterative;
    bool _ghost_removal;
  };
  

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ITERATIVECONSTITUENTSUBTRACTOR_HH__
