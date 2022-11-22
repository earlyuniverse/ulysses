// $Id: SimpleGhostRescaler.hh 859 2015-09-21 10:11:32Z gsalam $
//
// Copyright (c) 2012-, Matteo Cacciari, Jihun Kim, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_CONTRIB_SIMPLE_GHOST_RESCALER_HH__
#define __FASTJET_CONTRIB_SIMPLE_GHOST_RESCALER_HH__

#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/Error.hh"

FASTJET_BEGIN_NAMESPACE

namespace contrib{

// TODO: Things to think about for the future of this class:
//
// maybe we can have this class directly inside the GenericSubtraction
// one and cache a few things (like recombiner...)  once and for all
// for each jet. In that case, the rescaler would either be created
// with a jet (or have a "set_jet") --- we could even cache the
// constituents and potentially where the ghosts are --- and then the
// () call would take the scales as arguments.


/// \class SimpleGhostRescaler
/// construct, from an original jet, a composite jet where all the
/// ghosts have been rescaled
class SimpleGhostRescaler : public FunctionOfPseudoJet<PseudoJet>{
public:

  /// ctor with the following parameters
  ///  \param pt_scale     the 
  ///  \param mdelta_scale
  ///  \param ghost_scale
  SimpleGhostRescaler(double pt_scale, double mdelta_scale, double ghost_scale)
    : _pt_scale(pt_scale), _mdelta_scale(mdelta_scale), _ghost_scale(ghost_scale){};

  virtual PseudoJet result(const PseudoJet &jet) const{
    // directly remove the case os a composite jet
    if (!jet.has_associated_cluster_sequence()){
      // just handle the case of a composite jet here
      if (!jet.has_pieces())
	throw Error("Ghost rescaling can only be performed on jets with an associated ClusterSequence or composite jets (with pieces associated with a Clustersequence)");
      
      return join((*this)(jet.pieces()));
    }

    // make sure that the jet has an area
    if (!jet.has_area()){
      throw Error("Ghost rescaling can only be applied on jets with an area");
    }

    if (!jet.validated_csab()->has_explicit_ghosts()){
      throw Error("Ghost rescaling can only be applied on jets with explicit ghosts");
    }

    Selector ghost_selector = SelectorIsPureGhost();
    std::vector<PseudoJet> rescaled_constituents = jet.constituents();

    // loop over the constituents, rescale the ghosts and do not touch
    // the regular particles
    for (unsigned int i=0;i<rescaled_constituents.size(); i++){
      if (ghost_selector.pass(rescaled_constituents[i])){
	double pt = _pt_scale*(rescaled_constituents[i].perp()/_ghost_scale);

	// we want to choose a mass such that sqrt(m^2 + pt^2)-pt = _mdelta_scale.
	// m^2 = (pt+_mdelta_scale)^2 - pt^2 = mdelta_scale^2 + 2 pt*_mdelta_scale
	double m = sqrt(_mdelta_scale*(_mdelta_scale+2*pt));

	// Note that the line below not only updates the
	// momentum, it also kills all structure information. This is
	// mandatory because asking for a composite jet's constituents
	// recurses into its pieces!
	rescaled_constituents[i] = PtYPhiM(pt, rescaled_constituents[i].rap(), rescaled_constituents[i].phi(), m);
      }
    }

    return join(rescaled_constituents);
  }

private:
  double _pt_scale, _mdelta_scale, _ghost_scale;
};

}

FASTJET_END_NAMESPACE

#endif // __SIMPLE_GHOST_RESCALER_HH__
