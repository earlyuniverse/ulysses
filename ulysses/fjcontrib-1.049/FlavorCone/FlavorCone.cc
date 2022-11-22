// $Id: FlavorCone.cc 1056 2017-09-06 20:59:06Z philten $
//
// Copyright (c) 2017, Philip Ilten, Nicholas Rodd, Jesse Thaler, 
// and Michael Williams
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

#include "FlavorCone.hh"
#include <sstream>
#include <limits>

FASTJET_BEGIN_NAMESPACE // defined in fastjet/internal/base.hh

using namespace std;

namespace contrib {

// Main constructor just sets seeds and r value.
FlavorConePlugin::FlavorConePlugin(const std::vector<PseudoJet> &seeds, 
				   double rcut) : _rcut(rcut) {
  for (unsigned int iseed = 0; iseed < seeds.size(); ++iseed)
    _seeds.push_back(seeds[iseed]);
}

// Copy constructor.
FlavorConePlugin::FlavorConePlugin(const FlavorConePlugin & plugin) {
  *this = plugin;
}

// Return a description for the plugin.
string FlavorConePlugin::description() const {
  stringstream desc;
  desc << "FlavorCone plugin with " << _seeds.size() << " seeds and rcut = " 
       << _rcut;
  return desc.str();
}

// Run the jet clustering.
void FlavorConePlugin::run_clustering(ClusterSequence &cs) const {
  
  // Create extras for later storage of jet/seed mapping.
  FlavorConePlugin::Extras *extras(new FlavorConePlugin::Extras());
  int nprts(cs.jets().size());
  int nseeds(_seeds.size());
  int ijet = 0;
  vector<int> jets(nseeds, -1);

  // Combine the particles.
  for (int iprt = 0; iprt < nprts; ++iprt) {
    int imin = -1;
    double drmin(numeric_limits<double>::infinity());

    // Find nearest seed.
    for (int iseed = 0; iseed < nseeds; ++iseed) {
      double dr(cs.jets()[iprt].squared_distance(_seeds[iseed]));
      if (dr < drmin) {
        drmin = dr;
        imin = iseed;
      }
    }
    
    // If no seed closer than R_cut, abandon.
    if (drmin > _rcut*_rcut) continue;
    
    
    if (jets[imin] == -1) {
      // If no jet associated to that seed, make a new jet.
      jets[imin] = iprt;
    } else {
      // Otherwise assign this particle to exisiting jet.
      cs.plugin_record_ij_recombination(iprt, jets[imin], drmin, ijet);
      jets[imin] = ijet;
    }
  }

  // Flag jets for storage and associate seeds with extras.
  for (int iseed = nseeds - 1; iseed >= 0; --iseed) 
    if (jets[iseed] != -1) {
      cs.plugin_record_iB_recombination(jets[iseed], iseed);
      extras->_seeds[jets[iseed]] = _seeds[iseed];
    }
  cs.plugin_associate_extras(extras);
}

// Return the R-cut value.
double FlavorConePlugin::R() const {return _rcut;}
  
// Return if exclusive sequence is meaningful.
// It is not, since we do not change input particle order.
bool FlavorConePlugin::exclusive_sequence_meaningful() const {return false;}

// Return if for e+ e-.
// Currently, this only works for pp.
bool FlavorConePlugin::is_spherical() const {return false;}

// Main constructor for the FlavorCone Extras class.
LimitedWarning FlavorConePlugin::Extras::_warn_seed;
FlavorConePlugin::Extras::Extras() : _invalid_seed(0, 0, 0, -1) {;}

// Return the seed associated with a jet.
const PseudoJet &FlavorConePlugin::Extras::seed(const PseudoJet &jet) const {
  std::map<int, PseudoJet>::const_iterator itr = 
    _seeds.find(jet.cluster_hist_index());
  if (itr == _seeds.end()) {
    _warn_seed.warn("FlavorConePlugin::Extras::seed: No seed associated with "
		    "this jet, invalid seed with momentum (0, 0, 0, -1) "
		    "returned.");
    return _invalid_seed;
  } else return itr->second;
}

} // namespace contrib

FASTJET_END_NAMESPACE
