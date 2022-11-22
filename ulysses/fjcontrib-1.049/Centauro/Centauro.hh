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

#ifndef __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
#define __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  class CentauroPlugin : public JetDefinition::Plugin {
  public:

    /// Constructor for the Centauro Plugin class.

    /// Three floating point arguments are specified to set the parameters
    /// the radius parameter R, and gammaE and gammaPz are the energy and pz of virtual photon

    CentauroPlugin (double R, double gammaE, double gammaPz) : _R(R), _gammaE(gammaE), _gammaPz(gammaPz){}

    /// if only one argument is passed, assume that it runs in Breit frame so gammaE and gammaPz info not required so set
    /// so they are set to zero
    
    CentauroPlugin (double R) : _R(R), _gammaE(0), _gammaPz(0){}
    /// copy constructor
    CentauroPlugin (const CentauroPlugin & plugin) {
      *this = plugin;
    }

    // the things that are required by base class
    virtual std::string description () const;
    virtual void run_clustering(ClusterSequence &) const;


    virtual double R() const {return _R;}
    virtual double gammaE() const {return _gammaE;}
    virtual double gammaPz() const {return _gammaPz;}



    /// avoid the warning whenever the user requests "exclusive" jets
    /// from the cluster sequence
    virtual bool exclusive_sequence_meaningful() const {return true;}

  private:
    double _R;
    double _gammaE;
    double _gammaPz;
  };

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_CENTAUROJETALGORITHM_HH__
