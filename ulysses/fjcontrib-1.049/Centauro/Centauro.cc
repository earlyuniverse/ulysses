// This code is part of Fastjet contrib

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
//----------------------------------------------------------------------x

#include "Centauro.hh"
#include "fastjet/NNH.hh"

// strings and streams
#include <sstream>
#include <limits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  //----------------------------------------------------------------------
  /// class that contains the algorithm parameters R, photon energy and photon longitudinal momentum.
  class CentauroInfo {

  public:
    CentauroInfo(double Ri, double gammaEi, double gammaPzi)
    { R_ = Ri; gammaE_ = gammaEi; gammaPz_ = gammaPzi;}

    double gammaPz() { return gammaPz_; }
    double gammaE() { return gammaE_; }
    double R() { return R_; }

  private:
    double R_, gammaE_, gammaPz_;
  };

  class CentauroBriefJet {
  public:                                                                                                                                                                                                              //For definitions see https://arxiv.org/abs/2006.10751
    // n = (1,0,0,1)
    // nbar = (1,0,0,-1)
    // P = Q/2x(1,0,0,1)                //proton
    // q = Q(0,0,0,-1)                   //virtual photon
    // etabar = -2Q/(nbar*q)*pT/(n*p)  , so in the Breit frame:  // etabar = +2*pT/(n*p)  = 2*pT/(E-pz)
    // The distance is: (for f=x)
    // dij = [(etabar_i - etabar_j)^{2} + 2*etabar_i*etabar_j(1-cos(phi_i - phi_j))]/R^{2}

    void init(const PseudoJet & jet, CentauroInfo * info) {

      R = info->R();
      gammaE = info->gammaE();
      gammaPz = info->gammaPz();

      // photon 4-momentum is q = (gammaE, 0 , 0 , gammaPz);

      double norm = 1.0/sqrt(jet.modp2());
      // pseudo-jet information needed to calculate distance
      nx = jet.px() * norm;
      ny = jet.py() * norm;
      nz = jet.pz() * norm;
      pT = jet.perp();
      phi = jet.phi();
      if(gammaE!=0 and gammaPz!=0){ //gammaE and gammaPz passed, so not running in Breit frame
        Q = sqrt(-1.0*(gammaE*gammaE-gammaPz*gammaPz));
        etabar = -2.0*(Q/(gammaE+gammaPz))*(pT/(jet.E()-jet.pz()));
      }
      else{  //gammaE and gammaPz not passed, so assume that it is running in the Breit frame
        etabar = +2.0*pT/(jet.E()-jet.pz());
      }
      // beam distance
      diB = 1.0;
    }

    double distance(const CentauroBriefJet * jet) const {

      double dij = pow(etabar - jet->etabar, 2.0) + 2*etabar*jet->etabar*(1-cos(phi- jet->phi));
      dij = dij/pow(R,2.0);

      return dij;
    }

    double beam_distance() const {
      return diB;
    }

    double pT, phi, nx, ny, nz;
    double etabar;
    double diB;
    double R, gammaE, gammaPz, Q;
  };


  std::string CentauroPlugin::description () const {
    std::ostringstream desc;
    desc << "Centauro plugin with R = " << R();
    if(gammaE()==0 and gammaPz()==0){
      desc << " gamma E and gamma Pz parameters were not given --> assume you are giving particles momenta in Breit frame";
    }
    return desc.str();
  }

  void CentauroPlugin::run_clustering(fastjet::ClusterSequence & cs) const {
    int njets = cs.jets().size();
    CentauroInfo vinfo(R(), gammaE(), gammaPz());

    NNH<CentauroBriefJet,CentauroInfo> nnh(cs.jets(),&vinfo);

    while (njets > 0) {
      int i, j, k;
      double dij = nnh.dij_min(i, j);

      if (j >= 0) {
        cs.plugin_record_ij_recombination(i, j, dij, k);
        nnh.merge_jets(i, j, cs.jets()[k], k);
      } else {

        cs.plugin_record_iB_recombination(i, dij);
        nnh.remove_jet(i);
      }

      njets--;
    }
  }


} // namespace contrib

FASTJET_END_NAMESPACE
