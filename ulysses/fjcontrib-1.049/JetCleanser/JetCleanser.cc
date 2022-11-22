// JetCleanser Package
// Questions/Comments? dkrohn@physics.harvard.edu mattlow@uchicago.edu schwartz@physics.harvard.edu liantaow@uchicago.edu
//
// Copyright (c) 2013
// David Krohn, Matthew Low, Matthew Schwartz, and Lian-Tao Wang
// 
// $Id$
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

#include "JetCleanser.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  // Modification to satisfy C++11 (thanks to Gavin Salam)
  const double JetCleanser::jc_zero = 1.0e-6;

  /////////////////////////////
  // constructor
  JetCleanser::JetCleanser(JetDefinition subjet_def, cleansing_mode cmode, input_mode imode) {
    _subjet_def = subjet_def;
    _rsub = subjet_def.R();
    _cleansing_mode = cmode;
    _input_mode = imode;
    _SetDefaults();
  }

  /////////////////////////////
  // constructor
  JetCleanser::JetCleanser(double rsub, cleansing_mode cmode, input_mode imode) {
    JetDefinition subjet_def(kt_algorithm, rsub);
    _subjet_def = subjet_def;
    _rsub = rsub;
    _cleansing_mode = cmode;
    _input_mode = imode;
    _SetDefaults();
  }

  /////////////////////////////
  // _SetDefaults
  void JetCleanser::_SetDefaults(){

    // defaults
    _fcut = 0.0;
    _nsjmin = -1;

    // "undefined" values
    _linear_gamma0_mean = -1;
    _gaussian_gamma0_mean = -1;
    _gaussian_gamma1_mean = -1;
    _gaussian_gamma0_width = -1;
    _gaussian_gamma1_width = -1;
  }
 
  /////////////////////////////
  // SetGroomingParameters
  void JetCleanser::SetGroomingParameters(double fcut, int nsjmin) {
    if ( fcut < 0 || fcut > 1 ) throw Error("SetGroomingParameters(): fcut must be >= 0 and <= 1");
    _fcut = fcut;
    _nsjmin = nsjmin;
  }

  /////////////////////////////
  // SetLinearParameters
  void JetCleanser::SetLinearParameters(double g0_mean) {
    if ( g0_mean < 0 || g0_mean > 1 ) throw Error("SetLinearParameters(): g0_mean must be >= 0 and <= 1");
    _linear_gamma0_mean = g0_mean;
  }

  /////////////////////////////
  // SetGaussianParameters
  void JetCleanser::SetGaussianParameters(double g0_mean, double g1_mean, double g0_width, double g1_width) {
    if ( g0_mean < 0 || g0_mean > 1 ) throw Error("SetGaussianParameters(): g0_mean must be >= 0 and <= 1");
    if ( g1_mean < 0 || g1_mean > 1 ) throw Error("SetGaussianParameters(): g1_mean must be >= 0 and <= 1");
    if ( g0_width < 0 || g0_width > 1 ) throw Error("SetGaussianParameters(): g0_width must be >= 0 and <= 1");
    if ( g1_width < 0 || g1_width > 1 ) throw Error("SetGaussianParameters(): g1_width must be >= 0 and <= 1");
    _gaussian_gamma0_mean = g0_mean;
    _gaussian_gamma1_mean = g1_mean;
    _gaussian_gamma0_width = g0_width;
    _gaussian_gamma1_width = g1_width;
  }

  /////////////////////////////
  // description
  std::string JetCleanser::description() const {
    std::ostringstream oss;
    oss << "JetCleanser [";
    if ( _cleansing_mode == jvf_cleansing ) oss << "JVF mode, ";
    else if ( _cleansing_mode == linear_cleansing ) oss << "Linear mode, ";
    else if ( _cleansing_mode == gaussian_cleansing ) oss << "Gaussian mode, ";
    if ( _input_mode == input_nc_together ) oss << "input = neutral and charged together]" << std::endl;
    else if ( _input_mode == input_nc_separate ) oss << "input = neutral and charged separate]" << std::endl;

    if ( _nsjmin <= 0 ) oss << " Trimming: fcut = " << _fcut << std::endl;
    else if ( _fcut >= 1.0 ) oss << " Filtering: nsj = " << _nsjmin << std::endl;
    else oss << " Trimming + Filtering: fcut = " << _fcut << ", nsj = " << _nsjmin << std::endl;

    if ( _cleansing_mode == linear_cleansing ) oss << " g0_mean = " << _linear_gamma0_mean << std::endl;
    else if ( _cleansing_mode == gaussian_cleansing ) oss << " g0_mean = " << _gaussian_gamma0_mean
                                                         << ", g0_width = " << _gaussian_gamma0_width
                                                         << ", g1_mean = " << _gaussian_gamma1_mean
                                                         << ", g1_width = " << _gaussian_gamma1_width << std::endl;
    return oss.str();
  }


  /////////////////////////////
  // cleanse
  PseudoJet JetCleanser::operator()(const PseudoJet & jet, 
                                    const std::vector<PseudoJet> & tracks_lv, 
                                    const std::vector<PseudoJet> & tracks_pu) const {

    if ( _input_mode != input_nc_together ) throw Error("result(): This operator is only defined for input_nc_together mode");

    // get constituents
    if ( !(jet.has_constituents()) ) return PseudoJet();
    std::vector<PseudoJet> constituents_all = jet.constituents();

    // prepare for clustering
    std::vector< std::vector<PseudoJet> > follow_sets;
    follow_sets.push_back( constituents_all );
    follow_sets.push_back( tracks_lv );
    follow_sets.push_back( tracks_pu );

    // get subjets
    std::vector< std::vector<PseudoJet> > sets = ClusterSets(_subjet_def, constituents_all, follow_sets);
    std::vector<PseudoJet> subjets_all, subjets_tracks_lv, subjets_tracks_pu;
    subjets_all       = sets[0];
    subjets_tracks_lv = sets[1];
    subjets_tracks_pu = sets[2];

    // rescale subjets
    std::vector<PseudoJet> rescaled_subjets_all;
    for (unsigned i=0; i<subjets_all.size(); i++){
      double s_factor = _GetSubjetRescaling_nctogether(subjets_all[i].pt(), subjets_tracks_lv[i].pt(), subjets_tracks_pu[i].pt());
      PseudoJet rescaled_subjet_all = RescalePseudoJetConstituents( subjets_all[i], s_factor );
      if ( rescaled_subjet_all != 0.0 ) rescaled_subjets_all.push_back( rescaled_subjet_all );
    }

    // trim/filter subjets
    rescaled_subjets_all = sorted_by_pt( rescaled_subjets_all );
    std::vector<PseudoJet> trimmed_subjets_all;
    for (unsigned i=0; i<rescaled_subjets_all.size(); i++){
      bool pass_filtering = _nsjmin > 0 ? i < _nsjmin  : false;
      bool pass_trimming  = rescaled_subjets_all[i].pt() > _fcut*jet.pt();
      if ( pass_trimming || pass_filtering ) { trimmed_subjets_all.push_back( rescaled_subjets_all[i] ); }
    }

    return join( trimmed_subjets_all );
  }

  /////////////////////////////
  // cleanse
  PseudoJet JetCleanser::operator()(const std::vector<PseudoJet> & neutrals_all,
                                    const std::vector<PseudoJet> & tracks_lv, 
                                    const std::vector<PseudoJet> & tracks_pu) const {

    if ( _input_mode != input_nc_separate ) throw Error("result(): This operator is only defined for input_nc_separate mode");

    // assemble jet collections
    std::vector<PseudoJet> particles_all;
    for (unsigned i=0; i<neutrals_all.size(); i++){ particles_all.push_back(neutrals_all[i]); }
    for (unsigned i=0; i<tracks_lv.size(); i++){ particles_all.push_back(tracks_lv[i]); }
    for (unsigned i=0; i<tracks_pu.size(); i++){ particles_all.push_back(tracks_pu[i]); }
    PseudoJet jet = join(particles_all);

    // prepare for clustering
    std::vector< std::vector<PseudoJet> > follow_sets;
    follow_sets.push_back( particles_all );
    follow_sets.push_back( neutrals_all );
    follow_sets.push_back( tracks_lv );
    follow_sets.push_back( tracks_pu );

    // get subjets
    std::vector< std::vector<PseudoJet> > sets = ClusterSets(_subjet_def, particles_all, follow_sets);
    std::vector<PseudoJet> subjets_all, subjets_neutrals_all, subjets_tracks_lv, subjets_tracks_pu;
    subjets_all          = sets[0];
    subjets_neutrals_all = sets[1];
    subjets_tracks_lv    = sets[2];
    subjets_tracks_pu    = sets[3];

    // rescale neutral subjets and add to charged LV subjets
    std::vector<PseudoJet> rescaled_subjets_all;
    for (unsigned i=0; i<subjets_all.size(); i++){
      double s_factor = _GetSubjetRescaling_ncseparate(subjets_neutrals_all[i].pt(), subjets_tracks_lv[i].pt(), subjets_tracks_pu[i].pt());
      PseudoJet rescaled_subjet_ntrl = RescalePseudoJetConstituents( subjets_neutrals_all[i], s_factor );
      PseudoJet rescaled_subjet_all = join( rescaled_subjet_ntrl, subjets_tracks_lv[i] );
      if ( rescaled_subjet_all != 0.0 ) rescaled_subjets_all.push_back( rescaled_subjet_all );
    }

    // trim/filter subjets
    rescaled_subjets_all = sorted_by_pt( rescaled_subjets_all );
    std::vector<PseudoJet> trimmed_subjets_all;
    for (unsigned i=0; i<rescaled_subjets_all.size(); i++){
      bool pass_filtering = _nsjmin > 0 ? i < _nsjmin  : false;
      bool pass_trimming  = rescaled_subjets_all[i].pt() > _fcut*jet.pt();
      if ( pass_trimming || pass_filtering ) { trimmed_subjets_all.push_back( rescaled_subjets_all[i] ); }
    }

    return join( trimmed_subjets_all );
  }


  /////////////////////////////
  // helper: _CheckRescalingValues
  //  allow some leeway for detector effects
  void JetCleanser::_CheckRescalingValues(double & pt_all, const double & ptc_lv, const double & ptc_pu) const {
    double ratio = (ptc_lv + ptc_pu)/pt_all;
    if ( ratio > 1.05 ) throw Error("_CheckRescalingValues: ptc_lv + ptc_pu is more than 5\% larger than pt_all");
    else if ( ratio > 1.0 ){ pt_all = pt_all * ratio; }
  }


  /////////////////////////////
  // helper: _GetSubjetRescaling_nctogether
  double JetCleanser::_GetSubjetRescaling_nctogether(double pt_all, double ptc_lv, double ptc_pu) const {
    double scale;

    // jvf mode
    if ( _cleansing_mode == jvf_cleansing ) {
      scale = ptc_lv > jc_zero ? ptc_lv / ( ptc_lv + ptc_pu ) : 0.0;
    }

    // linear mode
    else if ( _cleansing_mode == linear_cleansing ) {
      if ( _linear_gamma0_mean < 0 ) throw Error("Linear cleansing parameters have not been set yet.");
      _CheckRescalingValues(pt_all,ptc_lv,ptc_pu);

      if ( ptc_pu > jc_zero && ptc_pu / ( pt_all - ptc_lv ) > _linear_gamma0_mean )
        scale = ptc_lv > jc_zero ? ptc_lv / ( ptc_lv + ptc_pu ) : 0.0;
      else scale = ptc_lv > jc_zero ? 1.0 - (1.0/_linear_gamma0_mean)*ptc_pu/pt_all : 0.0;
      //double linear_gamma1 = ptc_lv >= jc_zero ? ptc_lv / ( pt_all - (1.0/_linear_gamma0_mean)*ptc_pu ) : 0.0;
    }

    // gaussian mode
    else if ( _cleansing_mode == gaussian_cleansing ) {
      if ( _gaussian_gamma0_mean < 0 || _gaussian_gamma1_mean < 0 || _gaussian_gamma0_width < 0 || _gaussian_gamma1_width < 0 )
        throw Error("Gaussian cleansing parameters have not been set yet.");
      _CheckRescalingValues(pt_all,ptc_lv,ptc_pu);

      double _gaussian_gamma0 = _GaussianGetMinimizedGamma0(pt_all, ptc_lv, ptc_pu);
      //double _gaussian_gamma1 = _GaussianGetGamma1(_gaussian_gamma0, pt_all, ptc_lv, ptc_pu);
      scale = ptc_lv > jc_zero ? 1.0 - (1.0/_gaussian_gamma0)*ptc_pu/pt_all : 0.0;
    }
    else throw Error("_GetSubjetRescaling: Current cleansing mode is not recognized.");

    return scale > jc_zero? scale : 0.0;
  }

  /////////////////////////////
  // helper: _GetSubjetRescaling_ncseparate
  double JetCleanser::_GetSubjetRescaling_ncseparate(double ptn_all, double ptc_lv, double ptc_pu) const {
    double scale;

    // jvf mode
    if ( _cleansing_mode == jvf_cleansing ) {
      scale = ptc_lv > jc_zero && ptn_all > jc_zero ? ptc_lv / ( ptc_lv + ptc_pu ) : 0.0;
    }

    // linear mode
    else if ( _cleansing_mode == linear_cleansing ) {
      if ( _linear_gamma0_mean < 0 ) throw Error("Linear cleansing parameters have not been set yet.");
      double pt_all = ptn_all + ptc_lv + ptc_pu;
      _CheckRescalingValues(pt_all,ptc_lv,ptc_pu);

      if ( (ptc_pu > jc_zero && ptc_pu / ( pt_all - ptc_lv ) > _linear_gamma0_mean) || ptn_all < jc_zero ) 
        scale = ptc_lv > jc_zero && ptn_all > jc_zero? ptc_lv / ( ptc_lv + ptc_pu ) : 0.0;
      else scale = ptc_lv > jc_zero && ptn_all > jc_zero? 1.0 - (1.0/_linear_gamma0_mean-1.0)*ptc_pu/ptn_all : 0.0;
      //double linear_gamma1 = ptc_lv > jc_zero ? ptc_lv / ( pt_all - (1.0/_linear_gamma0_mean)*ptc_pu ) : 0.0;
    }

    // gaussian mode
    else if ( _cleansing_mode == gaussian_cleansing ) {
      if ( _gaussian_gamma0_mean < 0 || _gaussian_gamma1_mean < 0 || _gaussian_gamma0_width < 0 || _gaussian_gamma1_width < 0 )
        throw Error("Gaussian cleansing parameters have not been set yet.");
      double pt_all = ptn_all + ptc_lv + ptc_pu;
      _CheckRescalingValues(pt_all,ptc_lv,ptc_pu);

      double _gaussian_gamma0 = _GaussianGetMinimizedGamma0(pt_all, ptc_lv, ptc_pu);
      //double _gaussian_gamma1 = _GaussianGetGamma1(_gaussian_gamma0, pt_all, ptc_lv, ptc_pu);
      scale = ptc_lv > jc_zero && ptn_all > jc_zero ? 1.0 - (1.0/_gaussian_gamma0-1.0)*ptc_pu/ptn_all : 0.0;
    }
    else throw Error("_GetSubjetRescaling: Current cleansing mode is not recognized.");

    return scale > jc_zero? scale : 0.0;
  }


  /////////////////////////////
  // helper: _GaussianGetMinimizedGamma0
  double JetCleanser::_GaussianGetMinimizedGamma0(double pt_all, double ptc_lv, double ptc_pu) const {
    if ( pt_all == 0.0 && ptc_lv == 0.0 && ptc_pu == 0.0 ) return 0.0;
    if ( ptc_lv == 0.0 ) return ptc_pu / pt_all;

    double params[3] = {ptc_lv, ptc_pu, pt_all};
    std::map<double,double> map_fcn_to_x;

    for(double x0 = 0.0; x0 <= 1.0 + jc_zero; x0+=0.01) {
      map_fcn_to_x[_GaussianFunction(x0,params)] = x0;  
    }

    return map_fcn_to_x.begin()->second;
  }

  /////////////////////////////
  // helper: _GaussianGetGamma1
  double JetCleanser::_GaussianGetGamma1(double gamma0, double pt_all, double ptc_lv, double ptc_pu) const {
    if ( pt_all == 0.0 && ptc_lv == 0.0 && ptc_pu == 0.0 ) return 0.0;
    if ( gamma0 == 0.0 || fabs(pt_all - ptc_pu/gamma0) < jc_zero ) return 0.0;
    return ptc_lv / (pt_all - ptc_pu/gamma0);
  }

  /////////////////////////////
  // helper: _GaussianFunction
  double JetCleanser::_GaussianFunction(double x, void * params) const {
    double ptc_lv = ((double*)params)[0];
    double ptc_pu = ((double*)params)[1];
    double pt_all = ((double*)params)[2];
    double g1 = _GaussianGetGamma1(x, pt_all, ptc_lv, ptc_pu);

    if(g1 >= 1. || g1 <= 0. || x <= 0. || x >= 1.)
      return (x-1.)*(x-1.)+10.;

    return -exp(-(g1-_gaussian_gamma1_mean)*(g1-_gaussian_gamma1_mean)/2./_gaussian_gamma1_width/_gaussian_gamma1_width
           -(x-_gaussian_gamma0_mean)*(x-_gaussian_gamma0_mean)/2./_gaussian_gamma0_width/_gaussian_gamma0_width);
  }
  /////////////////////////////
  // helper: ClusterSets
  std::vector< std::vector<PseudoJet> > ClusterSets(const JetDefinition & jet_def, 
                                                    const std::vector<PseudoJet> & cluster_set,
                                                    const std::vector< std::vector<PseudoJet> > & follow_sets,
                                                    const double &ptmin) {
    // start set
    std::vector<PseudoJet> full_set;
    for (unsigned i=0; i<cluster_set.size(); i++){ full_set.push_back( cluster_set[i] ); }
   
    // convert to ghosts and add
    for (unsigned i=0; i<follow_sets.size(); i++){
      std::vector<PseudoJet> current_set = follow_sets[i];
      for (unsigned j=0; j<current_set.size(); j++){
        FollowSetGhostInfo* ghost_info = new FollowSetGhostInfo(i,j);
        PseudoJet ghost_pseudojet = 1.0e-60 * current_set[j];
        ghost_pseudojet.set_user_info( ghost_info );
        full_set.push_back( ghost_pseudojet );
      }
    }

    // cluster
    ClusterSequence *cs = new ClusterSequence( full_set, jet_def );
    std::vector<PseudoJet> jets = sorted_by_pt( cs->inclusive_jets(ptmin) );
    cs->delete_self_when_unused();

    // construct sets
    std::vector< std::vector<PseudoJet> > follow_jets;
    for (unsigned i=0; i<follow_sets.size(); i++){
      std::vector<PseudoJet> current_set;
      for (unsigned j=0; j<jets.size(); j++){ current_set.push_back( join( PseudoJet() ) ); }
      follow_jets.push_back( current_set );
    }
   
    // construct jet sets
    for (unsigned i=0; i<jets.size(); i++){
      std::vector<PseudoJet> constituents = jets[i].constituents();
      for (unsigned j=0; j<constituents.size(); j++){
        if ( constituents[j].has_user_info<FollowSetGhostInfo>() ){
          FollowSetGhostInfo ghost_info = constituents[j].user_info<FollowSetGhostInfo>();
          int set_id = ghost_info.GetSetId();
          int ind_id = ghost_info.GetIndId();

          std::vector<PseudoJet> current_set = follow_sets[set_id];
          if ( follow_jets[set_id][i] == 0.0 ) follow_jets[set_id][i] = join( current_set[ind_id] );
          else follow_jets[set_id][i] = join( follow_jets[set_id][i], current_set[ind_id] );
        }
      }
    }

    return follow_jets;
  }

  /////////////////////////////
  // helper: RescalePseudoJetVector
  std::vector<PseudoJet> RescalePseudoJetVector(const std::vector<PseudoJet> & jets, const double s_factor) {
    std::vector<PseudoJet> rescaled_jets;
    if ( s_factor == 0.0 ) return rescaled_jets;
    for (unsigned i=0; i<jets.size(); i++){ rescaled_jets.push_back( s_factor*jets[i] ); }
    //for (unsigned i=0; i<jets.size(); i++){ rescaled_jets.push_back( PseudoJet(s_factor*jets[i].px(),
    //                                                                           s_factor*jets[i].py(),
    //                                                                           s_factor*jets[i].pz(),
    //                                                                           s_factor*jets[i].e()) ); }
    return rescaled_jets;
  }

  /////////////////////////////
  // helper: RescalePseudoJetConstituents
  PseudoJet RescalePseudoJetConstituents(const PseudoJet & jet, const double s_factor) {
    if ( !(jet.has_constituents()) )
      return PseudoJet();
    std::vector<PseudoJet> constituents = jet.constituents();
    std::vector<PseudoJet> rconstituents = RescalePseudoJetVector( constituents, s_factor );
    return join( rconstituents );
  }

  /////////////////////////////
  // helper: _RunTests
  void JetCleanser::_RunTests() {
    std::cout << "----- Testing contrib::JetCleanser -----" << std::endl;
    std::cout << "Warning: All cleansing settings will be changed during the test." << std::endl;

    // INPUT MODE = input_nc_together
    _input_mode     = input_nc_together;
    _cleansing_mode = jvf_cleansing;
    std::cout << "Mode = [nc_together,jvf]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

    _cleansing_mode = linear_cleansing;
    _linear_gamma0_mean = 0.67;
    std::cout << std::endl << "Mode = [nc_together,linear]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

    _cleansing_mode = gaussian_cleansing;
    _gaussian_gamma0_mean  = 0.67;
    _gaussian_gamma0_width = 0.15;
    _gaussian_gamma1_mean  = 0.67;
    _gaussian_gamma1_width = 0.25;
    std::cout << std::endl << "Mode = [nc_together,gaussian]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

    // INPUT MODE = input_nc_together
    _input_mode     = input_nc_separate;
    _cleansing_mode = jvf_cleansing;
    std::cout << std::endl << "Mode = [nc_separate,jvf]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

    _input_mode     = input_nc_separate;
    _cleansing_mode = linear_cleansing;
    _linear_gamma0_mean = 0.67;
    std::cout << std::endl << "Mode = [nc_separate,linear]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

    _cleansing_mode = gaussian_cleansing;
    _gaussian_gamma0_mean  = 0.67;
    _gaussian_gamma0_width = 0.15;
    _gaussian_gamma1_mean  = 0.67;
    _gaussian_gamma1_width = 0.25;
    std::cout << std::endl << "Mode = [nc_separate,gaussian]" << std::endl;
    _RunTestRescaling(1.0 , 0.0 ,-0.5 );
    _RunTestRescaling(1.0 , 0.0 ,-0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.0 );
    _RunTestRescaling(1.0 , 0.0 , 0.01);
    _RunTestRescaling(1.0 , 0.0 , 0.50);
    _RunTestRescaling(1.0 , 0.0 , 0.99);
    _RunTestRescaling(1.0 , 0.0 , 1.0 );
    _RunTestRescaling(1.0 , 0.0 , 1.01);
    _RunTestRescaling(1.0 , 0.0 , 1.5 );

    _RunTestRescaling(1.0 , 0.01, 1.0 );
    _RunTestRescaling(1.0 , 0.5 , 0.3 );
    _RunTestRescaling(1.0 , 0.5 , 0.49);
    _RunTestRescaling(1.0 , 0.5 , 0.5 );
    _RunTestRescaling(1.0 , 0.5 , 0.51);
    _RunTestRescaling(1.0 , 0.5 , 0.7 );
    _RunTestRescaling(1.0 , 0.5 , 1.0 );

    _RunTestRescaling(1.0 ,-0.5 , 0.0 );
    _RunTestRescaling(1.0 ,-0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.01, 0.0 );
    _RunTestRescaling(1.0 , 0.50, 0.0 );
    _RunTestRescaling(1.0 , 0.99, 0.0 );
    _RunTestRescaling(1.0 , 1.0 , 0.0 );
    _RunTestRescaling(1.0 , 1.01, 0.0 );
    _RunTestRescaling(1.0 , 1.5 , 0.0 );

  }

  /////////////////////////////
  // helper: _RunTestRescaling
  void JetCleanser::_RunTestRescaling(double pt_all, double ptc_lv, double ptc_pu) const {
    double s_factor, ptn_all=0;
    
    if ( _input_mode == input_nc_separate ) ptn_all = pt_all - ptc_lv - ptc_pu;

    try { s_factor = _input_mode == input_nc_together? 
                     _GetSubjetRescaling_nctogether(pt_all, ptc_lv, ptc_pu) :
                     _GetSubjetRescaling_ncseparate(ptn_all, ptc_lv, ptc_pu); }
    catch(Error e) { s_factor = -1.0; }

    std::cout << " pt_all = " << pt_all
              << "   ptc_lv = " << ptc_lv
              << "   ptc_pu = " << ptc_pu;

    if ( _input_mode == input_nc_separate ) std::cout << "   ptn_all = " << ptn_all;
    if (s_factor<0.0) { std::cout << "   scale = error" << std::endl; }
    else                std::cout << "   scale = "  << s_factor << std::endl;

  }

} // namespace contrib
 
FASTJET_END_NAMESPACE
