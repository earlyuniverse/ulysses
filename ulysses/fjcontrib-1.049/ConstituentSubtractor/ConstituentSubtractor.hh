// $Id: ConstituentSubtractor.hh 1240 2020-02-23 13:51:05Z peter.berta $
//
// ConstituentSubtractor package
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu
//
// Copyright (c) 2014-, Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner
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

#ifndef __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__
#define __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__


#include <fastjet/internal/base.hh>
#include <fastjet/ClusterSequenceAreaBase.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/GridMedianBackgroundEstimator.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include "fastjet/tools/Transformer.hh" // to derive Subtractor from Transformer
#include "fastjet/LimitedWarning.hh"

#include <algorithm>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{


//------------------------------------------------------------------------
/// \class ConstituentSubtractor
/// A class to perform subtraction of background, e.g. pileup, from a set of particles at particle-level. The output is a jet or the whole event with corrected constituents.
///
/// This class corrects the input particles for background contamination with the algorithm described in:
/// Peter Berta, Martin Spousta, David W. Miller, Rupert Leitner [arXiv:1403.3108]
///
/// For individual jet background subtraction, see example_jet_by_jet.cc
/// For whole event background subtraction before jet clustering, see example_event_wide.cc
///
/// The distance used for matching between particle i and ghost k is defined as:
/// deltaR_{i,k}=pT_i * sin(theta_i)^polarAngleExp * sqrt((y_i-y_k)^2 + (phi_i-phi_k)^2)
///
/// The class accounts for position-dependent (in rapidity-azimuth plane) background densities, rho and rho_m. The user is encouraged to use them to get the maximal performance.
///
///
/// For the treatment of massive particles, choose one of the following options (currently the default option is 3 which seems to be one of the most optimal in our studies):      
/// 1. Keep the original mass and rapidity. Not recommended since jet mass is wrongly reconstructed. One can use it by adding before initialization:
///    subtractor.set_keep_original_masses();                
/// 2. Keep the original mass and pseudo-rapidity. Not recommended since jet mass is wrongly reconstructed. Use these two functions for this option:
///    subtractor.set_keep_original_masses();
///    subtractor.set_fix_pseudorapidity();
/// 3. Set all masses to zero. Keep rapidity unchanged. This is recommended and is the default option, since observed the best performance for high pt large-R jets. 
///    Nothing needs to be specified.
/// 4. Set all masses to zero. Keep pseudo-rapidity unchanged. Also recommended, almost the same performance as for option 3. Use function:
///    subtractor.set_fix_pseudorapidity();
/// 5. Do correction of m_delta=sqrt(p_T^2+m^2)-p_t. Not recommended. One can use it by adding before initialization:
///    subtractor.set_do_mass_subtraction();
///    One must additionally provide option for the rho_m background estimation. There are several possibilities:
///    a) Same background estimator as for rho. Use this function:  
///     subtractor.set_common_bge_for_rho_and_rhom();  // this must be used after set_background_estimator function. 
///    b) Use a separate background estimator - you need to specify it as an additional parameter in function set_background_estimator
///    c) Use scalar background density using function set_scalar_background_density(double rho, double rhom).
/// 6. Same as 5 just keep pseudo-rapidity unchanged. Not recommended. Additionally to what is written for option 5, use:                        
///    subtractor.set_fix_pseudorapidity();                 
/// 7. Keep rapidity and pseudo-rapidity fixed (scale fourmomentum). Recommended - observed better performance than the mass correction. Use:    
///    subtractor.set_scale_fourmomentum();



  namespace ConstituentSubtractorConstants{
    const double zero_pt=1e-50;  // This pt value is considered as zero. It is used for corrected particles with zero pt.
    const double zero_mass=1e-50; // This mass value is considered as zero. It is used for corrected particles with zero mass.
    // By default particles with zero mass and pt are discarded in case their corrected pt and mass are zero. The user can use the function "set_remove_zero_pt_and_mass_particles(false)" to keep such partiles.
  }


  class ConstituentSubtractor : public fastjet::Transformer{
  public:

    enum Distance {
      deltaR,     /// deltaR=sqrt((y_i-y_j)^2+(phi_i-phi_j)^2)), longitudinal Lorentz invariant
      angle   ///  angle between two momenta in Euclidean space                                          
    };
   
    ///
    /// default ctor
    ConstituentSubtractor();

    ///
    /// Constructor that takes a pointer to a background estimator for rho and optionally a pointer to a background estimator for rho_m.  If the latter is not supplied, rho_m is assumed to always be zero (this behaviour can be changed by calling use_common_bge_for_rho_and_rhom).
    ConstituentSubtractor(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom=0, double alpha=0, double max_distance=-1, Distance distance=deltaR);
    
    ///
    /// Constructor that takes an externally supplied value for rho and rho_m.
    ConstituentSubtractor(double rho, double rhom=0, double alpha=0, double max_distance=-1, Distance distance=deltaR);
    
    ///
    /// default dtor
    virtual ~ConstituentSubtractor(){}
    
    ///
    /// initialization. use it before event loop.
    virtual void initialize();
    
    ///
    /// Common description for jet-by-jet, event-wide, and iterative CS
    void description_common(std::ostringstream &descr) const;

    ///
    /// a description of what this class does
    virtual std::string description() const;
    
    ///
    /// action of the correction on a given jet. The output is PseudoJet object with subtracted constituents
    virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const;
    
    ///
    /// do the constituent subtraction for the input particles using the provided background proxies. The output is a vector with corrected particles - particles with zero corrected pt and mass are removed by default. The user can modify this behaviour by using functions set_remove_particles_with_zero_pt_and_mass.
    std::vector<fastjet::PseudoJet> do_subtraction(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const &backgroundProxies,std::vector<fastjet::PseudoJet> *remaining_backgroundProxies=0) const;
    

    ///
    /// this version should be not used. Use the version below.
    virtual std::vector<fastjet::PseudoJet> subtract_event(std::vector<fastjet::PseudoJet> const &particles, double max_eta);

    ///
    /// do the subtraction of the whole event - more user-friendly approach. The particles with |eta|>max_eta are discarded at the beginning, i.e. they are not used, nor returned. The ghosts are added automatically inside this function up to max_eta.
    virtual std::vector<fastjet::PseudoJet> subtract_event(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const *hard_proxies=0);
    
    ///
    /// do the subtraction of the whole event using the tracking information for charged particles, i.e. the 4-momenta of charged particles from signal vertex, and 4-momenta of charged particles from background. The user can set the scaling of charged particles from background and signal using parameters charged_background_scale (CBS) and charged_signal_scale (CSS). These scales are useful if one assumes correlation between charged and neutral particles or in case the inputs from calorimeter are miscalibrated wrt tracks. In case CBS=CSS=0, the input charged particles are not used. In case CBS=CSS=1, the input charged particles are not scaled. Recommending to try several combinations for CBS and CSS from range [0.8, 1.5]. It is no more necessary to provide background estimator. The GridMedianBackgroundEstimator is used - probably more flexibility will be added in the future. The rescaling function for background estimator can be also provided - the rescaling function will be used for the event after subtracting charged scaled particles.  The particles with |eta|>max_eta are discarded at the beginning, i.e. they are not used, nor returned.
    std::vector<fastjet::PseudoJet> subtract_event_using_charged_info(std::vector<fastjet::PseudoJet> const &particles, double charged_background_scale, std::vector<fastjet::PseudoJet> const &charged_background, double charged_signal_scale, std::vector<fastjet::PseudoJet> const &charged_signal, double max_eta);

    
    void set_rescaling(fastjet::FunctionOfPseudoJet<double> *rescaling);

    ///
    /// set the grid size (not area) for the background estimation with GridMedianBackgroundEstimator  when using charged info
    void set_grid_size_background_estimator(double const &grid_size_background_estimator);
    

    ///
    /// Set the pointer to a background estimator for rho and optionally a pointer to a background estimator for rho_m.  If the latter is not supplied, rho_m is assumed to always be zero (this behaviour can be changed by calling use_common_bge_for_rho_and_rhom).
    void set_background_estimator(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom=0);
    
    ///
    /// Set the scalar background densities rho and rho_m.
    void set_scalar_background_density(double rho, double rhom=0);
    
    ///
    /// When only one background estimator, bge_rho, is specified, calling this function with argument true, causes rho_m to be calculated from the same background estimator as rho, instead of being set to zero. Currently this only works if the estimator is a JetMedianBackgroundEstimator or other estimator which has such function.
    void set_common_bge_for_rho_and_rhom();

    ///
    /// This function should not be used anymore (instead the function above should be used without any parameter). When used with "true", then it has the same behaviour as set_common_bge_for_rho_and_rhom() above. By default, no mass subtraction is done, so it makes no sense to use this function with "false".
    void set_common_bge_for_rho_and_rhom(bool value);


    ///
    /// By using this function, the original masses of particles are kept. By default, the masses of all particles are set to zero. If you do not want to set the masses to zero but you want to do the subtraction for the mass part, use the function set_common_bge_for_rho_and_rhom or specify background estimator for rho_m.
    void set_keep_original_masses();

    ///
    /// When two background estimators are used (one for rho, the second for rho_m), setting this to true will result in rho_m being estimated using bge_rhom->rho_m() instead of bge_rhom->rho().
    ///  void set_use_bge_rhom_rhom(bool value=true);
    
    ///
    /// With this function, the user specifies if also the mass term sqrt(pT^2+m^2)-pT should be corrected during the subtraction procedure. This is automatically set when calling functions set_common_bge_for_rho_and_rhom, set_scalar_background_density or set_use_bge_rhom_rhom
    void set_do_mass_subtraction();
    
    ///
    /// set to true, if you want to remove particles which have zero both, pt and mass. By default, these particles are removed. 
    void set_remove_particles_with_zero_pt_and_mass(bool value=true);
    
    ///
    /// set to true, if you want to remove all zero pt particles - this means removing also particles with zero delta_m (massive particles with zero pt). By default, the zero pt particles with non-zero delta_m are not removed. 
    void set_remove_all_zero_pt_particles(bool value=true);
    
    ///
    /// function to change the alpha-parameter figuring in the distance measure deltaR. The larger the alpha, the more are preferred to be corrected the low pt particles. The default value is 0, i.e. by default the standard deltaR definition is used: deltaR=sqrt(deltay^2 + deltaphi^2)
    void set_alpha(double alpha);
    
    ///
    /// function to change the parameter polarAngleExp
    void set_polarAngleExp(double polarAngleExp);
    

    ///
    /// function to change the parameter ghost_area
    void set_ghost_area(double ghost_area);

    ///
    /// function to change the distance type
    void set_distance_type(Distance distance);

    ///
    /// function to change the free parameter max_distance. The distance is defined by enum Distance. The particle-ghost pairs with distance>max_distance are not used. When max_distance<=0, the max_distance parameter is not used (no upper limit on distance). The default value is -1, i.e. by default there is no upper limit for possible distance values. 
    void set_max_distance(double max_distance);
     
    ///
    /// function to change the free parameter max_distance. It is the same as the set_max_distance function. It is added back for backward-comptibility (since it was replaced in ConstituentSubtractor version 1.1.3 with function set_max_distance). It should be not used.
    void set_max_standardDeltaR(double max_distance);
    
    ///    
    /// function to return the maximal deltaR distance used in the subtraction
    double get_max_distance();
    
    ///
    /// This function returns true if the first argument is smaller than the second argument, otherwise returns false. The comparison is done only on the first element in the two pairs. This function is used to sort in ascending order the deltaR values for each pair particle-ghost while keeping track of particles and ghosts
    static bool _function_used_for_sorting(std::pair<double,std::pair<int,int> >  const &i,std::pair<double,std::pair<int,int> >  const &j);

    ///
    /// function used for sorting of Pseudojets according to eta
    static bool _rapidity_sorting(fastjet::PseudoJet const &i,fastjet::PseudoJet const &j);

    ///
    /// Function to construct massless ghosts in the whole eta-phi space up to max_eta
    void construct_ghosts_uniformly(double max_eta);

    ///
    /// Function to return vector of ghosts
    std::vector<fastjet::PseudoJet>  get_ghosts();

    ///
    /// Function to set the maximal eta value for subraction
    void set_max_eta(double max_eta);

    ///
    /// Function to return vector of areas associated to ghosts
    std::vector<double>  get_ghosts_area();


    ///
    /// Set ghost selector - this allows to perform Constituent Subtraction in certain phase space when using the event-wide correction. Only ghosts which pass the selector are used in the correction. No selection is done on particles (i.e. the user needs to perform the selection before applying ConstituentSubtractor if he/she wishes). The default behaviour is to not use any selector. The user should remember that the max_eta restriction is applied in any case, e.g. when using the function "subtract_event(particles)", then ghosts are distributed to maximal eta max_eta and particles outside this region are discarded.
    void set_ghost_selector(fastjet::Selector* selector);

    ///
    /// Set particle selector - if specified, only particles which pass this selector are corrected. The other particles are copied to the final corrected vector without modification
    void set_particle_selector(fastjet::Selector* selector);

    ///
    /// Use this function if you want to use the information about hard particle proxies in the subtraction procedure (typically the charged tracks from primary vertex and/or high pt calorimeter clusters). The parameters of this function are the "nearby_hard_radius" and "nearby_hard_factor". If there is at least one hard particle proxy within Distance specified by the nearby_hard_radius parameter, then the CS distance is multiplied by nearby_hard_factor.
    void set_use_nearby_hard(double const &nearby_hard_radius, double const &nearby_hard_factor);

    ///
    /// Use this function, if you want to fix the pseudo-rapidity instead of rapidity for the corrected particles. By default, the rapidity is fixed.
    void set_fix_pseudorapidity();

    ///
    /// Use this function if you want to fix both, pseudo-rapidity and rapidity. By default, rapidity is fixed and the mass is set to zero.
    void set_scale_fourmomentum();


  protected:
    std::vector<fastjet::PseudoJet> get_background_proxies_from_ghosts(std::vector<fastjet::PseudoJet> const &ghosts,std::vector<double> const &ghosts_area) const;
    void clear_ghosts();
    
    ///
    /// common function to initialize for Iterative and standard CS
    void _initialize_common();

    ///
    /// helper functions to find the minimal and maximal indeces for the vector of particles - useful to make the CS procedure faster
    unsigned int _find_index_after(double const &value, std::vector<double> const &vec) const;
    unsigned int _find_index_before(double const &value, std::vector<double> const &vec) const;

    ///
    /// function to get the trabsformed distance
    double _get_transformed_distance(double const &distance) const;

    fastjet::BackgroundEstimatorBase *_bge_rho, *_bge_rhom;
    bool _common_bge, _rhom_from_bge_rhom;
    double _rho, _rhom;
    bool _externally_supplied_rho_rhom, _do_mass_subtraction;
    Distance _distance;
    double _alpha;
    double _polarAngleExp;
    double _max_distance;
    bool _remove_particles_with_zero_pt_and_mass;
    bool _remove_all_zero_pt_particles;
    bool _use_max_distance;
    double _ghost_area;
    double _grid_size_phi;
    double _grid_size_rap;
    bool _ghosts_constructed;
    bool _ghosts_rapidity_sorted;
    int _n_ghosts_phi;
    double _max_eta;
    bool _masses_to_zero;
    bool _use_nearby_hard;
    double _nearby_hard_radius;
    double _nearby_hard_factor;
    bool _fix_pseudorapidity;
    bool _scale_fourmomentum;
    std::vector<fastjet::PseudoJet> const *_hard_proxies;

    std::vector<fastjet::PseudoJet> _ghosts;
    std::vector<double> _ghosts_area;
    std::vector<double> _ghosts_rapidities; // used for uniform grid
    double _grid_size_background_estimator;
    fastjet::Selector *_ghost_selector;
    fastjet::Selector *_particle_selector;
    fastjet::FunctionOfPseudoJet<double> *_rescaling;
    static LimitedWarning _warning_unused_rhom;
  };
  


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_CONSTITUENTSUBTRACTOR_HH__
