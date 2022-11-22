// $Id: ConstituentSubtractor.cc 1240 2020-02-23 13:51:05Z peter.berta $
//
// ConstituentSubtractor package
// Questions/comments: berta@ipnp.troja.mff.cuni.cz, Martin.Spousta@cern.ch, David.W.Miller@uchicago.edu, Rupert.Leitner@mff.cuni.cz
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

#include "ConstituentSubtractor.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

LimitedWarning ConstituentSubtractor::_warning_unused_rhom;


  ///
  /// default constructor
  ConstituentSubtractor::ConstituentSubtractor(){
    _bge_rho=0;
    _bge_rhom=0;
    _common_bge=false;
    _rhom_from_bge_rhom=false;
    _externally_supplied_rho_rhom=false;
    _do_mass_subtraction=false;
    _masses_to_zero=true;
    _distance=deltaR;
    _alpha=0;
    _polarAngleExp=0;
    _max_distance=-1;
    _remove_particles_with_zero_pt_and_mass=true;
    _remove_all_zero_pt_particles=false;
    _use_max_distance=false;
    _ghost_area=0.01;
    _grid_size_phi=-1;
    _grid_size_rap=-1;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _max_eta=-1;
    _grid_size_background_estimator=0.5;
    _rescaling=0;
    _use_nearby_hard=false;
    _fix_pseudorapidity=false;
    _scale_fourmomentum=false;
    _hard_proxies=0;
    _ghost_selector=0;
    _particle_selector=0;
  }



  /// Constructor that takes a pointer to a background estimator for
  /// rho and optionally a pointer to a background estimator for
  /// rho_m.
  ConstituentSubtractor::ConstituentSubtractor(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom, double alpha, double max_distance, Distance distance) :
    _bge_rho(bge_rho), _bge_rhom(bge_rhom), _common_bge(false),  _externally_supplied_rho_rhom(false), _distance(distance), _alpha(alpha), _max_distance(max_distance){
    if (_max_distance>0) _use_max_distance=true;
    else _use_max_distance=false;
    _polarAngleExp=0;
    _ghost_area=0.01;
    _remove_particles_with_zero_pt_and_mass=true;
    _remove_all_zero_pt_particles=false;
    _max_eta=-1;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _do_mass_subtraction=false;
    _masses_to_zero=true;
    _use_nearby_hard=false;
    _fix_pseudorapidity=false;
    _scale_fourmomentum=false;
    _hard_proxies=0;
    _ghost_selector=0;
    _particle_selector=0;
  }


  ///
  /// Constructor that takes an externally supplied value for rho and, optionally, for rho_m.
  ConstituentSubtractor::ConstituentSubtractor(double rho, double rhom, double alpha, double max_distance, Distance distance) :      
    _bge_rho(0), _bge_rhom(0), _common_bge(false), _rhom_from_bge_rhom(false), _rho(rho), _rhom(rhom), _externally_supplied_rho_rhom(true), _distance(distance), _alpha(alpha), _max_distance(max_distance) {
    if (_max_distance>0) _use_max_distance=true;
    else _use_max_distance=false;
    assert(_rho  >= 0);
    assert(_rhom >= 0);

    _do_mass_subtraction=false;
    _masses_to_zero=true;
    _polarAngleExp=0;
    _ghost_area=0.01;
    _remove_particles_with_zero_pt_and_mass=true;
    _remove_all_zero_pt_particles=false;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _max_eta=-1;
    _use_nearby_hard=false;
    _fix_pseudorapidity=false;
    _scale_fourmomentum=false;
    _hard_proxies=0;
    _ghost_selector=0;
    _particle_selector=0;
  }


  void ConstituentSubtractor::_initialize_common(){
    if (_max_eta<=0) throw Error("ConstituentSubtractor::initialize_common: The value for eta cut was not set or it is negative. It needs to be set before using the function initialize");
    if (_masses_to_zero && _do_mass_subtraction) throw Error("ConstituentSubtractor::initialize_common: It is specified to do mass subtraction and also to keep the masses at zero. Something is wrong.");
    if (_masses_to_zero && _scale_fourmomentum) throw Error("ConstituentSubtractor::initialize_common: It is specified to do scaling of fourmomenta and also to keep the masses at zero. Something is wrong.");
    if (_do_mass_subtraction && _scale_fourmomentum) throw Error("ConstituentSubtractor::initialize_common: It is specified to do mass subtraction and also to do scaling of fourmomenta. Something is wrong.");
    this->construct_ghosts_uniformly(_max_eta);
  }


  void ConstituentSubtractor::initialize(){
    this->_initialize_common();
  }



  void ConstituentSubtractor::set_background_estimator(fastjet::BackgroundEstimatorBase *bge_rho, fastjet::BackgroundEstimatorBase *bge_rhom){
    _bge_rho=bge_rho;
    _bge_rhom=bge_rhom;
  }


  void ConstituentSubtractor::set_scalar_background_density(double rho, double rhom){
    _rho=rho;
    _rhom=rhom;
    assert(_rho  >= 0);
    assert(_rhom >= 0);
    _externally_supplied_rho_rhom=true;
    _common_bge=false;
  }

  ///----------------------------------------------------------------------
  /// the action on a given jet
  fastjet::PseudoJet ConstituentSubtractor::result(const PseudoJet &jet) const{
    // make sure we have a BGE or a rho value
    if (!_bge_rho && !_externally_supplied_rho_rhom){
      throw Error("ConstituentSubtractor::result() constituent subtraction needs a BackgroundEstimator or a value for rho.");
    }
    if (_ghosts_constructed) throw Error("ConstituentSubtractor::result() The ghosts are constructed, but they are not needed when using this function. When you want to perform jet-by-jet correction, initialize a new ConstituentSubtractor without construction of ghosts.");

    ///----------------------------------------------------------------------
    /// sift ghosts and particles in the input jet
    std::vector<fastjet::PseudoJet> particles, ghosts;
    SelectorIsPureGhost().sift(jet.constituents(), ghosts, particles);

    std::vector<fastjet::PseudoJet> selected_particles,unselected_particles;
    if (_particle_selector) _particle_selector->sift(particles, selected_particles, unselected_particles);
    else selected_particles=particles;


    std::vector<double> ghosts_area;
    unsigned long nGhosts=ghosts.size();
    for (unsigned int j=0;j<nGhosts; ++j){
      ghosts_area.push_back(ghosts[j].area());
    }

    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(ghosts,ghosts_area);
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(particles,backgroundProxies);
    if (_particle_selector) subtracted_particles.insert(subtracted_particles.end(), unselected_particles.begin(), unselected_particles.end());
    fastjet::PseudoJet subtracted_jet=join(subtracted_particles);
    subtracted_jet.set_user_index(jet.user_index());

    return subtracted_jet;    
  }




  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::subtract_event(std::vector<fastjet::PseudoJet> const &particles, double max_eta){
    if (fabs(_max_eta/max_eta-1)>1e-5 && max_eta>0){
      _ghosts_constructed=false;
      _max_eta=max_eta;
    }
    if (!_ghosts_constructed) this->construct_ghosts_uniformly(_max_eta);
    return subtract_event(particles);
  }


  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::subtract_event(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const *hard_proxies){
    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    std::vector<fastjet::PseudoJet> selected_particles,unselected_particles;
    for (unsigned int iparticle=0;iparticle<particles.size();++iparticle){
      if (std::abs(particles[iparticle].eta())>_max_eta) continue;
      if (particles[iparticle].pt()<ConstituentSubtractorConstants::zero_pt && _remove_all_zero_pt_particles) continue;
      if (particles[iparticle].pt()<ConstituentSubtractorConstants::zero_pt && (_masses_to_zero || particles[iparticle].m()<ConstituentSubtractorConstants::zero_mass) && _remove_particles_with_zero_pt_and_mass) continue;
      if (_particle_selector){
	if (_particle_selector->pass(particles[iparticle])) selected_particles.push_back(particles[iparticle]);
	else unselected_particles.push_back(particles[iparticle]);
      }
      else selected_particles.push_back(particles[iparticle]);
    }
    if (_use_nearby_hard){
      if (hard_proxies) _hard_proxies=hard_proxies;
      else throw Error("ConstituentSubtractor::subtract_event: It was requested to use closeby hard proxies but they were not provided in this function!");
    }
    else if (hard_proxies) throw Error("ConstituentSubtractor::subtract_event: Hard proxies were provided but the set_use_hard_proxies function was not used before initialization. It needs to be called before initialization!");
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(selected_particles,backgroundProxies);
    if (_particle_selector) subtracted_particles.insert(subtracted_particles.end(), unselected_particles.begin(), unselected_particles.end());
    return subtracted_particles;
  }




  std::vector<fastjet::PseudoJet> ConstituentSubtractor::get_background_proxies_from_ghosts(std::vector<fastjet::PseudoJet> const &ghosts,std::vector<double> const &ghosts_area) const{
    unsigned long nGhosts=ghosts.size();
    std::vector<fastjet::PseudoJet> proxies;
    std::vector<double> pt;
    std::vector<double> mtMinusPt;
    proxies.reserve(nGhosts);
    pt.reserve(nGhosts);
    mtMinusPt.reserve(nGhosts);
    if (_externally_supplied_rho_rhom){
      for (unsigned int j=0;j<nGhosts; ++j){
	pt.push_back(_rho*ghosts_area[j]);
	mtMinusPt.push_back(_rhom*ghosts_area[j]);
      }
    }
    else{
      for (unsigned int j=0;j<nGhosts; ++j) pt.push_back(_bge_rho->rho(ghosts[j])*ghosts_area[j]);
      if (_bge_rhom){
	/*	if (_rhom_from_bge_rhom){
#if FASTJET_VERSION_NUMBER >= 30100
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rhom->rho_m(ghosts[j])*ghosts_area[j]);
#else
	  throw(Error("ConstituentSubtractor:: _rhom_from_bge_rhom not allowed for FJ<3.1"));
#endif  // end of code specific to FJ >= 3.1
	} else {
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rhom->rho(ghosts[j])*ghosts_area[j]);
	  }*/


	// since FJ 3.1.0, some background estimators have an automatic internal calculation of rho_m
#if FASTJET_VERSION_NUMBER >= 30100
	// check if the BGE has internal support for rho_m
	if (_bge_rhom->has_rho_m()){
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rhom->rho_m(ghosts[j])*ghosts_area[j]);
	} else {
	  throw Error("ConstituentSubtractor: The provided background estimator for rho_m has no support to compute rho_m, and other option to get it is not available in ConstituentSubtractor.");
	}
#endif  // end of code specific to FJ >= 3.1
#if FASTJET_VERSION_NUMBER < 30100
	throw Error("ConstituentSubtractor: You provided a background estimator for rho_m, but this is not supported in ConstituentSubtractor when using fastjet version smaller than 30100.");
#endif
      }
      else if (_common_bge){
	// since FJ 3.1.0, some background estimators have an automatic internal calculation of rho_m
#if FASTJET_VERSION_NUMBER >= 30100
	// check if the BGE has internal support for rho_m
	if (_bge_rho->has_rho_m()){
	  for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(_bge_rho->rho_m(ghosts[j])*ghosts_area[j]);
	} else {
#endif  // end of code specific to FJ >= 3.1
	  BackgroundJetPtMDensity m_density;
	  JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
	  const FunctionOfPseudoJet<double> * orig_density = jmbge->jet_density_class();
	  jmbge->set_jet_density_class(&m_density);
	  for (unsigned int j=0;j<nGhosts; ++j)  mtMinusPt.push_back(jmbge->rho(ghosts[j])*ghosts_area[j]);
	  jmbge->set_jet_density_class(orig_density);
#if FASTJET_VERSION_NUMBER >= 30100
	}
#endif
      }
      else { // a single bge, only rho requested
	for (unsigned int j=0;j<nGhosts; ++j) mtMinusPt.push_back(1e-200);
#if FASTJET_VERSION_NUMBER >= 30100
	// In FJ3.1 and BGE with rho_m support, add a warning, similar to that in Subtractor
	double const rho_m_warning_threshold = 1e-5;
	if (_bge_rho->has_rho_m() && _bge_rho->rho_m()>rho_m_warning_threshold*_bge_rho->rho() && !_masses_to_zero && !_scale_fourmomentum){
	  _warning_unused_rhom.warn("ConstituentSubtractor:: Background estimator indicates non-zero rho_m, but the ConstituentSubtractor does not use rho_m information, nor the masses are set to zero, nor the 4-momentum is scaled. Consider calling set_common_bge_for_rho_and_rhom() to include the rho_m information; or call set_keep_original_masses(false) to set masses for all particles to zero; or call set_scale_fourmomentum to scale the fourmomentum.");
	}
#endif     
      }
    }


    fastjet::PseudoJet proxy(0,0,0,1);
    for (unsigned int j=0;j<nGhosts;++j){
      double mass_squared=pow(mtMinusPt[j]+pt[j],2)-pow(pt[j],2);
      double mass=0;
      if (mass_squared>0) mass=sqrt(mass_squared);
      proxy.reset_momentum_PtYPhiM(pt[j],ghosts[j].rap(),ghosts[j].phi(),mass);
      proxies.push_back(proxy);
    }
    return proxies;
  }
  
  


  double ConstituentSubtractor::_get_transformed_distance(double const &distance) const{
    double max_distance_transformed=-1;
    if (_distance==ConstituentSubtractor::deltaR) max_distance_transformed=pow(distance,2); 
    if (_distance==ConstituentSubtractor::angle) max_distance_transformed=-cos(distance);
    return max_distance_transformed;
  }



  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::do_subtraction(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const &backgroundProxies,std::vector<fastjet::PseudoJet> *remaining_backgroundProxies) const{
    unsigned int nBackgroundProxies=backgroundProxies.size();
    unsigned int nParticles=particles.size();
    double max_distance_transformed=this->_get_transformed_distance(_max_distance);

    // sort particles according to rapidity to achieve faster performance for the whole event subtraction
    std::vector<fastjet::PseudoJet>  particles_sorted=particles;
    std::sort(particles_sorted.begin(),particles_sorted.end(),ConstituentSubtractor::_rapidity_sorting);

    // get the kinematic variables for particles and background proxies in advance to achieve faster performance
    std::vector<double> particles_phi,particles_rap,particles_pt,particles_mt,particles_factor,particles_px_normalized,particles_py_normalized,particles_pz_normalized;
    particles_phi.reserve(nParticles);
    particles_rap.reserve(nParticles);
    particles_pt.reserve(nParticles);
    particles_mt.reserve(nParticles);
    particles_factor.reserve(nParticles);
    particles_px_normalized.reserve(nParticles);
    particles_py_normalized.reserve(nParticles);
    particles_pz_normalized.reserve(nParticles);

    double pt_factor=1;
    double polarAngle_factor=1;
    double nearby_hard_factor=1;
    double max_distance_from_hard_transformed=this->_get_transformed_distance(_nearby_hard_radius);
    for (unsigned int i=0;i<nParticles; ++i){
      particles_phi.push_back(particles_sorted[i].phi());
      particles_rap.push_back(particles_sorted[i].rap());
      particles_pt.push_back(particles_sorted[i].pt());
      particles_mt.push_back(particles_sorted[i].mt());
      if (fabs(_alpha)>1e-5) pt_factor=pow(particles_pt[i],_alpha);
      double momentum=sqrt(particles_sorted[i].pt2()+particles_sorted[i].pz()*particles_sorted[i].pz());
      if (fabs(_polarAngleExp)>1e-5) polarAngle_factor=pow(particles_pt[i]/momentum,_polarAngleExp);
      if (_distance==ConstituentSubtractor::angle){
	particles_px_normalized.push_back(particles_sorted[i].px()/momentum);
	particles_py_normalized.push_back(particles_sorted[i].py()/momentum);
	particles_pz_normalized.push_back(particles_sorted[i].pz()/momentum);
      }
      if (_use_nearby_hard){
	nearby_hard_factor=1;
	double distance_from_hard_transformed=-1;
	for (unsigned int ihard=0;ihard<_hard_proxies->size();++ihard){
	  if (_distance==ConstituentSubtractor::deltaR){
	    double deltaPhi=fabs(_hard_proxies->at(ihard).phi()-particles_phi[i]);
	    if (deltaPhi>pi) deltaPhi=twopi-deltaPhi;
	    double deltaRap=_hard_proxies->at(ihard).rap()-particles_rap[i];
	    distance_from_hard_transformed = deltaPhi*deltaPhi+deltaRap*deltaRap;
	  }
	  if (_distance==ConstituentSubtractor::angle){
	    distance_from_hard_transformed=-(particles_px_normalized[i]*_hard_proxies->at(ihard).px()+particles_py_normalized[i]*_hard_proxies->at(ihard).py()+particles_pz_normalized[i]*_hard_proxies->at(ihard).pz())/sqrt(_hard_proxies->at(ihard).pt2()*_hard_proxies->at(ihard).pt2()+_hard_proxies->at(ihard).pz()*_hard_proxies->at(ihard).pz());
	  }
	  if (distance_from_hard_transformed<=max_distance_from_hard_transformed){
	    nearby_hard_factor=_nearby_hard_factor;
	    break;
	  }
	}
      }
      particles_factor.push_back(pt_factor*polarAngle_factor*nearby_hard_factor);
    }

    std::vector<double> backgroundProxies_phi,backgroundProxies_rap,backgroundProxies_pt,backgroundProxies_mt,backgroundProxies_px_normalized,backgroundProxies_py_normalized,backgroundProxies_pz_normalized;
    backgroundProxies_phi.reserve(nBackgroundProxies);
    backgroundProxies_rap.reserve(nBackgroundProxies);
    backgroundProxies_pt.reserve(nBackgroundProxies);
    backgroundProxies_mt.reserve(nBackgroundProxies);
    backgroundProxies_px_normalized.reserve(nBackgroundProxies);
    backgroundProxies_py_normalized.reserve(nBackgroundProxies);
    backgroundProxies_pz_normalized.reserve(nBackgroundProxies);
    for (unsigned int j=0;j<nBackgroundProxies;++j){
      backgroundProxies_phi.push_back(backgroundProxies[j].phi());
      backgroundProxies_rap.push_back(backgroundProxies[j].rap());
      backgroundProxies_pt.push_back(backgroundProxies[j].pt());
      backgroundProxies_mt.push_back(backgroundProxies[j].mt());
      if (_distance==ConstituentSubtractor::angle){
	double momentum=sqrt(backgroundProxies[j].pt2()+backgroundProxies[j].pz()*backgroundProxies[j].pz());
	backgroundProxies_px_normalized.push_back(backgroundProxies[j].px()/momentum);
	backgroundProxies_py_normalized.push_back(backgroundProxies[j].py()/momentum);
	backgroundProxies_pz_normalized.push_back(backgroundProxies[j].pz()/momentum);
      }
    }

    // finding the rapidity range of particles for each ghost to achieve faster performance in the double loop over particles and proxies below
    double max_number_pairs=0;
    std::vector<unsigned int> backgroundProxies_minParticleIndex,backgroundProxies_maxParticleIndex;
    if (_use_max_distance && _distance==ConstituentSubtractor::deltaR && _ghosts_rapidity_sorted && !_ghost_selector){
      for (unsigned int j=0;j<_ghosts_rapidities.size();++j){
	unsigned int min=this->_find_index_after(_ghosts_rapidities[j]-_max_distance,particles_rap);
	unsigned int max=this->_find_index_before(_ghosts_rapidities[j]+_max_distance,particles_rap);
	for (int k=0;k<_n_ghosts_phi;++k){
 	  backgroundProxies_minParticleIndex.push_back(min);
	  backgroundProxies_maxParticleIndex.push_back(max);
	}
	max_number_pairs+=max-min;
      }
      max_number_pairs=max_number_pairs*_n_ghosts_phi*_max_distance/3.1415;
    }
    else{
      for (unsigned int j=0;j<nBackgroundProxies; ++j){
	backgroundProxies_minParticleIndex.push_back(0);
	backgroundProxies_maxParticleIndex.push_back(nParticles);      
      }
      if (_use_max_distance && _ghosts_constructed && !_ghost_selector){
	if (_distance==ConstituentSubtractor::deltaR) max_number_pairs=nParticles*nBackgroundProxies*_max_distance*_max_distance/4./_max_eta;
	if (_distance==ConstituentSubtractor::angle) max_number_pairs=nParticles*nBackgroundProxies*_max_distance/3.141593;
      }
      else max_number_pairs=nParticles*nBackgroundProxies;
    }
    
    
    // computation of the CS distances
    std::vector<std::pair<double,std::pair<int,int> > > CS_distances;  // storing three elements: the CS distance, and corresponding particle and proxy indexes
    CS_distances.reserve(max_number_pairs);

    bool skip_particles_outside_phi_range=false; // used for speed optimization
    if (_distance==ConstituentSubtractor::deltaR && _use_max_distance && _max_distance<twopi/2.*0.9999 && _ghosts_constructed && !_ghost_selector) skip_particles_outside_phi_range=true;
    double distance_transformed = 0;
    bool switched=false;
    double particle_phi_max=0;
    double particle_phi_min=0;
    for (unsigned int j=0;j<nBackgroundProxies; ++j){
      if (skip_particles_outside_phi_range){
	switched=false;
	particle_phi_max=backgroundProxies_phi[j]+_max_distance;
	particle_phi_min=backgroundProxies_phi[j]-_max_distance;
	if (particle_phi_max>twopi){
	  particle_phi_min=particle_phi_max-twopi;
	  particle_phi_max=backgroundProxies_phi[j]-_max_distance;
	  switched=true;
	}
	if (particle_phi_min<0){
	  particle_phi_max=particle_phi_min+twopi;
	  particle_phi_min=backgroundProxies_phi[j]+_max_distance;
	  switched=true;
	}
      }
      for (unsigned int i=backgroundProxies_minParticleIndex[j];i<backgroundProxies_maxParticleIndex[j];++i){
	if (_distance==ConstituentSubtractor::deltaR){
	  if (skip_particles_outside_phi_range) if ((switched && particles_phi[i]>particle_phi_min && particles_phi[i]<particle_phi_max) || (!switched && (particles_phi[i]<particle_phi_min || particles_phi[i]>particle_phi_max))) continue;  // speed optimization only, this line has no effect on the subtraction performance
	  double deltaPhi=fabs(backgroundProxies_phi[j]-particles_phi[i]);
	  if (deltaPhi>pi) deltaPhi=twopi-deltaPhi;
	  double deltaRap=backgroundProxies_rap[j]-particles_rap[i];
	  distance_transformed = deltaPhi*deltaPhi+deltaRap*deltaRap;
	}
	if (_distance==ConstituentSubtractor::angle){
	  distance_transformed=-(particles_px_normalized[i]*backgroundProxies_px_normalized[j]+particles_py_normalized[i]*backgroundProxies_py_normalized[j]+particles_pz_normalized[i]*backgroundProxies_pz_normalized[j]);
	}
	if (!_use_max_distance || distance_transformed<=max_distance_transformed){
	  double CS_distance = distance_transformed*particles_factor[i];
	  CS_distances.push_back(std::make_pair(CS_distance,std::make_pair(i,j))); // have tried to use emplace_back and tuple here - did not lead to any speed improvement
	}
      }
    }

    
    // Sorting of the CS distances
    std::sort(CS_distances.begin(),CS_distances.end(),ConstituentSubtractor::_function_used_for_sorting);


    // The iterative process. Here, only finding the fractions of pt or delta_m=mt-pt to be corrected. The actual correction of particles is done later.
    unsigned long nStoredPairs=CS_distances.size();

    std::vector<double> backgroundProxies_fraction_of_pt(nBackgroundProxies,1.);
    std::vector<double> particles_fraction_of_pt(nParticles,1.);
    std::vector<double> backgroundProxies_fraction_of_mtMinusPt(nBackgroundProxies,1.);
    std::vector<double> particles_fraction_of_mtMinusPt(nParticles,1.);
    for (unsigned long iindices=0;iindices<nStoredPairs;++iindices){
      int particle_index=CS_distances[iindices].second.first;
      int proxy_index=CS_distances[iindices].second.second;

      if (backgroundProxies_fraction_of_pt[proxy_index]>0 && particles_fraction_of_pt[particle_index]>0 && particles_pt[particle_index]>0 && backgroundProxies[proxy_index].pt()>0){
	double ratio_pt=particles_pt[particle_index]*particles_fraction_of_pt[particle_index]/backgroundProxies_pt[proxy_index]/backgroundProxies_fraction_of_pt[proxy_index];
	if (ratio_pt>1){
	  particles_fraction_of_pt[particle_index]*=1-1./ratio_pt;
	  backgroundProxies_fraction_of_pt[proxy_index]=-1;
	}
	else {
	  backgroundProxies_fraction_of_pt[proxy_index]*=1-ratio_pt;
	  particles_fraction_of_pt[particle_index]=-1;
	}
      }
      if (_do_mass_subtraction && backgroundProxies_fraction_of_mtMinusPt[proxy_index]>0 && particles_fraction_of_mtMinusPt[particle_index]>0 && particles_mt[particle_index]>particles_pt[particle_index] && backgroundProxies_mt[proxy_index]>backgroundProxies_pt[proxy_index]){
	double ratio_mtMinusPt=(particles_mt[particle_index]-particles_pt[particle_index])*particles_fraction_of_mtMinusPt[particle_index]/(backgroundProxies_mt[proxy_index]-backgroundProxies_pt[proxy_index])/backgroundProxies_fraction_of_mtMinusPt[proxy_index];
	if (ratio_mtMinusPt>1){
	  particles_fraction_of_mtMinusPt[particle_index]*=1-1./ratio_mtMinusPt;
	  backgroundProxies_fraction_of_mtMinusPt[proxy_index]=-1;
	}
	else{
	  backgroundProxies_fraction_of_mtMinusPt[proxy_index]*=1-ratio_mtMinusPt;
	  particles_fraction_of_mtMinusPt[particle_index]=-1;
	}
      }
    }


    // do the actual correction for particles:
    std::vector<fastjet::PseudoJet> subtracted_particles;
    PseudoJet subtracted_const;
    for (unsigned int i=0;i<nParticles; ++i){
      // Particles with zero pt and zero mass are removed by default. The user can decide to keep them by using function "set_remove_particles_with_zero_pt_and_mass" - then the mass is set to zero and the pt is set to very small value (1e-50).
      // Particles with zero pt and non-zero mass (i.e. massive particles at rest) are kept with very small pt (1e-50). The user can use function "set_remove_all_zero_pt_particles" to remove them or he/she can decide in his/her code if to remove or not to remove them.

      bool particle_pt_larger_than_zero=true;
      double corrected_pt=ConstituentSubtractorConstants::zero_pt;
      if (particles_fraction_of_pt[i]>0) corrected_pt=particles_pt[i]*particles_fraction_of_pt[i];

      if (corrected_pt<=ConstituentSubtractorConstants::zero_pt){
	if (_remove_all_zero_pt_particles) continue;
	else{
	  corrected_pt=ConstituentSubtractorConstants::zero_pt; 
	  particle_pt_larger_than_zero=false;
	}
      }
      
      
      if (_scale_fourmomentum){
	if (particle_pt_larger_than_zero) subtracted_const=particles_sorted[i]*particles_fraction_of_pt[i];
	else{
	  double scale=1;
	  if (corrected_pt<particles_pt[i]) scale=corrected_pt/particles_pt[i];
	  subtracted_const=particles_sorted[i]*scale;
	}
	if (subtracted_const.m()<=ConstituentSubtractorConstants::zero_mass && !particle_pt_larger_than_zero && _remove_particles_with_zero_pt_and_mass) continue;
      }
      else{
	bool particle_mass_larger_than_zero=(!_masses_to_zero && particles_sorted[i].m()>ConstituentSubtractorConstants::zero_mass);
	if (!particle_pt_larger_than_zero && !particle_mass_larger_than_zero && _remove_particles_with_zero_pt_and_mass) continue;

	double new_mass=ConstituentSubtractorConstants::zero_mass;
	if (_do_mass_subtraction){
	  if (particles_fraction_of_mtMinusPt[i]>0){
	    double subtracted_mtMinusPt=(particles_mt[i]-particles_pt[i])*particles_fraction_of_mtMinusPt[i];
	    double mass_squared=pow(corrected_pt+subtracted_mtMinusPt,2)-pow(corrected_pt,2);
	    if (mass_squared>0) new_mass=sqrt(mass_squared);
	  }
	}
	else if (!_masses_to_zero) new_mass=particles_sorted[i].m();
	
	if (new_mass<=ConstituentSubtractorConstants::zero_mass){
	  if (!particle_pt_larger_than_zero && _remove_particles_with_zero_pt_and_mass) continue;
	  else new_mass=ConstituentSubtractorConstants::zero_mass; 
	}
	

	//		std::cout << "before correction: " << particles_sorted[i].pt() << "  " << particles_sorted[i].eta() << "  " << particles_sorted[i].rap() << "  " << particles_sorted[i].m() << std::endl;  
	if (_fix_pseudorapidity){
	  double scale=corrected_pt/particles_pt[i];
	  subtracted_const.reset(particles_sorted[i].px()*scale,particles_sorted[i].py()*scale,particles_sorted[i].pz()*scale,sqrt(pow(corrected_pt,2)+pow(scale,2)*pow(particles_sorted[i].pz(),2)+pow(new_mass,2)));
	}
	else subtracted_const.reset_PtYPhiM(corrected_pt,particles_rap[i],particles_phi[i],new_mass);
	//	std::cout << "after correction: " << subtracted_const.pt() << "  " << subtracted_const.eta() << "  " << subtracted_const.rap() << "  " << subtracted_const.m() << std::endl;  
      }
      subtracted_const.set_user_index(particles_sorted[i].user_index());
      subtracted_particles.push_back(subtracted_const);
    }
    
    // get the remaining background proxies if requested:
    if (remaining_backgroundProxies){
      PseudoJet subtracted_proxy;
      for (unsigned int i=0;i<nBackgroundProxies; ++i){
	// keeping all background proxies
	bool proxy_pt_larger_than_zero=(backgroundProxies_fraction_of_pt[i]>0 && backgroundProxies_pt[i]>0);
	bool proxy_mtMinusPt_larger_than_zero=(!_masses_to_zero && backgroundProxies_fraction_of_mtMinusPt[i]>0 && backgroundProxies_mt[i]>backgroundProxies_pt[i]);
	
	double scale=1e-100;
	if (proxy_pt_larger_than_zero) scale=backgroundProxies_fraction_of_pt[i];
	if (_scale_fourmomentum) subtracted_proxy=backgroundProxies[i]*scale;
	else{
	  double new_mass=1e-150;
	  if (proxy_mtMinusPt_larger_than_zero){
	    if (_do_mass_subtraction){
	      double subtracted_mtMinusPt=(backgroundProxies_mt[i]-backgroundProxies_pt[i])*backgroundProxies_fraction_of_mtMinusPt[i];
	      double mass_squared=pow(scale*backgroundProxies_pt[i]+subtracted_mtMinusPt,2)-pow(scale*backgroundProxies_pt[i],2);
	      if (mass_squared>0) new_mass=sqrt(mass_squared);
	    }
	    else new_mass=backgroundProxies[i].m();
	  }
	  if (_fix_pseudorapidity) subtracted_proxy.reset(backgroundProxies[i].px()*scale,backgroundProxies[i].py()*scale,backgroundProxies[i].pz()*scale,sqrt(pow(scale,2)*(backgroundProxies[i].pt2()+pow(backgroundProxies[i].pz(),2))+pow(new_mass,2)));
	  else subtracted_proxy.reset_PtYPhiM(scale*backgroundProxies_pt[i],backgroundProxies_rap[i],backgroundProxies_phi[i],new_mass);
	}
	remaining_backgroundProxies->push_back(subtracted_proxy);
      }
    }
    
    return subtracted_particles;
  }
  
  
  
  void ConstituentSubtractor::set_keep_original_masses(){
    _masses_to_zero=false;
  }
  



  std::vector<fastjet::PseudoJet> ConstituentSubtractor::subtract_event_using_charged_info(std::vector<fastjet::PseudoJet> const &particles, double charged_background_scale, std::vector<fastjet::PseudoJet> const &charged_background, double charged_signal_scale, std::vector<fastjet::PseudoJet> const &charged_signal, double max_eta){
    if (fabs(_max_eta/max_eta-1)>1e-5) _ghosts_constructed=false;
    if (!_ghosts_constructed) this->construct_ghosts_uniformly(max_eta);
    _ghosts_rapidity_sorted=false;  // no speed optimization implemented for this function yet

    std::vector<fastjet::PseudoJet>  scaled_charged_all;
    std::vector<fastjet::PseudoJet>  scaled_charged_signal;
    std::vector<fastjet::PseudoJet>  scaled_charged_background;
    for (unsigned int i=0;i<charged_background.size();++i){
      if (fabs(charged_background[i].eta())>max_eta) continue;
      scaled_charged_background.push_back(charged_background_scale*charged_background[i]);
      scaled_charged_all.push_back(scaled_charged_background[scaled_charged_background.size()-1]);
    }
    for (unsigned int i=0;i<charged_signal.size();++i){
      if (fabs(charged_signal[i].eta())>max_eta) continue;
      scaled_charged_signal.push_back(charged_signal_scale*charged_signal[i]);
      scaled_charged_all.push_back(scaled_charged_signal[scaled_charged_signal.size()-1]);
    }
    std::vector<fastjet::PseudoJet> selected_particles;
    for (unsigned int i=0;i<particles.size();++i){
      if (fabs(particles[i].eta())<max_eta) selected_particles.push_back(particles[i]);
    }
    std::vector<fastjet::PseudoJet> *remaining_charged_background= new std::vector<fastjet::PseudoJet>;
    double maxDeltaR=this->get_max_distance();
    if (maxDeltaR<=0) maxDeltaR=0.5;
    this->set_max_distance(0.2);
    std::vector<fastjet::PseudoJet> subtracted_particles_using_scaled_charged_signal=this->do_subtraction(selected_particles,scaled_charged_signal); 
    std::vector<fastjet::PseudoJet> subtracted_particles_using_scaled_charged_all=this->do_subtraction(subtracted_particles_using_scaled_charged_signal,scaled_charged_background,remaining_charged_background);  // remaining neutral background particles
    std::vector<fastjet::PseudoJet> scaled_charged_background_used_for_subtraction=this->do_subtraction(scaled_charged_background,*remaining_charged_background); 
    _bge_rho= new fastjet::GridMedianBackgroundEstimator(max_eta, _grid_size_background_estimator);
    if (_do_mass_subtraction) this->set_common_bge_for_rho_and_rhom();
    _bge_rho->set_rescaling_class(_rescaling);
    _bge_rho->set_particles(subtracted_particles_using_scaled_charged_all);

    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    backgroundProxies.insert(backgroundProxies.end(), scaled_charged_background_used_for_subtraction.begin(), scaled_charged_background_used_for_subtraction.end());

    this->set_max_distance(maxDeltaR);
    std::vector<fastjet::PseudoJet> subtracted_particles=this->do_subtraction(selected_particles,backgroundProxies);
    delete remaining_charged_background; 
    delete _bge_rho;
    return subtracted_particles;
  }



  void ConstituentSubtractor::set_grid_size_background_estimator(double const &grid_size_background_estimator){
    _grid_size_background_estimator=grid_size_background_estimator;
  }



  void ConstituentSubtractor::set_common_bge_for_rho_and_rhom(bool value){
    if (value) this->set_common_bge_for_rho_and_rhom();
    else throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom: This function should be not used with false! Read the instructions for mass subtraction in the header file.");
  }


  void ConstituentSubtractor::set_common_bge_for_rho_and_rhom(){ 
    if (!_bge_rho) throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed when _bge_rho is not set!");
    if (_bge_rhom) throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed in the presence of an existing background estimator for rho_m.");
    if (_externally_supplied_rho_rhom) throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is not allowed when supplying externally the values for rho and rho_m.");
    
#if FASTJET_VERSION_NUMBER >= 30100
    if (!_bge_rho->has_rho_m()){
#endif
      JetMedianBackgroundEstimator *jmbge = dynamic_cast<JetMedianBackgroundEstimator*>(_bge_rho);
      if (!jmbge)	throw Error("ConstituentSubtractor::set_common_bge_for_rho_and_rhom() is currently only allowed for background estimators of JetMedianBackgroundEstimator type.");
#if FASTJET_VERSION_NUMBER >= 30100
    }
#endif
    _common_bge=true;
  }


  // setting this to true will result in rho_m being estimated using bge_rhom->rho_m() instead of bge_rhom->rho()
  /*  void ConstituentSubtractor::set_use_bge_rhom_rhom(bool value){
    if (!value){
      _rhom_from_bge_rhom=false;
      return;
    }
 
#if FASTJET_VERSION_NUMBER < 30100
    throw Error("ConnstituentSubtractor::use_rhom_from_bge_rhom() can only be used with FastJet >=3.1.");
#else	
    if (!_bge_rhom) throw Error("ConstituentSubtractor::use_rhom_from_bge_rhom() requires a background estimator for rho_m.");
    
    if (!(_bge_rhom->has_rho_m())) throw Error("ConstituentSubtractor::use_rhom_from_bge_rhom() requires rho_m support for the background estimator for rho_m.");
#endif	
    _rhom_from_bge_rhom=true; 
    _do_mass_subtraction=true;
    }*/


  void ConstituentSubtractor::set_do_mass_subtraction(){
    _do_mass_subtraction=true;
    _masses_to_zero=false;
  }



  void ConstituentSubtractor::description_common(std::ostringstream &descr) const{
    if ( _externally_supplied_rho_rhom){
      descr << "       Using externally supplied rho = " << _rho << " and rho_m = " << _rhom << std::endl;
    } else {
      if (_bge_rhom && _bge_rho) {
	descr << "       Using rho estimation: " << _bge_rho->description() << std::endl;
	descr << "       Using rho_m estimation: " << _bge_rhom->description() << std::endl;
      } else {
	if (_bge_rho) descr << "       Using rho estimation: " << _bge_rho->description() << std::endl;
	else descr << "       No externally supplied rho, nor background estimator" << std::endl;
      }
    }  

    if (_do_mass_subtraction){
      descr << "       The mass part (delta_m) will be also corrected." << std::endl;
      if (_common_bge) descr << "       using the same background estimator for rho_m as for rho" << std::endl;
      else  descr << "       using different background estimator for rho_m as for rho" << std::endl;
    }
    else if (_masses_to_zero) descr << "       The masses of all particles will be set to zero." << std::endl;
    else if (_scale_fourmomentum) descr << "       The masses will be corrected by scaling the whole 4-momentum." << std::endl;
    else descr << "       The original mass of the particles will be kept." << std::endl;

    if (!_scale_fourmomentum){
      if (_fix_pseudorapidity) descr << "       The pseudo-rapidity of the particles will be kept unchanged (not rapidity)." << std::endl;
      else descr << "       The rapidity of the particles will be kept unchanged (not pseudo-rapidity)." << std::endl;
    }

    if (_use_nearby_hard) descr << "       Using information about nearby hard proxies with parameters _nearby_hard_radius=" << _nearby_hard_radius << " and _nearby_hard_factor=" << _nearby_hard_factor << std::endl;
    else descr << "       The information about nearby hard proxies will not be used." << std::endl;
  }





  std::string ConstituentSubtractor::description() const{
    std::ostringstream descr;
    descr << std::endl << "Description of fastjet::ConstituentSubtractor which can be used for event-wide or jet-by-jet correction:" << std::endl;
    this->description_common(descr);
    descr << "       Using parameters: max_distance = " << _max_distance << "   alpha = " << _alpha << std::endl;
    return descr.str();
  }



  void ConstituentSubtractor::set_distance_type(Distance distance){
    _distance=distance;
  }


  void ConstituentSubtractor::set_max_distance(double max_distance){
    if (max_distance>0){
      _use_max_distance=true;
      _max_distance=max_distance;
    }
    else _use_max_distance=false; 
  }



  void ConstituentSubtractor::set_max_standardDeltaR(double max_distance){
    this->set_max_distance(max_distance);
  }



  double ConstituentSubtractor::get_max_distance(){
    return _max_distance;
  }


  void ConstituentSubtractor::set_alpha(double alpha){
    _alpha=alpha;
  }

  void ConstituentSubtractor::set_polarAngleExp(double polarAngleExp){
    _polarAngleExp=polarAngleExp;
  }

 
  void ConstituentSubtractor::set_ghost_area(double ghost_area){
    _ghost_area=ghost_area;
    this->clear_ghosts();
  }


  void ConstituentSubtractor::set_max_eta(double max_eta){
    _max_eta=max_eta;
  }


  void ConstituentSubtractor::set_fix_pseudorapidity(){
    _fix_pseudorapidity=true;
  }


  void ConstituentSubtractor::set_scale_fourmomentum(){
    _scale_fourmomentum=true;
    _masses_to_zero=false;
  }


  void ConstituentSubtractor::set_remove_particles_with_zero_pt_and_mass(bool value){
    _remove_particles_with_zero_pt_and_mass=value;
  }


  void ConstituentSubtractor::set_remove_all_zero_pt_particles(bool value){
    _remove_all_zero_pt_particles=value;
  }


  bool ConstituentSubtractor::_function_used_for_sorting(std::pair<double,std::pair<int,int> >  const &i,std::pair<double,std::pair<int,int> >  const &j){
    return (i.first < j.first);
  }
  
  bool ConstituentSubtractor::_rapidity_sorting(fastjet::PseudoJet const &i,fastjet::PseudoJet  const &j){
    return (i.rap() < j.rap());
  }

  unsigned int ConstituentSubtractor::_find_index_after(double const &value, std::vector<double> const &vec) const{
    int size=vec.size();
    if (size==0) return -1;
    int nIterations=log(size)/log(2)+2;
    unsigned int lowerBound=0;
    unsigned int upperBound=size-1;
    if (value<=vec[0]) return 0;
    if (value>vec[size-1]) return size; // if the value is larger than the maximal possible rapidity, than setting min to size - then also max is set to size, and nothing will be run in the loop
      for (int i=0;i<nIterations;++i){
	unsigned int test=(upperBound+lowerBound)/2;
	if (value>vec[test]){
	  if (value<=vec[test+1]) return test+1;
	  lowerBound=test;
	}
	else{
	  if (value>vec[test-1]) return test;
	  upperBound=test;
	}
      }
    return lowerBound;
  }


  unsigned int ConstituentSubtractor::_find_index_before(double const &value, std::vector<double> const &vec) const{
    int size=vec.size();
    if (size==0) return -1;
    int nIterations=log(size)/log(2)+1;
    unsigned int lowerBound=0;
    unsigned int upperBound=size-1;
    if (value<vec[0]) return 0;  // if the value is lower than the minimal possible rapidity, than setting max to 0 - then also min is set to 0, and nothing will be run in the loop 
    if (value>=vec[size-1]) return size;  // it is higher by one to account for the "<" comparison in the for loop
    for (int i=0;i<nIterations;++i){
      unsigned int test=(upperBound+lowerBound)/2;
      if (value>=vec[test]){
	if (value<vec[test+1]) return test+1;  // it is higher by one to account for the "<" comparison in the for loop
	lowerBound=test;
      }
      else{
	if (value>=vec[test-1]) return test;  // it is higher by one to account for the "<" comparison in the for loop
	upperBound=test;
      }
    }
    return upperBound+1;
  }
  



  void ConstituentSubtractor::clear_ghosts(){
    _ghosts.clear();
    _ghosts_rapidities.clear();
    _ghosts_area.clear();
    _ghosts_rapidity_sorted=false;
    _ghosts_constructed=false;
  }


  void ConstituentSubtractor::construct_ghosts_uniformly(double max_eta){
    this->clear_ghosts();
    _max_eta=max_eta;
    double a=sqrt(_ghost_area);
    _n_ghosts_phi=(2*3.14159265/a+0.5); // rounding 
    int n_ghosts_rap=(2*max_eta/a+0.5); // rounding 
    _grid_size_phi=2*3.14159265/(double)_n_ghosts_phi;
    _grid_size_rap=2*max_eta/(double)n_ghosts_rap;
    double used_ghost_area=_grid_size_phi*_grid_size_rap;
    fastjet::PseudoJet ghost(0,0,0,1);
    for (int iRap=0;iRap<n_ghosts_rap;++iRap){
      double rapidity=_grid_size_rap*(iRap+0.5)-max_eta;
      _ghosts_rapidities.push_back(rapidity);
      for (int iPhi=0;iPhi<_n_ghosts_phi;++iPhi){
	ghost.reset_momentum_PtYPhiM(1,rapidity,_grid_size_phi*(iPhi+0.5),1e-200);
	if (_ghost_selector){
	  if (!_ghost_selector->pass(ghost)) continue;
	} 
	_ghosts.push_back(ghost);
	_ghosts_area.push_back(used_ghost_area);
      }
    }
    _ghosts_rapidity_sorted=true; // the ghosts are now sorted according to rapidity. This variable needs to be set to true to be able to use faster algorithm in "do_subtraction".
    _ghosts_constructed=true;
  }


  std::vector<fastjet::PseudoJet>  ConstituentSubtractor::get_ghosts(){
    return _ghosts;
  }


  std::vector<double>  ConstituentSubtractor::get_ghosts_area(){
    return _ghosts_area;
  }



  void ConstituentSubtractor::set_ghost_selector(fastjet::Selector* selector){
    _ghost_selector=selector;
    this->clear_ghosts();
  }

  void ConstituentSubtractor::set_particle_selector(fastjet::Selector* selector){
    _particle_selector=selector;
  }

  void ConstituentSubtractor::set_rescaling(fastjet::FunctionOfPseudoJet<double> *rescaling){
    _rescaling=rescaling;
  }


  void ConstituentSubtractor::set_use_nearby_hard(double const &nearby_hard_radius, double const &nearby_hard_factor){
    _nearby_hard_radius=nearby_hard_radius;
    _nearby_hard_factor=nearby_hard_factor;
    if (_nearby_hard_radius>0) _use_nearby_hard=true;
    else _use_nearby_hard=false;
  }

} // namespace contrib


FASTJET_END_NAMESPACE
