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

#include "IterativeConstituentSubtractor.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

  ///
  /// default constructor
  IterativeConstituentSubtractor::IterativeConstituentSubtractor(){
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
    _ghost_removal=true;
    _ghost_area=0.01;
    _grid_size_phi=-1;
    _grid_size_rap=-1;
    _ghosts_constructed=false;
    _ghosts_rapidity_sorted=false;
    _n_ghosts_phi=-1;
    _max_eta=-1;
    _rescaling=0;
    _use_nearby_hard_iterative=false;
    _fix_pseudorapidity=false;
    _scale_fourmomentum=false;
    _hard_proxies=0;
    _ghost_selector=0;
    _particle_selector=0;
  }


  void IterativeConstituentSubtractor::initialize(){
    if (_max_distances.size()==0) throw Error("IterativeConstituentSubtractor::initialize(): The vector for max_distances is empty. It should be specified before using the function initialize.");
    this->_initialize_common();
  }



  std::vector<fastjet::PseudoJet> IterativeConstituentSubtractor::subtract_event(std::vector<fastjet::PseudoJet> const &particles, double max_eta){
    throw Error("IterativeConstituentSubtractor::subtract_event(): This version of subtract_event should not be used. Use the version subtract_event(std::vector<fastjet::PseudoJet> const &particles)");
    std::vector<fastjet::PseudoJet> subtracted_particles;
    return subtracted_particles;
  }



  std::vector<fastjet::PseudoJet> IterativeConstituentSubtractor::subtract_event(std::vector<fastjet::PseudoJet> const &particles, std::vector<fastjet::PseudoJet> const *hard_proxies){
    bool original_value_of_ghosts_rapidity_sorted=_ghosts_rapidity_sorted;
    if (_use_nearby_hard_iterative){
      if (hard_proxies) _hard_proxies=hard_proxies;
      else throw Error("IterativeConstituentSubtractor::subtract_event: It was requested to use closeby hard proxies but they were not provided in this function!");
    }
    else if (hard_proxies) throw Error("IterativeConstituentSubtractor::subtract_event: Hard proxies were provided but the set_use_hard_proxies function was not used before initialization. It needs to be called before initialization!");

    std::vector<fastjet::PseudoJet> backgroundProxies=this->get_background_proxies_from_ghosts(_ghosts,_ghosts_area);
    std::vector<fastjet::PseudoJet> subtracted_particles,unselected_particles;
    for (unsigned int iparticle=0;iparticle<particles.size();++iparticle){
      if (std::abs(particles[iparticle].eta())>_max_eta) continue;
      if (particles[iparticle].pt()<ConstituentSubtractorConstants::zero_pt && _remove_all_zero_pt_particles) continue;
      if (particles[iparticle].pt()<ConstituentSubtractorConstants::zero_pt && (_masses_to_zero || particles[iparticle].m()<ConstituentSubtractorConstants::zero_mass) && _remove_particles_with_zero_pt_and_mass) continue;
      if (_particle_selector){
        if (_particle_selector->pass(particles[iparticle])) subtracted_particles.push_back(particles[iparticle]);
        else unselected_particles.push_back(particles[iparticle]);
      }
      else subtracted_particles.push_back(particles[iparticle]);
    }
    for (unsigned int iteration=0;iteration<_max_distances.size();++iteration){
      this->set_max_distance(_max_distances[iteration]);
      this->set_alpha(_alphas[iteration]);
      if (_use_nearby_hard_iterative) this->set_use_nearby_hard(_nearby_hard_radii[iteration],_nearby_hard_factors[iteration]);
      std::vector<fastjet::PseudoJet> *remaining_backgroundProxies=0;
      if (iteration!=_max_distances.size()-1) remaining_backgroundProxies=new std::vector<fastjet::PseudoJet>; // do not need to estimate the remaining background for the last iteration
      subtracted_particles=this->do_subtraction(subtracted_particles,backgroundProxies,remaining_backgroundProxies);

      if (iteration==_max_distances.size()-1){
	continue; // do not need to estimate the remaining background for the last iteration
      }

      double background_pt=0,background_mt=0, remaining_background_pt=0,remaining_background_mt=0;
      int backgroundProxies_original_size=backgroundProxies.size();
      for (int i=backgroundProxies_original_size-1;i>=0;--i){
	remaining_background_pt+=(remaining_backgroundProxies->at(i)).pt();
	remaining_background_mt+=(remaining_backgroundProxies->at(i)).mt();
	if (_ghost_removal && (remaining_backgroundProxies->at(i)).pt()>1e-10){
	  backgroundProxies[i]=backgroundProxies[backgroundProxies.size()-1];
	  backgroundProxies.pop_back();
	}
	else{
	  background_pt+=backgroundProxies[i].pt();
	  background_mt+=backgroundProxies[i].mt();
	}
      }
      delete remaining_backgroundProxies;
      if (_ghost_removal) _ghosts_rapidity_sorted=false; // After erasing some elements from vector backgroundProxies, we cannot use the speed optimization in do_subtraction function, therefore setting _ghosts_rapidity_sorted to false.
      for (unsigned int i=0;i<backgroundProxies.size();++i){
	double pt= backgroundProxies[i].pt()*remaining_background_pt/background_pt;
	double mtMinusPt=0;
	if (background_mt>background_pt+1e-20) mtMinusPt=(backgroundProxies[i].mt()-backgroundProxies[i].pt())*(remaining_background_mt-remaining_background_pt)/(background_mt-background_pt);
	double mass=0;
	if (mtMinusPt>1e-20) mass=sqrt(pow(mtMinusPt+pt,2)-pow(pt,2));
	backgroundProxies[i].reset_momentum_PtYPhiM(pt,backgroundProxies[i].rap(),backgroundProxies[i].phi(),mass);
      } 
    }    
    _ghosts_rapidity_sorted=original_value_of_ghosts_rapidity_sorted; // Setting _ghosts_rapidity_sorted to its original value
    if (_particle_selector) subtracted_particles.insert(subtracted_particles.end(), unselected_particles.begin(), unselected_particles.end());
    return subtracted_particles;
  }


  std::string IterativeConstituentSubtractor::description() const{
    std::ostringstream descr;
    descr << std::endl << "Description of fastjet::IterativeConstituentSubtractor:" << std::endl;
    this->description_common(descr);
    descr << "       IterativeConstituentSubtractor parameters: " << std::endl;
    for (unsigned int iteration=0;iteration<_max_distances.size();++iteration){
      descr << "            Iteration " << iteration+1 << ":  max_distance = " << _max_distances[iteration] << "  alpha = " << _alphas[iteration] << std::endl;
    }
    return descr.str();
  }



  void IterativeConstituentSubtractor::set_parameters(std::vector<double> const &max_distances, std::vector<double> const &alphas){
    if (max_distances.size()!=alphas.size()) throw Error("IterativeConstituentSubtractor::set_parameters(): the provided vectors have different size. They should have the same size.");
    if (max_distances.size()==0 || alphas.size()==0) throw Error("IterativeConstituentSubtractor::set_parameters(): One of the provided vectors is empty. They should be not empty.");
    _max_distances=max_distances;
    _alphas=alphas;
  }


  void IterativeConstituentSubtractor::set_nearby_hard_parameters(std::vector<double> const &nearby_hard_radii, std::vector<double> const &nearby_hard_factors){
    if (nearby_hard_radii.size()!=nearby_hard_factors.size()) throw Error("IterativeConstituentSubtractor::set_use_nearby_hard(): the provided vectors have different size. They should have the same size.");
    if (nearby_hard_radii.size()==0 || nearby_hard_factors.size()==0) throw Error("IterativeConstituentSubtractor::set_use_nearby_hard(): One of the provided vectors is empty. They should be not empty.");
    _nearby_hard_radii=nearby_hard_radii;
    _nearby_hard_factors=nearby_hard_factors;
    _use_nearby_hard_iterative=true;
  }



  void IterativeConstituentSubtractor::set_ghost_removal(bool ghost_removal){
    _ghost_removal=ghost_removal;
  }



} // namespace contrib


FASTJET_END_NAMESPACE
