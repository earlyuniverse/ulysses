// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2020 The HepMC collaboration (see AUTHORS for details)
//
#ifndef _TauolaHepMC3Particle_h_included_
#define _TauolaHepMC3Particle_h_included_

#warning "HepMC3 interface is available in the latest version of TAUOLA, see http://tauolapp.web.cern.ch/tauolapp/. This interface will be removed in the future HepMC3 versions."

/**
 * @class TauolaHepMC3Particle
 *
 * @brief Interface to GenParticle objects
 *
 * This class implements the virtual methods of
 * TauolaParticle. In this way it provides an
 * interface between the generic TauolaParticle class
 * and a GenParticle object.
 *
 * This code is licensed under GNU General Public Licence.
 * For more informations, see: http://www.gnu.org/licenses/
 */

#include <iostream>
#include <vector>

#include "HepMC3/GenParticle.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"

//#include "DecayList.h"
#include "Tauola/TauolaParticle.h"
#include "Tauola/f_Decay.h"
namespace Tauolapp
{
using namespace HepMC3;

class TauolaHepMC3Particle: public TauolaParticle {

public:
    /** General constructor */
    TauolaHepMC3Particle();

    ~TauolaHepMC3Particle();

    /** Constructor which keeps a pointer to the GenParticle*/
    TauolaHepMC3Particle(GenParticlePtr particle);

    /** Constructor which creates a new GenParticle and
         sets the properties pdg_id, statu and mass. */
    TauolaHepMC3Particle(int pdg_id, int status, double mass);

    /** Returns the GenParticlePtr */
    GenParticlePtr getHepMC3();

    /** Remove the decay branch from the event record and reset the particle status code to stable. */
    void undecay();

    /** Set the mothers of this particle via a vector of TauolaParticle*/
    void setMothers(std::vector<TauolaParticle*> mothers);

    /** Set the daughters of this particle via a vector of TauolaParticle*/
    void setDaughters(std::vector<TauolaParticle*> daughters);

    /** Returns the mothers of this particle via a vector of TauolaParticle */
    std::vector<TauolaParticle*> getMothers();

    /** Returns the daughters of this particle via a vector of TauolaParticle */
    std::vector<TauolaParticle*> getDaughters();

    /** Set the PDG ID code of this particle */
    void setPdgID(int pdg_id);

    /** Set the status of this particle */
    void setStatus(int statu);

    /** Set the mass of this particle */
    void setMass(double mass);

    /** Get the PDG ID code of this particle */
    int getPdgID();

    /** Get the status of this particle */
    int getStatus();

    /** Get the barcode of this particle */
    int getBarcode();

    /** Check that the 4 momentum in conserved at the vertices producing
        and ending this particle */
    void checkMomentumConservation();

    /** Overriding of TauolaParticle decayEndgame method.
        Converts the momentum and length units
        and sets the vector (X,T) position */
    void decayEndgame();

    /** Create a new particle of type TauolaHepMC3Particle, with the given
        properties. The new particle bares no relations to this
        particle, but it provides a way of creating a instance of
        this derived class. eg. createNewParticle() is used inside
        filhep_() so that a TauolaHepMC3Particle can be created without
        the method having explicit knowledge of the TauolaHepMC3Particle
        class */
    TauolaHepMC3Particle * createNewParticle(int pdg_id, int status, double mass,
            double px, double py,
            double pz, double e);

    /** Print some information about this particle to standard output */
    void print();

    /** Returns the px component of the four vector*/
    double getPx();

    /** Returns the py component of the four vector */
    double getPy();

    /** Returns the pz component of the four vector */
    double getPz();

    /** Returns the energy component of the four vector */
    double getE();

    /** Set the px component of the four vector */
    void setPx( double px );

    /** Set the px component of the four vector */
    void setPy( double py );

    /** Set the pz component of the four vector */
    void setPz( double pz );

    /** Set the energy component of the four vector */
    void setE( double e );


private:

    /** Sets the position for whole decay tree starting from given particle */
    void recursiveSetPosition(GenParticlePtr p,FourVector pos);

    /** A pointer to the GenParticle particle */
    GenParticlePtr m_particle;

    /** A list of mothers */
    std::vector<TauolaParticle*> m_mothers;

    /** A list of daughters */
    std::vector<TauolaParticle*> m_daughters;

    /** List to keep track of new particles which have been
        created from this one, so we can call their destructor later */
    std::vector<TauolaParticle*> m_created_particles;

};

} // namespace Tauolapp
#endif
