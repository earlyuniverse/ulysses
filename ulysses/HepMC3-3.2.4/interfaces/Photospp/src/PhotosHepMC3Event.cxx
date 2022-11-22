// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#include <vector>
#include "Photos/PhotosHepMC3Particle.h"
#include "Photos/PhotosHepMC3Event.h"
#include "Photos/Log.h"

#include "HepMC3/Print.h"

namespace Photospp
{
PhotosHepMC3Event::PhotosHepMC3Event(GenEvent * event)
{
    m_event=event;
    for(auto p: m_event->particles() )
    {
        PhotosParticle *particle = new PhotosHepMC3Particle(p);
        particles.push_back(particle);
    }
}

PhotosHepMC3Event::~PhotosHepMC3Event()
{
    while(particles.size())
    {
        PhotosParticle *p = particles.back();
        particles.pop_back();
        if(p) delete p;
    }
}

GenEvent * PhotosHepMC3Event::getEvent()
{
    return m_event;
}

void PhotosHepMC3Event::print()
{
    if(!m_event) return;
    Print::listing(*m_event);
}

vector<PhotosParticle*> PhotosHepMC3Event::getParticleList()
{
    return particles;
}

} // namespace Photospp
