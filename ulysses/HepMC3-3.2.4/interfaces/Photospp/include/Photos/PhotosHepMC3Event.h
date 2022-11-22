// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
#ifndef _PhotosHepMC3Event_h_included_
#define _PhotosHepMC3Event_h_included_

#warning "HepMC3 interface is available in the latest version of PHOTOS, see http://photospp.web.cern.ch/photospp/. This interface will be removed in the future HepMC3 versions."

/**
 * @class PhotosHepMC3Event
 *
 * @brief Interface to GenEvent objects
 *
 * This class implements the virtual methods of
 * PhotosEvent. In this way it provides an
 * interface between the generic PhotosEvent class
 * and a GenEvent object.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 *
 * This code is licensed under GNU General Public Licence.
 * For more informations, see: http://www.gnu.org/licenses/
 */

#include <vector>
#include "HepMC3/GenEvent.h"
#include "PhotosEvent.h"
#include "PhotosParticle.h"

namespace Photospp
{
using namespace HepMC3;
class PhotosHepMC3Event : public PhotosEvent
{
public:
    ~PhotosHepMC3Event();

    /** Constructor which keeps a pointer to the GenEvent*/
    PhotosHepMC3Event(GenEvent * event);

    /** Returns the GenEvent */
    GenEvent * getEvent();

    /** Returns the list of particles */
    std::vector<PhotosParticle*> getParticleList();

    /** Prints event summary */
    void print();
private:
    /** The event */
    GenEvent * m_event;
    /** Particle list */
    std::vector<PhotosParticle *> particles;
};

} // namespace Photospp
#endif
