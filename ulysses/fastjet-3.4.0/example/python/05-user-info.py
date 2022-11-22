#!/usr/bin/env python
"""Illustration of the assignment of pythonic user information to PseudoJets:

  - the script has a ParticleInfo class, to store information about particles

  - the read_event(...) function creates a ParticleInfo for each
    particle and assigns it to the corresponding PseudoJet, using the
    PseudoJet.set_python_info(...) call.

  - the print_jets(...) function gets the jet constituents, examines
    the ParticleInfo for each one and uses it to determine additional
    information about each jet. It uses the PseudoJet.python_info(...)
    call.

For this script to work, make sure that the installation location for
the fastjet python module (cf. fastjet-config --pythonpath) is
included in your PYTHONPATH environment variable.

"""
from __future__ import print_function

import fastjet as fj
import gzip

def main():

    # get the banner out of the way early on
    fj.ClusterSequence.print_banner()
    print()

    # set up our jet definition and a jet selector
    jet_def = fj.JetDefinition(fj.antikt_algorithm, 0.4)
    selector = fj.SelectorPtMin(15.0) & fj.SelectorAbsRapMax(4.5)
    print("jet definition is:",jet_def)
    print("jet selector is:", selector,"\n")

    filename = '../data/Pythia-dijet-ptmin100-lhc-pileup-1ev.dat'
    f = open(filename,'r')
    #filename = '/Users/gsalam/work/fastjet/data/Pythia-PtMin50-LHC-10kev.dat.gz'
    #f = gzip.GzipFile(filename,'rb')
    
    # get the event
    iev = 0
    while True:
        event = read_event(f)
        iev += 1
        if (len(event) == 0): break
        
        # cluster it
        jets = selector(jet_def(event))

        # print some info
        npileup = 0
        for p in event:
            if (p.python_info().subevent_index > 0): npileup += 1
        print("Event {0} has {1} particles (of which {2} from pileup)".format(iev, len(event), npileup))
        print_jets(jets)

#----------------------------------------------------------------------
def print_jets(jets):
    print("{0:>5s} {1:>10s} {2:>10s} {3:>10s} {4:>12s} {5:>12s} {6:>12s}".format(
        "jet #", "pt", "rap", "phi", "primary pt", "N particles", "N photons"))

    for ijet in range(len(jets)):
        jet = jets[ijet]

        # figure out how many particles and how many photons the jet contains
        # and how much pt comes from the primary vertex
        constituents = jet.constituents()
        nphotons = 0
        primary_pt = 0
        for c in constituents:
            if (c.python_info().pdg_id == 22): nphotons += 1
            if (c.python_info().subevent_index <= 0): primary_pt += c.pt()
            
        print("{0:5d} {1:10.3f} {2:10.4f} {3:10.4f} {4:10.3f} {5:12d} {6:12d}".format(
            ijet, jet.pt(), jet.rap(), jet.phi(), primary_pt, len(constituents), nphotons))
        
#----------------------------------------------------------------------            
class ParticleInfo(object):
    """illustrative class for use in assigning pythonic user information
    to a PseudoJet.
    """
    def __init__(self, subevent_index, particle_index, pdg_id=0):
        self.subevent_index = subevent_index
        self.particle_index = particle_index
        self.pdg_id = pdg_id

    def set_pdg_id(self, pdg_id):
        self.pdg_id = pdg_id

    def __str__(self):
        return "subevent_index={0}, particle_index={1}, pdg_id={2}".format(
            self.subevent_index, self.particle_index, self.pdg_id)
        
#----------------------------------------------------------------------
event_index = 0
def read_event(file_or_filename):
    """
Routine that can take either an existing opened file object, or a
filename (which it will open itself) and then reads an event from that
file. An event is deemed to end when the file comes to an end or when
the reader encounters the string "#END".

The event is converted to a python list of PseudoJets
    """

    # open the file if necessary
    if (isinstance(file_or_filename,str)) : f = open(file_or_filename, 'r')
    else                                  : f = file_or_filename

    # create an empty list
    event = []
    particle_index = 0
    subevent_index = -1
    while True:
        line = f.readline()
        # exit if the file has come to an end
        if   (not line): break
        # or if we reach the string "#END"
        if   (len(line) >=4 and line[0:4] == '#END'): break

        # #SUBSTART means that we start a new subevent
        if (len(line)>=9 and line[0:9] == '#SUBSTART'):
            subevent_index += 1
            continue
        
        # ignore comment lines or empty lines
        elif (line[0] == '#' or len(line) <= 1): continue

        # assume we have a good line and split it into px, py, pz, E
        p = line.split()
        # create the PseudoJet
        pj = fj.PseudoJet(float(p[0]),float(p[1]),float(p[2]),float(p[3]))

        # then create its user info
        info = ParticleInfo(subevent_index, particle_index)
        # optionally with a pdg_id (if it was in the file)
        if (len(p) > 4): info.set_pdg_id(int(p[4]))
        
        # finally assign the user info and append the particle the event
        pj.set_python_info(info)
        event.append(pj);
        
        # don't forget to increment it
        particle_index += 1

    global event_index
    event_index += 1
    return event
    
if __name__ == '__main__':
    main()
