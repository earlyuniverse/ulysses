#!/usr/bin/env python
"""Simple example to illustrate use of FastJet areas from Python.

For this script to work, make sure that the installation location for
the fastjet python module (cf. fastjet-config --pythonpath) is
included in your PYTHONPATH environment variable.

"""

from __future__ import print_function
from fastjet import *

def main():
    #----------------------------------------------------------------------
    # read the event
    event = read_event("../data/single-event.dat")
    print("Event has {0} particles".format(len(event)))
    
    #----------------------------------------------------------------------
    # cluster the event
    jet_def = JetDefinition(antikt_algorithm, 0.4)
    area_def = AreaDefinition(active_area, GhostedAreaSpec(5.0))
    cs = ClusterSequenceArea(event, jet_def, area_def)
    jets = SelectorPtMin(5.0)(sorted_by_pt(cs.inclusive_jets()))

    print("jet def:", jet_def)
    print("area def:", area_def)
    print("#-------------------- initial jets --------------------")
    print_jets(jets)

    #----------------------------------------------------------------------
    # estimate the background
    maxrap       = 4.0
    grid_spacing = 0.55
    gmbge = GridMedianBackgroundEstimator(maxrap, grid_spacing)
    gmbge.set_particles(event)
    print("#-------------------- background properties --------------------")
    print("rho   = ", gmbge.rho())
    print("sigma = ", gmbge.sigma())
    print()
    
    #----------------------------------------------------------------------
    # subtract the jets
    subtractor = Subtractor(gmbge)
    subtracted_jets = subtractor(jets)
    print("#-------------------- subtracted jets --------------------")
    print_jets(subtracted_jets)
    

#----------------------------------------------------------------------
def read_event(filename):
    f = open(filename, 'r')
    event = []
    for line in f:
        p = line.split()
        event.append(PseudoJet(float(p[0]),float(p[1]),float(p[2]),float(p[3])));

    return event

#----------------------------------------------------------------------
def print_jets(jets):
    print("{0:>5s} {1:>10s} {2:>10s} {3:>10s} {4:>10s}".format(
        "jet #", "pt", "rap", "phi", "area"))

    for ijet in range(len(jets)):
        jet = jets[ijet]
        print("{0:5d} {1:10.3f} {2:10.4f} {3:10.4f} {4:10.4f}".format(
            ijet, jet.pt(), jet.rap(), jet.phi(), jet.area()))
    

if __name__ == '__main__':
    main()


