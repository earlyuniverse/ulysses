#!/usr/bin/env python
"""Simple example to illustrate use of various FastJet tools in python.

For this script to work, make sure that the installation location for
the fastjet python module (cf. fastjet-config --pythonpath) is
included in your PYTHONPATH environment variable.

"""

from __future__ import print_function
from fastjet import *
import re
import copy

def main():
    event = read_event("../data/single-event.dat")
    print("Event has {} particles and is of type {}".format(len(event), type(event)))
    
    jet_def = JetDefinition(antikt_algorithm, 0.4)
    jets = SelectorPtMin(5.0)(jet_def(event))
    print(jets)
    print(len(jets))
    for jet in jets:
        print(jet.pt(), jet.rap())

    # select the hardest jet to play with
    print("---------- Selecting hardest jet ----------")
    hard_jet = SelectorNHardest(1)(jets)[0]
    print("Hard jet: ", hard_jet)
    print()
    
    # try boosting wrt the 2nd jet
    print("---------- testing (Un)boost ----------")
    prest = jets[1]
    boost=Unboost(prest)
    boosted_hard = boost(hard_jet)
    print("After (un)boosting ", boosted_hard)

    # try boosting all jets
    boosted_all = boost(jets)
    print("After (un)boosting all: ", boosted_all)
    print("Hardest became ", boosted_all[0])
    print("Prest became ", boosted_all[1])
    print()

    # try reclustering and run MDT
    print("---------- testing Recluster + MDT ----------")
    
    recluster = Recluster(cambridge_algorithm)
    print("Using recluster: ",recluster)
    rec_hard= recluster(hard_jet)
    print("reclustered hard jet: ", rec_hard)
    rec_all = recluster(jets)
    print("reclustered jets: ", rec_all)
    print("reclustered hard jet: ", rec_all[0])
    print()

    mdt = MassDropTagger(1.0, 0.09)
    print("Using MDT: ", mdt)
    mdt_hard = mdt.result(rec_hard)
    print("MDT'd hard jet (after recluster): ", mdt_hard)
    mdt_all = mdt(rec_all)
    print("MDT'd jets (after recluster): ", mdt_all)
    print("MDT'd hard jet (after recluster): ", mdt_all[0])
    mdt_hard_norec = mdt.result(hard_jet)
    print("MDT'd hard jet (no recluster --- warning expected): ", mdt_hard_norec)
    print()

#----------------------------------------------------------------------
def read_event(filename):
    f = open(filename, 'r')
    event = []
    for line in f:
        p = line.split()
        event.append(PseudoJet(float(p[0]),float(p[1]),float(p[2]),float(p[3])));

    return event
    
if __name__ == '__main__':
    main()

