# Translated into python from testPolarization.cc
# garren@fnal.gov, Oct. 2010
# andrii.verbytskyi@mpp.mpg.gov, Nov. 2018
#
# In this example we will place the following event into HepMC "by hand"
#
#     name status pdg_id  parent Px       Py    Pz       Energy      Mass
#  1  !p+!    3   2212    0,0    0.000    0.000 7000.000 7000.000    0.938
#  2  !p+!    3   2212    0,0    0.000    0.000-7000.000 7000.000    0.938
# =========================================================================
#  3  !d!     3      1    1,1    0.750   -1.569   32.191   32.238    0.000
#  4  !u~!    3     -2    2,2   -3.047  -19.000  -54.629   57.920    0.000
#  5  !W-!    3    -24    1,2    1.517   -20.68  -20.605   85.925   80.799
#  6  !gamma! 1     22    1,2   -3.813    0.113   -1.833    4.233    0.000
#  7  !d!     1      1    5,5   -2.445   28.816    6.082   29.552    0.010
#  8  !u~!    1     -2    5,5    3.962  -49.498  -26.687   56.373    0.006

from pyHepMC3TestUtils import update_path, python_label
import sys

sys.path = update_path()

import random, math
from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm


def test_Polarization():
    xout1 = hm.WriterAscii(python_label() + "testPolarization1.dat")
    xout2 = hm.WriterAscii(python_label() + "testPolarization2.dat")
    xout4 = hm.WriterAsciiHepMC2(python_label() + "testPolarization4.out")
    xout5 = hm.WriterAscii(python_label() + "testPolarization5.out")

    # build the graph, which will look like
    #                       p7                   #
    # p1                   /                     #
    #   \v1__p3      p5---v4                     #
    #         \_v3_/       \                     #
    #         /    \        p8                   #
    #    v2__p4     \                            #
    #   /            p6                          #
    # p2                                         #
    #
    # define a flow pattern as  p1 . p3 . p6
    #                       and p2 . p4 . p5
    #

    # First create the event container, with Signal Process 20, event number 1
    #
    evt = hm.GenEvent(hm.Units.GEV, hm.Units.MM)
    evt.set_event_number(1)
    evt.add_attribute("signal_process_id", hm.IntAttribute(20))
    v1 = hm.GenVertex()
    evt.add_vertex(v1)
    p1 = hm.GenParticle(hm.FourVector(0, 0, 7000, 7000), 2212, 3)
    evt.add_particle(p1)
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("flow1", hm.IntAttribute(231))
    p1.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p1.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))

    v2 = hm.GenVertex()
    evt.add_vertex(v2)
    p2 = hm.GenParticle(hm.FourVector(0, 0, -7000, 7000), 2212, 3)
    evt.add_particle(p2)
    p2.add_attribute("flow1", hm.IntAttribute(243))
    p2.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p2.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v2.add_particle_in(p2)

    p3 = hm.GenParticle(hm.FourVector(0.750, -1.569, 32.191, 32.238), 1, 3)
    evt.add_particle(p3)
    p3.add_attribute("flow1", hm.IntAttribute(231))
    p3.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p3.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v1.add_particle_out(p3)
    p4 = hm.GenParticle(hm.FourVector(-3.047, -19.0, -54.629, 57.920), -2, 3)
    evt.add_particle(p4)
    p4.add_attribute("flow1", hm.IntAttribute(243))
    p4.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p4.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v2.add_particle_out(p4)

    v3 = hm.GenVertex()
    evt.add_vertex(v3)
    v3.add_particle_in(p3)
    v3.add_particle_in(p4)
    p6 = hm.GenParticle(hm.FourVector(-3.813, 0.113, -1.833, 4.233), 22, 1)
    evt.add_particle(p6)
    p6.add_attribute("flow1", hm.IntAttribute(231))
    p6.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p6.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v3.add_particle_out(p6)
    p5 = hm.GenParticle(hm.FourVector(1.517, -20.68, -20.605, 85.925), -24, 3)
    evt.add_particle(p5)
    p5.add_attribute("flow1", hm.IntAttribute(243))
    p5.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p5.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi * 2))
    v3.add_particle_out(p5)
    # create v4
    v4 = hm.GenVertex(hm.FourVector(0.12, -0.3, 0.05, 0.004))
    evt.add_vertex(v4)
    v4.add_particle_in(p5)
    p7 = hm.GenParticle(hm.FourVector(-2.445, 28.816, 6.082, 29.552), 1, 1)
    evt.add_particle(p7)
    v4.add_particle_out(p7)
    p8 = hm.GenParticle(hm.FourVector(3.962, -49.498, -26.687, 56.373), -2, 1)
    evt.add_particle(p8)
    v4.add_particle_out(p8)

    evt.add_attribute("signal_process_vertex", hm.IntAttribute(v3.id()))
    # the event is complete, we now print it out
    hm.Print.content(evt)
    hm.Print.listing(evt, 8)
    print(hm.version())
    print(hm.Print.line(v4, True))
    a = 0
    print(evt.particles())
    print(len(evt.particles()))
    print(evt.particles()[0])
    for ip in evt.particles():
        print(hm.Print.line(ip, True))
    xout1.write_event(evt)

    # write event in old format
    xout4.write_event(evt)
    # make a copy and write it
    xout5.write_event(hm.GenEvent(evt))
    # try changing polarization
    p2.add_attribute("theta", hm.DoubleAttribute(random.random() * math.pi))
    p2.add_attribute("phi", hm.DoubleAttribute(random.random() * math.pi))
    xout2.write_event(evt)
    xout1.close()
    xout2.close()
    xout4.close()
    xout5.close()
    # now clean-up by deleteing all objects from memory
    #
    # deleting the event deletes all contained vertices, and all particles
    # contained in those vertices
    evt.clear()

    assert (
        COMPARE_ASCII_FILES(python_label() + "testPolarization1.dat", python_label() + "testPolarization5.out") == 0
    ) and (COMPARE_ASCII_FILES(python_label() + "testPolarization1.dat", python_label() + "testPolarization2.dat") != 0)
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Polarization()
    except:
        result = 1
    sys.exit(result)
