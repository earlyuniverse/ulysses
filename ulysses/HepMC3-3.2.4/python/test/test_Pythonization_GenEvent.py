from pyHepMC3TestUtils import update_path
import sys

sys.path = update_path()

from pyHepMC3TestUtils import COMPARE_ASCII_FILES
from pyHepMC3 import HepMC3 as hm
import random


def test_Pythonization_GenEvent():
    evt = hm.GenEvent()
    evt.add_particle(hm.GenParticle())
    print(evt.particles)
    # evt.Particles[0].Momentum = (1, 2, 3, 4)
    # evt.Particles[0].Pid = 5
    evt.add_vertex(hm.GenVertex())
    # evt.Vertices[0].Position = (1, 2, 3, 4)
    # assert len(evt.Particles) == 1
    # assert evt.Particles[0].Id == 1
    # assert evt.Particles[0].Pid == 5
    # assert evt.Particles[0].Momentum == (1, 2, 3, 4)
    # assert len(evt.Vertices) == 1
    # assert evt.Vertices[0].Id == -1
    # assert evt.Vertices[0].Position == (1, 2, 3, 4)
    # assert repr(evt) == "GenEvent(momentum_unit=1, length_unit=0, event_number=0, particles=[GenParticle(FourVector(1, 2, 3, 4), status=0, id=1, production_vertex=0, end_vertex=-1)], vertices=[GenVertex(FourVector(1, 2, 3, 4), status=0, id=-1, particles_in=[], particles_out=[])], run_info=None)"
    return 0


if __name__ == "__main__":
    result = 1
    try:
        result = test_Pythonization_GenEvent()
    except:
        result = 1
    sys.exit(result)
