from leptomts.tools import *

from leptomts.leptocalc import LeptoCalc
from leptomts.leptomts import LeptoOld

from leptomts.etab1DME                 import EtaB_1DME
from leptomts.etab2DME                 import EtaB_2DME
from leptomts.etab3DME                 import EtaB_3DME
from leptomts.etab1BE                  import EtaB_1BE
from leptomts.etab2BE                  import EtaB_2BE
from leptomts.etab3DMEscattering       import EtaB_3DME_Scattering
from leptomts.etab3DMEscatteringRHtaur import EtaB_3DS_Scattering_RHtaur

from leptomts.etab2resonant            import EtaB_2Resonant # Buggy

testpars = {
        'delta'  :270,
        'a21'      :0,
        'a31'      :0,
        't23':48.7,
        't12':33.63,
        't13': 8.52,
        'x1'    :45,
        'y1'    :45,
        'x2'    :45,
        'y2'    :45,
        'x3'    :45,
        'y3'    :45,
        'm'     :-0.60206,
        'M1'     :11,
        'M2'     :12,
        'M3'     :15
        }
