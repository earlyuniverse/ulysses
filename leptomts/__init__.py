from leptomts.tools import *

from leptomts.leptocalc import LeptoCalc
from leptomts.leptomts import LeptoOld

from leptomts.etab1ds                  import EtaB_1DS
from leptomts.etab2ds                  import EtaB_2DS
from leptomts.etab3ds                  import EtaB_3DS
from leptomts.etab1dsapprox            import EtaB_1DS_Approx
from leptomts.etab2dsapprox            import EtaB_2DS_Approx
from leptomts.etab1dszerowidth         import EtaB_1DS_ZeroWidth
from leptomts.etab2dsresonant          import EtaB_2DS_Resonant
from leptomts.etab3dsscattering        import EtaB_3DS_Scattering
from leptomts.etab3dsscatteringooetaur import EtaB_3DS_Scattering_OOEtauR
from leptomts.etab3dsblanchett         import EtaB_3DS_Blanchett

testpars = {
        'delta'  :270,
        'a'      :0,
        'b'      :0,
        'theta23':48.7,
        'theta12':33.63,
        'theta13': 8.52,
        'x1'    :45,
        'y1'    :45,
        'x2'    :45,
        'y2'    :45,
        'x3'    :45,
        'y3'    :45,
        'ordering':0,
        'm1'     :-0.60206,
        'M1'     :11,
        'M2'     :12,
        'M3'     :15
        }
