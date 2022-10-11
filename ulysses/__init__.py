from ulysses import numba

from ulysses.tools import *

from ulysses.ulsbase import ULSBase

from ulysses.etab1DME                 import EtaB_1DME
from ulysses.etab2DME                 import EtaB_2DME
from ulysses.etab3DME                 import EtaB_3DME
from ulysses.etab1BE1F                import EtaB_1BE1F
from ulysses.etab1BE2F                import EtaB_1BE2F
from ulysses.etab1BE3F                import EtaB_1BE3F
from ulysses.etab2BE1F                import EtaB_2BE1F
from ulysses.etab2BE2F                import EtaB_2BE2F
from ulysses.etab2BE3F                import EtaB_2BE3F
from ulysses.etab3DMEsct              import EtaB_3DMEsct
from ulysses.etab1BE1Fsf              import EtaB_1BE1Fsf

from ulysses.etab2RES                 import EtaB_2RES # Buggy
from ulysses.etab2RESmix              import EtaB_2RESmix
from ulysses.etab2RESsp               import EtaB_2RESsp

#shipped with version 2
from ulysses.etab1BE1F_Case2          import EtaB_1BE1F_Case2
from ulysses.etab1BE1F_Case3          import EtaB_1BE1F_Case3
from ulysses.etab1BE1F_Case4          import EtaB_1BE1F_Case4

from ulysses.etabARS                  import EtaB_ARS
from ulysses.etabARS                  import EtaB_ARS_INTERP

from ulysses.etabPBH                  import EtaB_PBH

testpars = {
        'delta' :270,
        'a21'   :0,
        'a31'   :0,
        't23'   :48.7,
        't12'   :33.63,
        't13'   : 8.52,
        'x1'    :45,
        'y1'    :45,
        'x2'    :45,
        'y2'    :45,
        'x3'    :45,
        'y3'    :45,
        'm'     :-0.60206,
        'M1'    :11,
        'M2'    :12,
        'M3'    :15
        }

BMpoint3 = {
        'delta' :31.713030 ,
        'a21'   :130.953483,
        'a31'   :649.655874,
        't12'   :33.630000 ,
        't23'   :46.633046 ,
        't13'   :8.520000  ,
        'x1'    :-72.335979,
        'y1'    :170.549206,
        'x2'    :86.969063 ,
        'y2'    :2.223559  ,
        'x3'    :-1.862141 ,
        'y3'    :178.312158,
        'm'     :-0.942835 ,
        'M1'    :6.500000  ,
        'M2'    :7.200000  ,
        'M3'    :7.900000
        }


















