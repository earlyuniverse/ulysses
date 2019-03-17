# # Solver for Boltzmann equation with one decaying sterile, three flavours and CI parametrisation

################################################################################################################## 
# Input the parameters delta,a,b,theta12,theta23,theta13,x1,y1,x2,y2,x3,y3,ordering,m1,M1,M2,M3 into etaB where: #
# 														 #
# PMNS = (delta, a, b - PMNS phases, theta12, theta23, theta13 - mixing angles)					 #
# 														 #
# xys = x1,y1,x2,y2,x3,y3 - angles in orthogonal R matrix							 #
# 														 #
# lightMs = ordering - 1 normal ordering, 0 inverted ordering, m1 - free mass parameter in the light mass matrix #
# 														 #
# heavyMs = M1, M2, M3 - heavy mass parameters (masses are 10^Mi)						 #
# 														 #
# Output of etaB is the asymmetry and True or False for the Yukawa perturbativity test.				 #
##################################################################################################################

import cmath

def readConfig(fname):
    from collections import OrderedDict

    fixed  = OrderedDict()
    ranges = OrderedDict()
    with open(fname) as f:
        for line in f:
            l=line.strip()
            if l.startswith("#"):
                continue
            if len(l) == 0:
                continue

            fields = l.split()
            if len(fields)==1:
                print( "Warning, not understood instruction:", l)
                continue

            elif len(fields)==2:
                try:
                    fixed[fields[0]] = float(fields[1])
                except:
                    fixed[fields[0]] = complex(fields[1])

            elif len(fields)==3:
                ranges[fields[0]] = (float(fields[1]), float(fields[2]))

            else:
                print ("Warning, not understood instruction:", l)
                continue
        return ranges, fixed
