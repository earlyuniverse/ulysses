import cmath
import leptomts

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

def selectLepto(model):
    import leptomts
    from leptomts.etab1DME                 import EtaB_1DME
    from leptomts.etab2DME                 import EtaB_2DME
    from leptomts.etab3DME                 import EtaB_3DME
    from leptomts.etab1BE                  import EtaB_1BE
    from leptomts.etab2BE                  import EtaB_2BE
    from leptomts.etab2resonant            import EtaB_2Resonant
    from leptomts.etab3dsscattering        import EtaB_3DS_Scattering
    from leptomts.etab3dsscatteringooetaur import EtaB_3DS_Scattering_OOEtauR
    from leptomts.etab3dsblanchett         import EtaB_3DS_Blanchett
    if   model=="1DME":                    return leptomts.EtaB_1DME()
    elif model=="2DME":                    return leptomts.EtaB_2DME()
    elif model=="3DME":                    return leptomts.EtaB_3DME()
    elif model=="1BE":                     return leptomts.EtaB_1BE()
    elif model=="2BE":                     return leptomts.EtaB_2BE()
    elif model=="2resonant":               return leptomts.EtaB_2Resonant()
    elif model=="3dsscattering":           return leptomts.EtaB_3DS_Scattering()
    elif model=="3dsscatteringooetaur":    return leptomts.EtaB_3DS_Scattering_OOEtauR()
    elif model=="3dsblanchett":            return leptomts.EtaB_3DS_Blanchett()
    else:
        raise Exception("Specified model '{}' unknown".format(model))
