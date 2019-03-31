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

def selectLepto(model, approx=False):
    import leptomts
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
    if model=="1ds":
        return leptomts.EtaB_1DS() if not approx else leptomts.EtaB_1DS_Approx()
    elif model=="2ds":
        return leptomts.EtaB_2DS() if not approx else leptomts.EtaB_2DS_Approx()
    elif model=="3ds":
        return leptomts.EtaB_3DS() if not approx else None
    elif model=="1dszerowidth":
        return leptomts.EtaB_1DS_ZeroWidth() if not approx  else None
    elif model=="2dsresonant":
        return leptomts.EtaB_2DS_Resonant() if not approx  else None
    elif model=="3dsscattering":
        return leptomts.EtaB_3DS_Scattering() if not approx  else None
    elif model=="3dsscatteringooetaur":
        return leptomts.EtaB_3DS_Scattering_OOEtauR() if not approx  else None
    elif model=="3dsblanchett":
        return leptomts.EtaB_3DS_Blanchett() if not approx  else None
    else:
        raise Exception("Specified model '{}' unknown".format(model))
