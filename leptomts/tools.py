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

def selectLepto(model, **kwargs):
    avail = ["1DME", "2DME", "3DME", "1BE", "2BE", "2resonant", "3DMEsct", "3DMErhtau"]
    import leptomts
    if   model=="1DME":                    return leptomts.EtaB_1DME(**kwargs)
    elif model=="2DME":                    return leptomts.EtaB_2DME(**kwargs)
    elif model=="3DME":                    return leptomts.EtaB_3DME(**kwargs)
    elif model=="1BE":                     return leptomts.EtaB_1BE(**kwargs)
    elif model=="2BE":                     return leptomts.EtaB_2BE(**kwargs)
    elif model=="2resonant":               return leptomts.EtaB_2Resonant(**kwargs)
    elif model=="3DMEsct":                 return leptomts.EtaB_3DME_Scattering(**kwargs)
    elif model=="3DMErhtau":               return leptomts.EtaB_3DS_Scattering_RHtaur(**kwargs)
    else:
        raise Exception("Specified model '{}' unknown.\n Select from: {}".format(model, avail))
