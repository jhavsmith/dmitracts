
import sys
import numpy as np
import math
exec(open("tracts_mod.py").read())

# Read in the parameters
rate = sys.argv[1]
tstart = int(sys.argv[2])
npts = int(sys.argv[3])
maxlen = int(sys.argv[4])
pop = int(sys.argv[5])
Ls = [1]
thefile = open("psivec.txt","r")
psivec = list(np.loadtxt("psivec.txt"))

# Run tracts code
bins = np.arange(0,maxlen*(1+.5/npts),float(maxlen)/npts)
mig = continuous_mig_hybridzone(rate, tstart)
model = demographic_model(mig, psivec)
nDist = model.popNdist(1)
thefile = open("/home/joelsmith/Projects/dmis/code/dmitracts/nDist.txt", 'w+')
for item in nDist:
    thefile.write("%s\n" % item)
thefile.close()
