"""
Plot up attenuation curve
"""

from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import sys
python3 = sys.version_info > (3,0)
if python3:
    import pickle
else:
    import cPickle as pickle


class atten_spline:
    def __init__(self, m0, pgvv, pgvh, pgvm, pgav, pgah, pgam):
        self.m0 = m0
        self.pgvv = pgvv
        self.pgvh = pgvh
        self.pgvm = pgvm
        self.pgav = pgav
        self.pgah = pgah
        self.pgam = pgam

# Read in attenuation curve and sample on distarray
# model = 'DWAK'
# model = 'EH45Tcold'
# model = 'MAAK'
# model = '3d_1dmean'
depth_in_km = 10.0
f1 = 0.1
f2 = 1.0
fig = plt.figure()
lines = ['k-', 'k--']
models = ['EH45Tcold', 'EH45TcoldCrust1b']
for i, model in enumerate(models):
    print("Working on model {}: {}".format(i, model))
    attenfile = ("atten_curves/" +
                 "{}_{:.1f}km_filt_{:.3f}_{:.3f}".format(model, depth_in_km,
                                                         f1, f2) +
                 "_atten_interp.pkl")
    with open(attenfile, 'rb') as f:
        atten_interp = pickle.load(f, encoding='latin1')

    diststep = 5.0
    distarray = np.arange(diststep, 180.0 + diststep, diststep)
    amparray = np.zeros_like(distarray)
    for j, dist in enumerate(distarray):
        amparray[j] = interpolate.splev(dist, atten_interp.pgav, ext=3)

    plt.plot(distarray, amparray, lines[i], label=model)
    
plt.xlabel(r'Distance ($^{\circ}$)')
plt.ylabel('Peak acceleration ($m/s^2$)')
plt.legend()
# pngfile = "{}_{:.1f}km_filt_{:.3f}_{:.3f}_atten_interp.png".format(model,
#                                                                   depth_in_km,
#                                                                   f1, f2)
pngfile = 'atten_curve.png'
plt.savefig(pngfile)
