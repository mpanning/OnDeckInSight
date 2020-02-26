"""
Use an instaseis database to generate an approximate attenuation relationship
for peak ground displacement, velocity, and acceleration as function of 
distance from source in a 1D model
"""

import math
import numpy as np
import obspy
from tqdm import tqdm
import instaseis
# import matplotlib.pyplot as plt
# import gutenbergrichter as gr
import csv
from scipy import interpolate
import sys
python3 = sys.version_info > (3,0)

if python3:
    import pickle
else:
    import cPickle as pickle
    from requests.exceptions import ConnectionError

def wtcoef(t,t1,t2,t3,t4):
    """
    Function to calculate cosine taper

    returns weight coefficient between 0 and 1

    cosine taper from 0 to 1 t1 < t < t2
    1 for t2 < t < t3
    cosine taper from 1 to 0 t3 < t < t4
    0 for t < t1 or t > t2
    """

    if t3 > t4:
        raise ValueError('wtcoef: t3>t4')
    if t1 > t2:
        raise ValueError('wtcoef: t1>t2')

    if (t >= t2) and (t <= t3):
        wt = 1.0
    elif (t >= t4) or (t <= t1):
        wt = 0.0
    elif (t > t3) and (t < t4):
        wt = 0.5 * (1.0 + math.cos(math.pi * (t - t3)/(t4 - t3)))
    elif (t > t1) and (t < t2):
        wt = 0.5 * (1.0 + math.cos(math.pi * (t - t2)/(t2 - t1)))
    else:
        print(t, t1, t2, t3, t4)
        raise ValueError('wtcoef: this should be impossible')
    return wt

class atten_spline:
    def __init__(self, m0, pgvv, pgvh, pgvm, pgav, pgah, pgam):
        self.m0 = m0
        self.pgvv = pgvv
        self.pgvh = pgvh
        self.pgvm = pgvm
        self.pgav = pgav
        self.pgah = pgah
        self.pgam = pgam


# Set planet to 'mars' or 'europa'
planet = 'mars'
# planet = 'europa'

if (planet == 'mars'):
    # db_short = 'DWAK'
    # db_short = 'EH45Tcold'
    db_short = 'EH45TcoldCrust1b'
    # db_short = 'MAAK'
    instaseisDB = "http://instaseis.ethz.ch/blindtest_1s/{}_1s/".format(db_short)
    depth_in_km = 10.0
elif (planet == 'europa'):
    db_short = 'ice20'
    instaseisDB = 'http://instaseis.ethz.ch/icy_ocean_worlds/Eur020km-00pMS-hQ_hyd30km_2s'
    depth_in_km = 1.0
else:
    raise ValueError("Set a planet: currently should be mars or europa")
    
maxRetry = 25
db = instaseis.open_db(instaseisDB)
ifFilter = True
f1 = 0.1
f2 = 1.0
# depth_in_km = 10.0

taperFrac = 0.05 #end taper length as fraction of db record length
endCutFrac = 0.0 #Allows cutting of end of records to remove numerical probs
# endCutFrac = 0.5 #Allows cutting of end of records to remove numerical probs

dt = db.info['dt']
dbnpts = db.info['npts']

# Create taper windowing function. Can be done once, since all
# seismograms should have the same length
# Taper the end of the data to avoid abrupt endings and cut off some numerical
# problems 
t1 = 0
t2 = 0
t4 = int(dbnpts * (1 - endCutFrac)) - 1
t3 = int(t4 - taperFrac * db.info.npts)

wt = np.zeros(dbnpts)
for s in range(dbnpts):
    wt[s] = wtcoef(s, t1, t2, t3, t4)

# Define three sources [dipslip, strikeslip, oblique] all at pole
src_strikes = [0.0, 0.0, 0.0]
src_dips = [45.0, 90.0, 60.0]
src_rakes = [90.0, 0.0, 60.0]
src_lat = 90.0
src_lon = 0.0
src_m0 = 1.0e13
src_depth = 1000. * depth_in_km

distStep = 5.0
azStep = 30.0
distMax = 180.0
# distStep = 30.0
# azStep = 30.0

sources = []
for i in range(len(src_strikes)):
    sources.append(instaseis.Source.from_strike_dip_rake(latitude=src_lat,
                                                         longitude=src_lon,
                                                         depth_in_m=src_depth,
                                                         strike=src_strikes[i],
                                                         rake=src_rakes[i],
                                                         dip=src_dips[i],
                                                         M0=src_m0))

#Now loop on receiver distance in degrees
dists = []
pgv_vert_dist = []
pgv_horiz_dist = []
pgv_max_dist = []
pga_vert_dist = []
pga_horiz_dist = []
pga_max_dist = []

print('Looping on distances')
for distance in tqdm(np.arange(0.5*distStep, distMax, distStep)):
    dists.append(distance)
    rec_lat = 90.0 - distance

    pga_vert_sum = 0.
    pga_horiz_sum = 0.
    pgv_vert_sum = 0.
    pgv_horiz_sum = 0.
    pga_max = 0.
    pgv_max = 0.
    azimuths = np.arange(-180, 180, azStep)
    ntraces = len(azimuths)*len(sources)
    for azimuth in azimuths:
        rec_lon = azimuth
        receiver = instaseis.Receiver(latitude=rec_lat, longitude=rec_lon,
                                      network='XX', station='EURP')
        for source in sources:
            for i in range(maxRetry):
                try:
                    st = db.get_seismograms(source=source, receiver=receiver,
                                            remove_source_shift=False,
                                            kind='velocity')
                except ConnectionError:
                    continue
                break
            else:
                print("Could not connect after max retries")
                continue
                
            # Filter if desired
            if(ifFilter):
                st.filter('bandpass', freqmin=f1, freqmax=f2) 
            # Taper and cut
            for tr in st:
                tr.data = wt * tr.data
            
            # st.plot()
            pgv_vert = np.amax(st[0].data)
            if (pgv_vert > pgv_max):
                pgv_max = pgv_vert
            pgv_vert_sum = pgv_vert_sum + pgv_vert
            pgv_horiz = max(np.amax(st[1].data), np.amax(st[2].data))
            if (pgv_horiz > pgv_max):
                pgv_max = pgv_horiz
            pgv_horiz_sum = pgv_horiz_sum + pgv_horiz
            st.detrend('demean')
            st.detrend()
            st.differentiate()
            pga_vert = np.amax(st[0].data)
            if (pga_vert > pga_max):
                pga_max = pga_vert
            pga_vert_sum = pga_vert_sum + pga_vert
            pga_horiz = max(np.amax(st[1].data), np.amax(st[2].data))
            if (pga_horiz > pga_max):
                pga_max = pga_horiz
            pga_horiz_sum = pga_horiz_sum + pga_horiz
           
    pgv_vert_dist.append(pgv_vert_sum/ntraces)
    pgv_horiz_dist.append(pgv_horiz_sum/ntraces)
    pgv_max_dist.append(pgv_max)
    pga_vert_dist.append(pga_vert_sum/ntraces)
    pga_horiz_dist.append(pga_horiz_sum/ntraces)
    pga_max_dist.append(pga_max)

# print(pgv_vert_dist)
# print(pgv_horiz_dist)
# print(pgv_max_dist)
# print(pga_vert_dist)
# print(pga_horiz_dist)
# print(pga_max_dist)

# Write out data points to a csv file
if(ifFilter):
    filename = 'atten_curves/{}_{:.1f}km_filt_{:.3f}_{:.3f}_atten_curve.csv'.format(db_short, depth_in_km, f1, f2)
else:
    filename = 'atten_curves/{}_{:1f}km_atten_curve.csv'.format(db_short,
                                                                depth_in_km)
print(filename)
with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(['Distance', 'PGV vertical (average)',
                        'PGV horizontal (average)', 'PGV max',
                        'PGA vertical (average)', 'PGA horizontal (average)',
                        'PGA max'])
    for i in range(len(dists)):
        csvwriter.writerow([dists[i], pgv_vert_dist[i], pgv_horiz_dist[i],
                            pgv_max_dist[i], pga_vert_dist[i],
                            pga_horiz_dist[i], pga_max_dist[i]])

# Create spline interpolants and write to pickle file
splpgvv = interpolate.splrep(dists, pgv_vert_dist)
splpgvh = interpolate.splrep(dists, pgv_horiz_dist)
splpgvm = interpolate.splrep(dists, pgv_max_dist)
splpgav = interpolate.splrep(dists, pga_vert_dist)
splpgah = interpolate.splrep(dists, pga_horiz_dist)
splpgam = interpolate.splrep(dists, pga_max_dist)

interp_atten = atten_spline(src_m0, splpgvv, splpgvh, splpgvm, splpgav,
                            splpgah, splpgam)
if(ifFilter):
    filename = 'atten_curves/{}_{:.1f}km_filt_{:.3f}_{:.3f}_atten_interp.pkl'.format(db_short, depth_in_km, f1, f2)
else:
    filename = 'atten_curves/{}_{:.1f}km_atten_interp.pkl'.format(db_short,
                                                              depth_in_km)
print(filename)
with open(filename, 'wb') as f:
    pickle.dump(interp_atten, f, -1)





                   
