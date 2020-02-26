"""
Code to calculate odds of observing an event by calculating noise psds from a
supplied mseed file compared with attenuation curves for events using an 
assumed Gutenberg-Richter relationship
"""

from obspy import read, read_inventory
from obspy.signal import PPSD
from obspy.core.stream import Stream
from obspy.core import UTCDateTime
import gutenbergrichter as gr
import numpy as np
import math
import csv
from tqdm import tqdm
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from dateutil.rrule import MINUTELY, SECONDLY
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

# Define a flat response file
paz = {'gain': 1.0,
       'poles': [],
       'zeros': [],
       'sensitivity': 1.0}

# Overall parameters
# Should set up to switch to VBB Sci mode after WTS
snr = 50.0
diststep = 5.0
distarray = np.arange(diststep, 180.0 + diststep, diststep)
f1 = 0.1
f2 = 1.0
psdperiodmin = 1.0/f2 # The period range chosen to extract psd values
psdperiodmax = 1.0/f1
psdperiodstep = 1.0
psdperiodrange = np.arange(psdperiodmin, psdperiodmax, psdperiodstep)
makeplots = False

# Define a matrix of assumed Gutenberg-Richter parameters
max_m0 = math.pow(10.0, 19.5) # Assuming a fixed corner, will calculate m0total
bvals = np.array([0.9, 0.95, 1.0, 1.05, 1.1])
avals = np.arange(2.5, 6.5, 0.25)
magmin = 1.0
magmax = 8.0
magstep = 0.5
magarray = np.arange(magmin, magmax, magstep)
m0array = gr.calc_m0(magarray)

nevents = np.zeros((len(avals), len(bvals), len(magarray)))
m0total = np.zeros((len(avals), len(bvals)))

# Fill in matrix of nevents and m0total
for i, a in enumerate(avals):
    for j, b in enumerate(bvals):
        gro = gr.GutenbergRichter(a=a, b=b, max_m0=max_m0)
        gro.calc_m0total()
        m0total[i, j] = gro.m0total
        nevents[i, j, :] = gro.get_N(magarray)
        # subtract events from next bin up to avoid double counting
        for k in range(len(magarray)): 
            if (k + 1 < len(magarray)):
                nevents[i, j, k] = nevents[i, j, k] - nevents[i, j, k + 1]
        
secday = 24.0*3600.0
secyear = secday*365.0

# Read in attenuation curve and sample on distarray
# model = 'DWAK'
# model = 'EH45Tcold'
# model = 'MAAK'
# model = '3d_1dmean'
model = 'ice20'
depth_in_km = 1.0
attenfile = "atten_curves/{}_{:.1f}km_filt_{:.3f}_{:.3f}_atten_interp.pkl".format(model, depth_in_km, f1, f2)
with open(attenfile, 'rb') as f:
    # try:
    #     atten_interp = pickle.load(f)
    # except UnicodeDecodeError:
    atten_interp = pickle.load(f, encoding='latin1')
        

amparray = np.zeros_like(distarray)
afrac = np.zeros_like(distarray)
for i, dist in enumerate(distarray):
    amparray[i] = interpolate.splev(dist, atten_interp.pgav, ext=3)
    afrac[i] = 0.5*(1.0 - math.cos(math.radians(dist))) # Fraction of surface area

amp_mag_dist = []
for m0 in m0array:
    amp_mag_dist.append(m0/atten_interp.m0 * amparray)

# Get data from mseed files and process
# datafiles = ["combined.SP.1201-1220.merged.mseed",
#              "../SP_data/surf_noWTS/ELYSE.SP.1220-0114.dl0114.mseed"]
datadir = "/Users/panning/work_local/Insight/datafiles/"
# datafiles = [datadir + "weeklies/2019-02-03-2019-02-09.ELYSE.vbbrotate.mseed",
#              datadir + "weeklies/2019-02-10-2019-02-16.ELYSE.vbbrotate.mseed",
#              datadir + "weeklies/2019-02-17-2019-02-23.ELYSE.vbbrotate.mseed"]
# metafiles = [None,
#              None,
#              None]
# channels = ['VBBZ', 'VBBZ', 'VBBZ']
# locations = ['XX', 'XX', 'XX']
# datafiles = [datadir + "ELYS0.allseispress.dl0204.mseed",
#              datadir + "ELYSE.allseispress.1220-1231.dl0204.mseed", 
#              datadir + "ELYSE.allseispress.0101-0131.dl0212.mseed",
#              datadir + "ELYSE.allseispress.0201-0202.dl0211.mseed",
#              datadir + "weeklies/2019-02-03-2019-02-09.ELYSE.vbbrotate.mseed",
#              datadir + "weeklies/2019-02-10-2019-02-16.ELYSE.vbbrotate.mseed",
#              datadir + "weeklies/2019-02-17-2019-02-23.ELYSE.vbbrotate.mseed"]
# metafiles = [datadir + "ELYS0.dl0129.response.xml",
#              datadir + "ELYSE.1220-1231.dl0129.response.xml",
#              datadir + "ELYSE.0101-0131.dl0129.response.xml",
#              datadir + "ELYSE.all.dl0205.response.xml",
#              None,
#              None,
#              None]
# channels = ['SHU', 'SHU', 'SHU', 'SHU', 'VBBZ', 'VBBZ', 'VBBZ']
# locations = ['68', '68', '68', '68', 'XX', 'XX', 'XX']
datafiles = [datadir + "ELYS0.allseispress.dl0204.mseed"]
metafiles = [datadir + "ELYS0.dl0129.response.xml"]
channels = ['SHU']
locations = ['68']

# Acceptable VBBZ channels (low and high gain)
VBBZchannels = np.array(['MLZ', 'MHZ'])
VBBZlocations = ['07', '02']
# datafiles = ["combined.SP.1201-1220.merged.mseed"]
# metafiles = ["ELYS0_dl12-26.xml"]

cumlen = 0
nev = np.zeros((0, len(avals), len(bvals), len(magarray)))
stime = []
cum_prob = np.zeros((0, len(avals), len(bvals)))
cum_prob_ij = np.zeros((len(avals), len(bvals)))
cum_nev = np.zeros((0, len(avals), len(bvals)))
cum_nev_ij = np.zeros((len(avals), len(bvals)))

for i, datafile in enumerate(datafiles):
    st = read(datafile)
    if metafiles[i] is None:
        inv = paz
    else:
        inv = read_inventory(metafiles[i])
    if channels[i] == 'VBBZ':
        sts = st.select(channel=VBBZchannels[0], location=VBBZlocations[0])
        for j, channel in enumerate(VBBZchannels[1:]):
            sts += st.select(channel=channel, location=VBBZlocations[j+1])
    else:
        sts = st.select(channel=channels[i], location=locations[i])
    # Fix to remove overlaps, but not mask the data
    sts = sts.merge()
    sts = sts.split()
    sts.sort(keys=['starttime', 'endtime', 'channel'])

    
    print(sts)
    for j, tr in enumerate(sts):
        print("Working on trace {}".format(j))
        print(tr)
        length = tr.stats['endtime'] - tr.stats['starttime']
        cumlen = cumlen + length
        nevents_tr = nevents*length/secyear
        ppsd = PPSD(tr.stats, metadata=inv, ppsd_length=200.0)
        ppsd.add(Stream(tr))
        psdmean = 0
        for period in psdperiodrange:
            psds = ppsd.extract_psd_values(period)[0]
            psdmean = psdmean + math.pow(10.0, 0.05*np.mean(psds))
        psdamp = psdmean/len(psdperiodrange)    
        threshold = psdamp*snr
        print("{} Threshold: {}".format(j,threshold))
        nev_tr = np.zeros_like(nevents)
        for k, mag in enumerate(magarray):
            idx = next((x for x, v in enumerate(amp_mag_dist[k][::-1])
                        if v>threshold), None)
            if idx is not None:
                idx = len(distarray)-idx-1
                nev_tr[:, :, k] = afrac[idx]*nevents_tr[:, :, k]

        # We run into problems if cum_nev > ~0.1, so if this happens, we
        # divide the number of events evenly over smaller chunks
        nev_tr_chk = np.zeros((len(avals), len(bvals)))
        for ii in range(len(avals)):
            for jj in range(len(bvals)):
                nev_tr_chk[ii, jj] = np.sum(nev_tr[ii, jj, :])
        nev_max = np.max(nev_tr_chk)
        nev_cut = 0.1
        if nev_max < nev_cut:
            nev = np.append(nev, [nev_tr], axis=0)
        else:
            nseg = int(nev_max/nev_cut) + 1
            fac = 1.0/nseg
            nev = np.append(nev, np.tile(nev_tr*fac, (nseg, 1, 1, 1)), axis=0)
        stime.append(tr.stats['starttime'])
        for ii in range(len(avals)):
            for jj in range(len(bvals)):
                probs_not = 1.0 - np.sum(nev[:, ii, jj, :], 1)
                cum_prob_ij[ii, jj] = 1.0 - np.prod(probs_not)
                cum_nev_ij[ii, jj] = np.sum(nev[:, ii, jj, :])
        cum_prob = np.append(cum_prob, [cum_prob_ij], axis=0)
        cum_nev = np.append(cum_nev, [cum_nev_ij], axis=0)
        
print("Total length: {} hrs of data".format(cumlen/3600.0))

# Total number of expected events (roughly equal to probability if << 1)
# nev = np.array(nev)
cumulative_nev = np.zeros((len(avals), len(bvals)))
prob_final = cum_prob[-1, :, :]
for ii in range(len(avals)):
    for jj in range(len(bvals)):
        cumulative_nev[ii, jj] = np.sum(nev[:, ii, jj, :])

print("Estimated cumulative number of events:")
print(cumulative_nev)
print("Estimated probability of seeing one event:")
print(prob_final)

# Output the data as a csv file
print("Outputting data to csv files")
filename = "{}_{}_{:.1f}km_{:.3f}_{:.3f}_nev.csv".format(model, snr,
                                                         depth_in_km, f1, f2)
with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(np.append(['a values/b values'], bvals))
    for ii in range(len(avals)):
        csvwriter.writerow(np.append([avals[ii]], cumulative_nev[ii, :]))

filename = "{}_{}_{:.1f}km_{:.3f}_{:.3f}_prob.csv".format(model, snr,
                                                          depth_in_km, f1, f2)
with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(np.append(['a values/b values'], bvals))
    for ii in range(len(avals)):
        csvwriter.writerow(np.append([avals[ii]], prob_final[ii, :]))

filename = "{}_{}_{:.1f}km_{:.3f}_{:.3f}_m0total.csv".format(model, snr,
                                                             depth_in_km,
                                                             f1, f2)
with open(filename, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(np.append(['a values/b values'], bvals))
    for ii in range(len(avals)):
        csvwriter.writerow(np.append([avals[ii]], m0total[ii, :]))

if makeplots:
    # Make some plots
    print("Plotting up the data")
    trel = [] # trel is relative time in seconds from first starttime
    for i in range(len(stime)):
        trel.append(stime[i] - stime[0])
    trel = np.array(trel)
    t = (trel/secday) + mdates.date2num(stime[0].datetime)
    for ii in range(len(avals)):
        for jj in range(len(bvals)):
            fig = plt.figure()
            plt.plot(t, cum_prob[:, ii, jj], color='b',
                     label="Cumulative probability")
            plt.plot(t, cum_nev[:, ii, jj], color='r', label="Expected events")
            ax1 = plt.gca()
            ax1.xaxis_date()
            locator = mdates.AutoDateLocator(minticks=3, maxticks=6)
            locator.intervald[MINUTELY] = [1, 2, 5, 10, 15, 30]
            locator.intervald[SECONDLY] = [1, 2, 5, 10, 15, 30]
            ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
            ax1.xaxis.set_major_locator(locator)
            ax1.set_ylim(0, 1.0)
            plt.legend()
            plt.xlabel("Time (UTC)")
            plt.ylabel("Probability or number of events")
            plt.title("Threshold defined by model {} and snr {}".format(model,
                                                                        snr))
            fig.savefig("{}_{}_{}_{}_probability.png".format(model, snr,
                                                             avals[ii],
                                                             bvals[jj]))
            plt.close(fig)

            fig = plt.figure()
            plt.bar(magarray, np.sum(nev[:, ii, jj, :], 0), 0.9*magstep)
            plt.xlabel("Magnitude Mw")
            plt.ylabel("Expected number of events")
            fig.savefig("{}_{}_{}_{}_mag_events.png".format(model, snr,
                                                            avals[ii],
                                                            bvals[jj]))
            plt.close(fig)
