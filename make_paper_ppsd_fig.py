"""
Make ppsd plots for SP channels on deck and compare to cruise noise estimate
"""

from obspy import read, read_inventory
from obspy.signal import PPSD
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy.imaging.cm import pqlx
import matplotlib.pyplot as plt
import numpy as np

# datafile = "SP_EH_12-1_12-13.mseed"
# datafile = "SP_SHU_MH_12-1_12-19.mseed"
# metadata = "metadata.xml"
# ondeckdatafile = "datafiles/combined.SP.1201-1220.merged.mseed"
# ondeckmetadata = "datafiles/ELYS0_dl12-26.xml"
ondeckdatafile = "datafiles/ELYS0.allseispress.dl0204.mseed"
# ondeckmetadata = "datafiles/ELYS0.dl0129.response.xml"
ondeckmetadata = "datafiles/ELYS0.all.dl0226.response.xml"
cMHdata = "datafiles/CRUI3.SP.mseed"
cMHmeta = "datafiles/CRUI3.xml"
cEHdata = "datafiles/CRUI1-2.mseed"
cEHmeta = "datafiles/CRUI1.xml"

#first get cruise ppsd info
print("Working on cruise data")
stMHc = read(cMHdata)
invMHc = read_inventory(cMHmeta)
stMHc_sel = stMHc.select(channel='MHW')
trc = stMHc_sel[0]
ppsdMHc = PPSD(trc.stats, metadata=invMHc, ppsd_length=600.0,
               skip_on_gaps=True, period_limits=(0.02, 100.0),
               db_bins=(-200, -50, 1.))
ppsdMHc.add(stMHc_sel)
(cMHpd, cMHpsd) = ppsdMHc.get_mode()
stEHc = read(cEHdata)
invEHc = read_inventory(cEHmeta)
stEHc_sel = stEHc.select(channel='EHW')
trc = stEHc_sel[0]
ppsdEHc = PPSD(trc.stats, metadata=invEHc, ppsd_length=200.0,
               skip_on_gaps=True, period_limits=(0.02, 100.0),
               db_bins=(-200, -50, 1.))
ppsdEHc.add(stEHc_sel)
(cEHpd, cEHpsd) = ppsdEHc.get_mode()

# For reference, earth low and high noise models
(nlnmpd, nlnmpsd) = get_nlnm()
(nhnmpd, nhnmpsd) = get_nhnm()

# channels = ['EHU', 'EHV', 'EHW']
# channels = ['EHU']
# channels = ['SHU', 'MHV', 'MHW']

st = read(ondeckdatafile)
inv = read_inventory(ondeckmetadata)

# On deck SP data
print("Working on on-deck data")
chn = 'EHU'
tr = st.select(channel=chn)[1] #first one may have metadata problem
ppsd = PPSD(tr.stats, metadata=inv, ppsd_length=600.0, skip_on_gaps=True,
            period_limits=(0.02, 100.0), db_bins=(-200,-50, 1.))
st_select = st.select(channel=chn)
ppsd.add(st_select)
(ondeckpd, ondeckpsd) = ppsd.get_mean()
(ondeck05pd, ondeck05psd) = ppsd.get_percentile(percentile=5)
(ondeck95pd, ondeck95psd) = ppsd.get_percentile(percentile=95)

weeklydir = 'datafiles/weeklies/'
metadata = 'datafiles/ELYSE.all.dl0226.response.xml'

# On ground SP
print("Working on on-ground data")
chn = 'SHU'
loc = '68'
weeks = np.array(['2019-01-06-2019-01-12', '2019-01-13-2019-01-19',
                  '2019-01-20-2019-01-26', '2019-01-27-2019-02-02'])

# Do first file
week = weeks[0]
datafile = weeklydir + week + '.ELYSE.allseispress.mseed'
st = read(datafile)
inv = read_inventory(metadata)
tr = st.select(channel=chn, location=loc)[0]
ppsd = PPSD(tr.stats, metadata=inv, ppsd_length=600.0, skip_on_gaps=True,
            period_limits=(0.02, 100.0), db_bins=(-200,-50, 1.))
st_select = st.select(channel=chn, location=loc)
ppsd.add(st_select)

# Do remaining weeks
for week in weeks[1:]:
    datafile = weeklydir + week + '.ELYSE.allseispress.mseed'
    st = read(datafile)
    st_select = st.select(channel=chn, location=loc)
    ppsd.add(st_select)

(ongroundpd, ongroundpsd) = ppsd.get_mean()
(onground05pd, onground05psd) = ppsd.get_percentile(percentile=5)
(onground95pd, onground95psd) = ppsd.get_percentile(percentile=95)
   
# Under WTS VBB
print("Working on under WTS data")
chn = 'BZC'
loc = '58'
weeks = ['2019-02-10-2019-02-16', '2019-02-17-2019-02-23', '2019-02-24-2019-03-02', '2019-03-03-2019-03-09', '2019-03-10-2019-03-16', '2019-03-17-2019-03-23', '2019-03-24-2019-03-30', '2019-03-31-2019-04-06', '2019-04-07-2019-04-13'] 

# Do first file
week = weeks[0]
datafile = weeklydir + week + '.ELYSE.allseispress.mseed'
st = read(datafile)
inv = read_inventory(metadata)
tr = st.select(channel=chn, location=loc)[0]
ppsd = PPSD(tr.stats, metadata=inv, ppsd_length=600.0, skip_on_gaps=True,
            period_limits=(0.02, 100.0), db_bins=(-200,-50, 1.))
st_select = st.select(channel=chn, location=loc)
ppsd.add(st_select)

# Do remaining weeks
for week in weeks[1:]:
    datafile = weeklydir + week + '.ELYSE.allseispress.mseed'
    st = read(datafile)
    st_select = st.select(channel=chn, location=loc)
    ppsd.add(st_select)

(underWTSpd, underWTSpsd) = ppsd.get_mean()
(underWTS05pd, underWTS05psd) = ppsd.get_percentile(percentile=5)
(underWTS95pd, underWTS95psd) = ppsd.get_percentile(percentile=95)

print("Plotting")
fig = plt.figure()
plotfile = 'noise_evolution.png'
plt.plot(nlnmpd, nlnmpsd, color='0.5', linestyle='-', linewidth=2)
plt.plot(nhnmpd, nhnmpsd, color='0.5', linestyle='-', linewidth=2,
         label='Earth noise')
plt.plot(ondeckpd, ondeckpsd, 'r-', label='On deck')
plt.plot(ondeck05pd, ondeck05psd, 'r--')
plt.plot(ondeck95pd, ondeck95psd, 'r--')
plt.plot(ongroundpd, ongroundpsd, 'g-', label='On ground')
plt.plot(onground05pd, onground05psd, 'g--')
plt.plot(onground95pd, onground95psd, 'g--')
plt.plot(underWTSpd, underWTSpsd, 'b-', label='Under WTS')
plt.plot(underWTS05pd, underWTS05psd, 'b--')
plt.plot(underWTS95pd, underWTS95psd, 'b--')
# Add Viking points
vikingMidnightpd = np.array([1./3., 1./3.])
vikingMidnightpsd = np.array([-103., -97.])
vikingFullpd = vikingMidnightpd
vikingFullpsd = np.array([-85., -79.])
plt.plot(vikingMidnightpd, vikingMidnightpsd, 'm+', label='Viking midnight')
plt.plot(vikingFullpd, vikingFullpsd, 'c*', label='Viking full sol')
plt.xlim((0.02, 100))
plt.xscale('log')
plt.xlabel('Period (s)')
plt.ylim((-200, -50))
plt.ylabel(r'Acceleration power (m$^2$/s$^4$/Hz [dB])')
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(plotfile)


    # plotfile = "{}_lpsc_ppsd.png".format(chn)
    # ppsd.plot(show=False, show_coverage=False, max_percentage=10,
    #           period_lim=[0.02, 100], cmap=pqlx)
    # ax = plt.gca()
    # fig = plt.gcf()
    # ax.plot(cEHpd, cEHpsd, linewidth=2, color='darkgreen')
    # ax.plot(cMHpd[1:], cMHpsd[1:], linewidth=2, color='darkgreen',
    #         label='Cruise PSD')
    # ax.plot(0., 0., linewidth=2, color='darkgrey', label='Earth noise model')
    # ax.legend(loc=1)
    # plt.savefig(plotfile)




