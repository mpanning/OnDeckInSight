"""
Make ppsd plots for SP channels on deck and compare to cruise noise estimate
"""

from obspy import read, read_inventory
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx
import matplotlib.pyplot as plt

# datafile = "SP_EH_12-1_12-13.mseed"
# datafile = "SP_SHU_MH_12-1_12-19.mseed"
# metadata = "metadata.xml"
datafile = "combined.SP.1201-1220.merged.mseed"
metadata = "ELYS0_dl12-26.xml"
cMHdata = "CRUI3.SP.mseed"
cMHmeta = "CRUI3.xml"
cEHdata = "CRUI1-2.mseed"
cEHmeta = "CRUI1.xml"

#first get cruise ppsd info
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


# channels = ['EHU', 'EHV', 'EHW']
# channels = ['SHU', 'MHV', 'MHW']
channels = ['EHU']

st = read(datafile)
inv = read_inventory(metadata)

for chn in channels:
    tr = st.select(channel=chn)[1] #first one may have metadata problem
    ppsd = PPSD(tr.stats, metadata=inv, ppsd_length=600.0, skip_on_gaps=True,
                period_limits=(0.02, 100.0), db_bins=(-200,-50, 1.))
    st_select = st.select(channel=chn)
    ppsd.add(st_select)
    plotfile = "{}_ppsd_panelA.png".format(chn)
    ppsd.plot(show=False, show_coverage=False, max_percentage=10,
              period_lim=[0.02, 100], cmap=pqlx)
    ax = plt.gca()
    fig = plt.gcf()
    ax.plot(cEHpd, cEHpsd, linewidth=2, color='darkgreen')
    ax.plot(cMHpd[1:], cMHpsd[1:], linewidth=2, color='darkgreen',
            label='Cruise PSD')
    ax.plot(0., 0., linewidth=2, color='darkgrey', label='Earth noise model')
    ax.legend(loc=1)
    plt.savefig(plotfile)


