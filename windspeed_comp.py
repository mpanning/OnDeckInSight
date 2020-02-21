"""
Compare windspeed and vertical component on-deck SP measurements correcting
for Viking response

Also make wind and SP time domain comparison figure
"""
from obspy import read, read_inventory
from obspy.core import Stream, UTCDateTime
from mars_tools import insight_time
import math
import matplotlib.pyplot as plt
import scipy.fftpack
import numpy as np
from tqdm import tqdm
import csv

SPdatafile = "datafiles/SHU68.ELYS0.mseed"
windsdatafile = "datafiles/SHU_windspeed.mseed"
metadata = "datafiles/ELYS0.xml"

# Read SP data and remove instrument response and plot acceleration
st_SP = read(SPdatafile)[1:] # First trace has bad response for some reason
inv = read_inventory(metadata)

pre_filt = (0.01, 0.02, 4.5, 5.0) # Only looking above 0.01 Hz

st_SP.attach_response(inv)
st_SP.remove_response(output="ACC", pre_filt=pre_filt)

st_SP.plot(outfile='SP_acc.png', method='full')
st_zoom = st_SP.copy()
zoomst = UTCDateTime(2018,12,13,5,53)
zoomet = UTCDateTime(2018,12,13,12,59)
st_zoom.trim(starttime=zoomst, endtime=zoomet)
st_zoom.plot(method='full', outfile='SP_zoom.png')

# Reread and convert to Viking DU
st_SP = read(SPdatafile)[1:] # First trace has bad response for some reason
pre_filt = (0.01, 0.05, 4.5, 5.0) # Only looking above 0.01 Hz

st_SP.attach_response(inv)
st_SP.remove_response(output="DISP", pre_filt=pre_filt)

# Viking response
mag3hz=218000;	
w0=2.0*math.pi*4.0;
Q=.6;
QF=.707;
p0=complex(0.,2.*math.pi*2.8);
wc=w0;
GN=0.76e-3
HN=0.44e-3
alpha=2;

VIK0=1./(w0*w0+w0/Q*p0+p0*p0)/((1.+p0/wc/QF+(p0/wc)**2)**2.25/p0*(1.+math.pi*.5/p0))**alpha;
# The following snippet makes the frequency dependent magnification
        # for i=1:200
	# freqvik(i)=10^(-4+(i-1)/200*6);
	# p=complex(0.,2*pi*freqvik(i));
        # VIK(i)=mag3hz/VIK0/(w0*w0+w0/Q*p+p*p)/((1.d0+p/wc/QF+(p/wc)^2)^2.25/p*(1.+pi*.5/p))^alpha;
	# accvik(i)=abs(GN/VIK(i)*p^2);
        # depvik(i)=abs(GN/VIK(i));
        # end
	# loglog(freqvik,1./depvik)

# Plot a sample Viking response
# tr = st_SP[0]
# y = tr.data
# dt = tr.stats['delta']
# ts = np.arange(0, len(y)*dt, dt)

# yf = scipy.fftpack.fft(y)
# freqs = scipy.fftpack.fftfreq(len(y), dt)

# VIK = np.zeros_like(yf)
df = 0.001
freqs = np.arange(df, 100000*df, df)
VIK = np.zeros_like(freqs, dtype=np.complex)
for i, freq in enumerate(freqs):
    if freq == 0:
        VIK[i] = 0 # ignore DC term
    else:
        p = complex(0., 2.*math.pi*freq)
        VIK[i] = mag3hz/VIK0/GN/(w0*w0+w0/Q*p+p*p)/((1.+p/wc/QF+(p/wc)**2)**2.25/p*(1.+math.pi*0.5/p))**alpha

fig = plt.figure()
plt.loglog(freqs,np.abs(VIK))
fig.savefig("vikingresponse.png")
plt.close(fig)

# Try to do fft of data and multiply by Viking response
st_SP_IR = st_SP.copy()
for tr in st_SP:
    y = tr.data
    dt = tr.stats['delta']
    ts = np.arange(0, len(y)*dt, dt)

    yf = scipy.fftpack.fft(y)
    freqs = scipy.fftpack.fftfreq(len(y), dt)

    VIK = np.zeros_like(yf)
    for i, freq in enumerate(freqs):
        if freq == 0:
            VIK[i] = 0 # ignore DC term
        else:
            p = complex(0., 2.*math.pi*freq)
            VIK[i] = mag3hz/VIK0/GN/(w0*w0+w0/Q*p+p*p)/((1.+p/wc/QF+(p/wc)**2)**2.25/p*(1.+math.pi*0.5/p))**alpha

    # # Plot Viking response
    # plt.plot(freqs, np.abs(VIK))
    # plt.show()

    # Convert displacement to Viking DU
    vikf = yf*VIK
    yvik = scipy.fftpack.ifft(vikf)

    # # plot original
    # plt.plot(ts, y)
    # plt.show()

    # # plot viking du
    # plt.plot(ts, yvik.real)
    # plt.show()

    # Replace data
    tr.data = yvik.real

#st_SP.plot(method='full')
# Read windspeed data
st_WS = read(windsdatafile).select(channel='VWS')

(st_SP + st_WS).plot(equal_scale=False, method='full')

# Create a stream just with winds that overlap SP data
st_windoverlap = Stream()
for tr in st_SP:
    stime = tr.stats['starttime']+2. # strip out first two seconds with transient
    etime = tr.stats['endtime']
    st_temp = read(windsdatafile, starttime=stime, endtime=etime).select(channel='VWS')
    st_windoverlap += st_temp

st_windoverlap = st_windoverlap.merge().split()
(st_windoverlap+st_SP).plot(method='full', equal_scale=False)
# One more plot with original
(st_windoverlap+st_SP_IR).plot(method='full', equal_scale=False)
st_windoverlap.decimate(15) #reduce sampling to 30 seconds
(st_windoverlap+st_SP).plot(method='full', equal_scale=False)
print(st_windoverlap)


# Test comparison
# Initiate figure instance and clear csv file for output
fig = plt.figure()
csvfile = 'WS_SPviking_rms.csv'
with open(csvfile, 'w') as f: # This should empty the file and add headers
    writer = csv.writer(f)
    writer.writerow(["RMS wind speed", "RMS SP Viking output"])
    
for tr in tqdm(st_windoverlap):
    stime = tr.stats['starttime']
    etime = tr.stats['endtime']
    st_SPtemp = st_SP.copy().trim(starttime=stime, endtime=etime)
    st_WStemp = st_WS.copy().trim(starttime=stime, endtime=etime)

    dt = tr.stats['delta']
    meanSPsig = np.zeros_like(tr.data)
    meanWSsig = np.zeros_like(tr.data)

    try:
        SP_times = st_SPtemp[0].times(reftime=stime)
        WS_times = st_WStemp[0].times(reftime=stime)
    except IndexError:
        continue
    # tr_times = tr.times()

    for i, sample in enumerate(tr.data):
        meanSPsig[i] = np.sqrt(np.mean(st_SPtemp[0].data[np.logical_and(SP_times>i*dt, SP_times<(i+1)*dt)]**2))
        meanWSsig[i] = np.sqrt(np.mean(st_WStemp[0].data[np.logical_and(WS_times>i*dt, WS_times<(i+1)*dt)]**2))

    plt.loglog(meanWSsig, meanSPsig, 'ko')
    rows = np.transpose(np.array([meanWSsig, meanSPsig]))
    with open(csvfile, 'a') as f: # Append to csv file for later plotting
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)
        
# Redo the last one for label purposes
plt.loglog(meanWSsig, meanSPsig, 'ko',
           label='InSight SP converted to Viking response')
xlimits = (1, 20)
ylimits = (0.1, 120)
plt.xlim(xlimits)
plt.xlabel("Mean wind speed (m/s)")
plt.ylim(ylimits)
plt.ylabel("Mean seismic noise (Viking DU)")

# Reread all data to calculate best fit to data assuming a slope of 2 in
# log log space
WSall = []
SPall = []
with open(csvfile, 'r') as f:
    reader = csv.reader(f)
    csv_headers = next(reader) # header row
    for row in reader:
        if (not np.isnan(float(row[0]))) and (not np.isnan(float(row[1]))):
            WSall.append(float(row[0]))
            SPall.append(float(row[1]))

# With fixed slope m, least squares intercept is b = mean(y) - m*mean(x)
WSall = np.array(WSall)
SPall = np.array(SPall)
m_fixed = 2.0
b_fit = np.mean(np.log10(SPall)) - m_fixed*np.mean(np.log10(WSall))
print("Best fit intercept: {}".format(b_fit))
yISfit = (math.pow(10., (m_fixed*math.log10(xlimits[0]) + b_fit)),
          math.pow(10., (m_fixed*math.log10(xlimits[1]) + b_fit)))
plt.plot(np.array(xlimits), np.array(yISfit), 'g-',
         label='InSight fit')

# Add in empirical fit from Viking
# Line on Viking fig is log(DU) = 2(log(ws/2.5))
yVfit = (math.pow(10., 2.*math.log10(xlimits[0]/2.5)),
         math.pow(10., 2.*math.log10(xlimits[1]/2.5)))
plt.plot(np.array(xlimits), np.array(yVfit), 'r-',
         label='Viking-2 fit (Anderson et al., 1977)')
plt.legend(loc='lower right')
fig.savefig("WS_SPviking_rms.png")



# # Loop over wind traces and calculate mean SP_viking over each sample
# for tr in st_windowoverlap:
#     stime = tr.stats['starttime']
#     etim = tr.stats['endtime']
#     st_temp = st_SP.copy().trim(starttime=stime, endtime=etime)
    







