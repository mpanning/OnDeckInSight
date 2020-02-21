# OnDeckInSight
Processing software for "On-deck seismology: Lessons from InSight for future planetary seismology"

All raw data used in the paper is available for download either at the IRIS Data Management Center, or through the [Mars SEIS Data service MSDS](https://www.seis-insight.eu/en/science/science-summary).  This data is read using the [obspy python package](https://www.obspy.org).

The authoritative source for the archived InSight data, including all instruments is the Planetary Data System.  SEIS data is archived at the [Geosciences node](https://pds-geosciences.wustl.edu/missions/insight/index.htm).  Wind data from the TWINS sensor is archived at the [Atmospheres node](https://atmos.nmsu.edu/data_and_services/atmospheres_data/INSIGHT/insight.html).

THe following python codes are included for producing figures in the paper:

`make_ppsd_panelA_fig.py` - Panel A of figure 2

`make_paper_ppsd_fig.py` - Panel B of figure 2

`windspeed_comp.py` - Figure 3 and panel C of figure 1 (Note that the summary RMS data produced by this script is saved in the file WS_SPviking_rms.csv)

`plot_atten_curve.py` - Figure 5 (Note that this uses interpolated attenuation curves from the `atten_curves` subdirectory)

The specific data used in these codes including specific filenames can be downloaded at the archive preserved at [this archive](http://doi.org/FILLINDOI).  This data should be saved in a subdirectory called `datafiles` in order to use the python codes as written.

All software in this repository is licensed using GNU General Public License version 3 as described in the file `COPYING`


