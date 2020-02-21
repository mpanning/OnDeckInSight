# OnDeckInSight
Processing software for "On-deck seismology: Lessons from InSight for future planetary seismology"

All raw data used in the paper is available for download either at the IRIS Data Management Center, or through the [Mars SEIS Data service MSDS](https://www.seis-insight.eu/en/science/science-summary).  This data is read using the [obspy python package](https://www.obspy.org).

As an example of typical processing done in this study, Panel A of figure 2 can be reproduced with `make_ppsd_panelA_fig.py`, while Panel B of Figure 2 can be reproduced with `make_paper_ppsd_fig.py`.  Equivalent files to the data and metadata files can be obtained through the webservices links at the above MSDS site.

The specific data used in these codes including specific filenames can be downloaded at the archive preserved at [this archive](http://doi.org/FILLINDOI)

All software in this repository is licensed using GNU General Public License version 3 as described in the file `COPYING`


