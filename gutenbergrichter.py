"""
Module for calculating and using Gutenberg-Richter relationships

Calculations based on Golombek et al., 1992, with moment magnitude
defined by Aki and Richards converted to Newton-meters
log m0 = 1.5 Mw + 9.1

a and b are the coefficients in the classic Gutenberg-Richter relationship
log N(M) = a - bM

m0total is the total cumulative moment released in a year
max_m0 is the assumed maximum possible event size

mA and mB are internal values used in calculations based on the moment and
magnitude relationship used in Golombek et al., 1992
N = A m0^(-B)

All N values are assumed to be normalized to 1 Earth year
"""

import math
import numpy as np
import random
from tqdm import tqdm
import csv

secday = 24.0*60.0*60.0
secyear = secday*365.0
rad2deg = 180.0/math.pi

# Default values for catalogs
m0total = 1.0e17
max_m0 = math.pow(10.0, 19.5)
slope = 1.0
default_catlength = 2.0*secday


class Catalog(object):
    """
    An object with an event catalog generated from a GutenbergRichter object
    """

    def __init__(self, data=None, max_dep=None, length=None, Mws=None,
                 Ns=None, id_dict=None):
        """
        An object to hold event catalog data

        data should be an array to hold the catalog data for each event
        catlength is the duration of the catalog in seconds
        Mws and Ns are the magnitude and N values for the G-R plot
        """

        self.data = data
        self.length = length
        self.max_dep = max_dep
        self.Mws = Mws
        self.Ns = Ns
        self.id_dict = id_dict

class GutenbergRichter(object):
    """
    An object containing Gutenberg-Richter relationship information
    """

    def __init__(self, a=None, b=None, m0total=None, max_m0=None, min_m0=None,
                 max_dep=None):
        self.a = a
        self.b = b
        self.m0total = m0total
        self.max_m0 = max_m0
        self.min_m0 = min_m0
        if max_dep is None:
            self.max_dep = 10.0
        else:
            self.max_dep = max_dep

    def calc_a(self):
        """
        function to determine a from the G-R relationship

        Follows Golombek et al., 1992 to calculate a in a standard G-R 
        relationship when supplied b, the total seismic moment released,
        and the moment of the maximum size event
        """
        if (self.b is None) or (self.m0total is None) or (self.max_m0 is None):
            raise ValueError("b, m0total, and max_m0 must be set")
        self.mB = 2.0*self.b/3.0
        self.mA = (1.0 - self.mB) * self.m0total/(self.mB * 
                                                  self.max_m0**(1.0 - self.mB))
        self.a = math.log10(self.mA) - 9.1*self.mB

    def calc_m0total(self):
        """
        Function to determine the total moment release
        """
        if (self.a is None) or (self.b is None) or (self.max_m0 is None):
            raise ValueError("a, b, and max_m0 must be set")
        self.mB = 2.0*self.b/3.0
        self.mA = math.pow(10.0, (self.a + 9.1*self.mB))
        self.m0total = ((self.mA*self.mB/(1.0 - self.mB))*
                        self.max_m0**(1.0 - self.mB))

    def calc_max_m0(self):
        """
        Function to determine the maximum event size
        """
        if (self.a is None) or (self.b is None) or (self.m0total is None):
            raise ValueError("a, b, and m0total must be set")
        self.mB = 2.0*self.b/3.0
        self.mA = math.pow(10.0, (self.a + 9.1*self.mB))
        self.max_m0 = math.pow((1.0 - self.mB)/(self.mA*self.mB)*self.m0total,
                               1.0/(1.0 - self.mB))

    def get_N(self, Mw):
        """
        Function to retrieve the number of events at a magnitude level

        Mw can be a single value or an np array
        """
        if (self.a is None) or (self.b is None):
            raise ValueError("a and b must be set")
        N = np.power(10.0, (self.a - self.b*Mw))
        return N

    def generate_catalog(self, length, max_dep=None, Mws=None, Msamp=None):
        """
        Function to generate a catalog of events of a specified length in 
        seconds

        Returns a tuple of catalog,Ns,Mws

        catalog: an array with 6 columns 
         [time (s), magnitude (Mw), dist (deg), backaz (deg),
          strike, rake, dip (all in deg)]
        Ns: an array of number of events greater than or equal to a given mag
        Mws: the magnitude values for Ns
        """
        if Mws is None:
            if Msamp is None:
                Msamp = 0.25
            if self.max_m0 is None:
                maxM = 9.0
            else:
                maxM = calc_Mw(self.max_m0)
            if self.min_m0 is None:
                minM = 0.0
            else:
                minM = calc_Mw(self.min_m0)

            Mws = np.arange(minM, maxM + Msamp, Msamp)

        if max_dep is None:
            max_dep = 10.0

        Nsec = self.get_N(Mws)/secyear

        Ns = np.zeros_like(Nsec)
        catalog = []
        for sec in tqdm(range(0, int(length))):
            ran = random.random()
            mask = (Nsec > ran).astype(int)
            if (mask.max() == 1):
                Ns = Ns + mask
                index=np.amax(np.where(mask == 1)[0])
                time = sec + random.random()
                mag = Mws[index] + random.random()*Msamp
                theta = math.acos(random.uniform(-1.0,1.0))
                delta = rad2deg * theta
                backaz = random.uniform(0,360)
                depth = random.uniform(0,max_dep)
                strike = random.uniform(0,360)
                rake = random.uniform(0,360)
                dip = random.uniform(0,90)
                catalog.append(np.array([time, mag, delta, backaz, depth, 
                                         strike, rake, dip]))
                

        id_dict = {'time' : 0, 'magnitude' : 1, 'delta' : 2, 'backaz' : 3,
                   'depth' : 4, 'strike' : 5, 'rake' : 6, 'dip' : 7}
        catalog = np.array(catalog)
        self.catalog = Catalog(data=catalog, length=length, max_dep=max_dep,
                               Ns=Ns, Mws=Mws, id_dict=id_dict)

        #return (catalog,Ns,Mws)

    def write_csv(self, filename=None):
        """
        Function to output a csv file with all catalog details
        """

        if filename is None:
            filename='gr_catalog.csv'
        with open(filename, 'w') as f:
            output = csv.writer(f)
            output.writerow(['a', 'b', 'M0 total', 'Max M0', 'Min M0',
                             'Max Depth'])
            output.writerow([str(self.a), str(self.b), str(self.m0total),
                             str(self.max_m0), str(self.min_m0),
                             str(self.max_dep)])
            output.writerow(['Length (s)', self.catalog.length])
            output.writerow(['Magnitude and numbers stats'])
            output.writerow(self.catalog.Mws)
            output.writerow(self.catalog.Ns)
            output.writerow(['Time', 'Mw', 'Colatitude', 'Longitude', 'Depth',
                             'Strike', 'Rake', 'Dip'])
            output.writerows(self.catalog.data)
        
def calc_m0(Mw):
    """
    Function to calculate seismic moment (in Nm) from a moment magnitude
    """

    m0 = 10.0**(1.5*Mw + 9.1)
    return m0

def calc_Mw(m0):
    """
    Function to calculate moment magnitude from a moment given in Nm
    """

    Mw = 2.0*(math.log10(m0) - 9.1)/3.0
    return Mw
    
