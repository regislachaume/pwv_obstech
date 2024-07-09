from astropy.coordinates import EarthLocation
from astropy import units as u
from astropy.units import Quantity

import pygeodesy 
from pygeodesy.ellipsoidalKarney import LatLon

class GeoidLocation(EarthLocation):

    @property
    def geoid_height(self):

        geoid_filename = str(get_resource('data/geoids/egm2008-5.pgm'))
    
        pos = LatLon(self.lat.deg, self.lon.deg)
        ginterpolator = pygeodesy.GeoidKarney(geoid_filename)
        
        return self.height - ginterpolator(pos) * u.m 

