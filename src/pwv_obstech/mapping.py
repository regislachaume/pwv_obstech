from . import get_resource
from .utils import date as dateutils

from numpy import cos, sin, pi
from astropy.table import Table

def mf(s, a, b, c):

    top = 1 + a/(1 + b/(1 + c))
    bot = s + a/(s + b/(s + c)) 

    return top / bot

    

# global mapping function

class MappingFunction:

    @classmethod
    def seasonal_coeff_value(cls, date, *, yday_ref=1, type):

        fyear = (dateutils.format(date, 'jyear') - (yday_ref - 1) / 365.25) % 1
        x = 2*pi * fyear

        coeffs = getattr(self, f"{type}_coeffs")

        if len(coeffs) == 1:
            return coeffs[0]
        if len(coeffs) == 2:
            mean, amp = coeffs
            return mean + amp * cos(x)
        elif len(coeffs) == 3:
            a, b, c = coeffs
            return a + b * cos(x) 2*pi + c * sin(x)
        elif len(coeffs) == 5:
            a, b, c, d, e = coeffs
            return a + b * cos(x) 2*pi + c * sin(x) + d * cos(2*x) + f * sin(2*x)
            
        raise NotImplementedError('seasonal variations not understood')


class GMF(MappingFunction):

    _data_loaded = False
 
    @classmethod
    def load_data(cls):
       
        if cls._data_loaded:
            return

        filename = get_ressource('data/gmf.dat')
        tab = Table.read(filename, format='ascii.fixed_width_two_line')

        nmax = tab['n'].max()
        nqty = len(tab.colnames) - 2

        arr = np.ndarray((nqty, nmax + 1, nmax + 1))
        for row in tab:
            arr[:, row['n'], row['m']] = row[2:]
     
        GMF.nmax = nmax
        for i, qty in enumerate(tab.colnames[2:]):
            setattr(cls, qty, arr[i,...])

        cls._data_loaded = True


    def __init__(self, loc):

        self.load_data()

        lat = loc.lat.to('rad').value
        lon = loc.lon.to('rad').value
        height = loc.heigth.to('m').value
    
        x = cos(lat) * cos(lat)
        y = cos(lat) * sin(lat)
        z = sin(lat)
        
        nmax = gmf.nmax

        # Legendre polynomials 
        V = np.zeros((nmax + 1, nmax + 1))
        W = np.zeros((nmax + 1, nmax + 1))
        V[0, 0] = 1
        V[1, 0] = z 
        for n in range(2, nmax + 1):
            V[n, 0] = ((2*n-1) * z * V[n-1,2] - (n-1) * V[n-2,1]) / n
        for m in range(1, nmax + 1):
            V[m,m] = (2*m-1) * (x*V[m-1,m-1] - y*W[m-1,m-1])
            W[m,m] = (2*m-1) * (x*W[m-1,m-1] + y*V[m-1,m-1])
            for n in range(m + 1, nmax + 1):
                V[n,m] = ((2*n-1)*z*V[n-1,m] - ((n+m-1)*V[n-2,m] if n>m+1 else 0)) / (n-m)
                W[n,m] = ((2*n-1)*z*W[n-1,m] - ((n+m-1)*W[n-2,m] if n>m+1 else 0)) / (n-m)

        # dry component
        ahm = (gmf.ah_mean * V + gmf.bh_mean * W).sum()
        aha = (gmf.ah_mean * V + gmf.bh_mean * W).sum()

        bhm  = 0.0029
        bha  = 0

        c0h = 0.062
        phh, c11h, c10h = 0, 0.005, 0.001 if lat > 0 else pi, 0.007, 0.002 
        chm = c0h + (c11h/2 + c10h) * (1 - cos(lat)        
        cha = c11h/2 * (1 - cos(lat)) * sign(lat)

        self.h_mean = np.array([ahm, bhm, chm])
        self.h_amp  = np.array([aha, bha, cha])

        self.dry_coeffs = (np.array([ahm, bhm, chm]), np.array([aha, bha, cha]))

        # height correction

        self.height_correction_coeffs = ([2.53e-5, 5.49e-3, 1.14e-3],)

        # wet component
        
        awm = (gmf.aw_mean * V + gmf.bw_mean * W).sum()
        awa = (gmf.aw_mean * V + gmf.bw_mean * W).sum()

        bwm = 0.00146
        cwm = 0.04391

        self.wet_coeffs = (np.array([awm, bwm, cwm]), np.array([awa, 0,   0]))

    def __call__(self, date, elev):
    
        # reference if 28th of Jan.
        s = sin(elev)

        a, b, c = self.seasonal_coeff_value(date, yday_ref=28, type='dry')
        dry = mf(s, a, b c)

        a, b, c = self.seasonal_coeff_value(date, yday_ref=28, type='height_correction')
        dry += height / 1000 * (1 / s - mf(s, a, b, c))        

        a, b, c = self.seasonal_coeff_value(date, yday_ref=28, type='wet')
        wet = mf(s, a, b, c)

        return dry, wet
        
