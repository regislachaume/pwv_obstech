import numpy as np

from typing import TypeAlias
from astropy.table import Table
from astropy.coordinates import EarthLocation
from astropy.time import Time
from scipy.integrate import trapezoid

from .meteoutils import (zenithal_hydrostatic_delay,
                        saturated_vapour_pressure,
                        absolute_humidity,
                        zwd_to_pwv_factor)

def wvp_weighted_mean_temperature(
    T: list[float], # ⁰C
    P: list[float], # hPa
    h: list[float], # %
    z: list[float], # m
) -> float: # ⁰C

    delta_z = 20 # m

    zi = np.arange(z[0], z[-1], delta_z)
    Ti = np.interp(zi, z, T)
    hi = np.interp(zi, z, h)
    Pi = np.interp(zi, z, P)

    ei = hi * saturated_vapour_pressure(Ti, P=Pi)
    Tm = (ei / (Ti + 273.15))[::-1].sum() / (ei / (Ti + 273.15)**2)[::-1].sum()

    return Tm - 273.15

def precipitable_water_vapour(
    T: list[float], # ⁰C 
    P: list[float], # hPa
    h: list[float], # %
    z: list[float], # m
) -> float: # kg/m² = mm

    delta_z = 20 # m
    zi = np.arange(z[0], z[-1], delta_z)
    Ti = np.interp(zi, z, T)
    hi = np.interp(zi, z, h)
    Pi = np.interp(zi, z, P)

    Hi = absolute_humidity(Pi, Ti, hi)
    pwv = trapezoid(Hi, dx=delta_z)

    return pwv

Profile: TypeAlias = tuple[list[float], list[float], list[float], list[float]]

class VerticalProfile:
    """VerticalProfile blends a vertical profile prediction (model) with
 measured data above the ground. """

    @classmethod
    def __init__(
        self,
        T0: float,
        P0: float,
        h0: float,
        *,
        site: EarthLocation,
        t: Time,
        profile: Profile = [[], [], [], []],
    ) -> None:
        
        z, T, P, h = np.array(profile)
        # keep profile at and above ground conditions
        z0 = site.height.to_value('m')

        above = z > z0
        z = np.hstack([z0, z[above]])
        T = np.hstack([T0, T[above]])
        P = np.hstack([P0, P[above]])
        h = np.hstack([h0, h[above]])
        
        names = ['height', 'temperature', 'pressure', 'humidity']
        tab = Table([z, T, P, h], names=names)
        tab.columns['height'].unit = 'm'
        tab.columns['temperature'].unit = '⁰C'
        tab.columns['pressure'].unit = 'hPa'
        tab.columns['humidity'].unit = '%'
        for col in tab.colnames:
            tab[col].format = '.1f'

        self._profile = tab
 
        self._site = site
        self._t = t

        if len(z) > 1: 
            self._Tm = wvp_weighted_mean_temperature(T, P=P, h=h, z=z)
            self._pwv = precipitable_water_vapour(T, P=P, h=h, z=z)
        else:
            # print('NotImplementedError Vienna mapping function')
            self._Tm = T[0] - 15 # ad hoc
            self._pwv = np.nan
            
        self._zhd = zenithal_hydrostatic_delay(
            P[0], site.lat.deg, height=self._site.height.to_value('m')
        )
        self._pwv_factor = zwd_to_pwv_factor(self._Tm)
        self._zwd = self._pwv / self._pwv_factor

    def __repr__(self):

        time = self._t.isot[0:19]
        s = self._site
        loc = f"{s.lon.deg:+.5f}{s.lat.deg:+.5f}{s.height.to_value('m'):+.0f}"
        head = f"<{type(self).__name__} at coordinate {loc} ({time} UTC)>"

        rows = repr(self._profile).split("\n")[1:]

        return head + "\n" + "\n".join(rows)

    def __str__(self):
        
        s = self._site
        loc = f"({s.lon.deg:.5f}, {s.lat.deg:.5f}, {s.height.to_value('m'):.0f})"
        head = f"# location: {loc}\n# time: {self._t.isot[0:19]}\n"
        
        rows = repr(self._profile) 

        return head + rows

    @property
    def profile(self) -> Table: return self._profile.copy()

    @property
    def latitute(self) -> float: return self._site.lat.deg
    
    @property
    def longitude(self) -> float: return self._site.lon.deg

    @property
    def height(self) -> float: return self._profile['height'][0]
    
    @property    
    def temperature(self) -> float: return self._profile['temperature'][0]

    @property
    def time(self) -> Time: return self._t

    @property    
    def pressure(self) -> float: return self._profile['pressure'][0]

    @property    
    def humidity(self) -> float: return self._profile['humidity'][0]

    @property     
    def wvp_weighted_temperature(self) -> float: return self._Tm
    
    @property     
    def Tm(self) -> float: return self._Tm
 
    @property
    def zenithal_hydrostatic_delay(self) -> float: return self._zhd
    
    @property
    def zhd(self) -> float: return self._zhd

    @property
    def zenithal_wet_delay(self) -> float: return self._zwd
    
    @property
    def zwd(self) -> float: return self._zwd

    @property
    def precipitable_water_vapour(self) -> float: return self._pwv
    
    @property
    def pwv(self) -> float: return self._pwv
    
    def pwv_to_zwd(self, pwv): return pwv / self._pwv_factor 
    def zwd_to_pwv(self, zwd): return zwd * self._pwv_factor 

