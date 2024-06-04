import openmeteo_requests
from astropy.time import Time, TimeDelta
from astropy.table import Table, Row
from astropy.coordinates import EarthLocation
from collections import defaultdict
import numpy as np
import requests_cache
from retry_requests import retry

from . import meteoutils 
from .vertical import VerticalProfile
from .meteomodel import VerticalProfileClient

HOUR = TimeDelta('1hr', scale='tai')
SECOND = TimeDelta('1s', scale='tai')
DAY = TimeDelta('24hr')

class OpenMeteoVerticalProfileClient(VerticalProfileClient):

    # sampling of meteo data at pressure levels in hPa

    PRESSURE_LEVELS = [
        1000, 975, 950, 925, 900, 850, 800, 700, 600, 500, 400, 300, 
        250, 200, 150, 100, 70, 50, 30
    ]
    TEMPERATURE_PROFILE_PARAMS = [
        f'temperature_{p:.0f}hPa' for p in PRESSURE_LEVELS
    ]
    HUMIDITY_PROFILE_PARAMS = [
        f'relative_humidity_{p:.0f}hPa' for p in PRESSURE_LEVELS
    ]
    HEIGHT_PROFILE_PARAMS = [
        f'geopotential_height_{p:.0f}hPa' for p in PRESSURE_LEVELS
    ]
    VERTICAL_PROFILE_PARAMS = [
        *TEMPERATURE_PROFILE_PARAMS,
        *HUMIDITY_PROFILE_PARAMS,
        *HEIGHT_PROFILE_PARAMS, 
    ]


    def __init__(self, cache_directory='.http_cache'):

        cache_session = requests_cache.CachedSession(
            cache_directory, expire_after=3600 * 24
        )
        retry_session = retry(
            cache_session, 
            retries=5, backoff_factor=0.2
        )
        self._client = openmeteo_requests.Client(session=retry_session)
        self._last_call = Time('2000-01-01T00:00:00')
    
    def query(
        self, 
        site: EarthLocation, 
        t: Time | None = None, 
        *,
        measured_conditions: dict[str, float] = {},
    ) -> Table:
       
        # avoid requests in very quick succession
        if (dt := (self._last_call - Time.now()).sec + 10) > 0:
            sleep(dt)

        now = Time.now() 
        if t is None:
            t = now
        
        # possible rounding error and leap seconds

        start = Time(int(t.mjd * 24) / 24, format='mjd') + SECOND
        end = start + HOUR 

        # determine if archive or forecast

        if (now - start).jd > 19:
            # msg = 'no vertical profile data after 19 days'
            # print(NotImplementedError(msg))
            ...
        elif (now - start).jd < -15:
            msg = 'no prediction beyond 15 days'
            raise NotImplementedError(msg)

        url = f'https://api.open-meteo.com/v1/forecast'

        # full vertical profile
        hourly = [
            'surface_pressure', 'temperature_2m', 'relative_humidity_2m',
            *self.VERTICAL_PROFILE_PARAMS
        ]

        # query
        params = dict(
            longitude=site.lon.deg,
            latitude=site.lat.deg,
            elevation=site.height.to_value('m'),
            start_hour=start.isot[0:16],
            end_hour=end.isot[0:16],
            hourly=hourly
        )

        # single site API
        response = self._client.weather_api(url, params=params)[0]

        lon = response.Longitude()
        lat = response.Latitude()
        height = response.Elevation()
        
        data = response.Hourly()

        # transform into astropy Table

        ndata = data.VariablesLength() 
        time = np.arange(data.Time(), data.TimeEnd(), data.Interval())
        time = Time(time, format='unix')    
        time.format = 'isot'
        names = ['time', *params['hourly']]
        cols = [
            time, 
            *[data.Variables(i).ValuesAsNumpy() for i in range(ndata)]
        ]
        tab = Table(cols, names=names)

        # remame columns 
 
        tab['temperature_2m'].name = 'T'
        tab['surface_pressure'].name = 'P'
        tab['relative_humidity_2m'].name= 'h'

        # interpolate to requested time

        x = (t - start) / HOUR
        
        for name in tab.colnames:
            v0 = tab[name][0]
            v1 = tab[name][-1]
            tab[name][0] = v0 +  (v1 - v0) * x
        
        row = tab[0]

        # vertical profile 
        
        z = [row[par] for par in self.HEIGHT_PROFILE_PARAMS]
        T = [row[par] for par in self.TEMPERATURE_PROFILE_PARAMS]
        P = self.PRESSURE_LEVELS
        h = [row[par] for par in self.HUMIDITY_PROFILE_PARAMS]

        # ground conditions
        # if local conditions are known, supersed surface quantities

        z0 = height
        P0 = measured_conditions.get('P', row['P'])
        T0 = measured_conditions.get('T', row['T'])
        h0 = measured_conditions.get('h', row['h'])

        profile = VerticalProfile(
            T0, P0, h0, site=site, t=t, profile=[z, T, P, h]
        )

        return profile
 
