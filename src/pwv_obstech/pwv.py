from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.coordinates import EarthLocation

from pathlib import Path
from dataclasses import dataclass
from time import sleep

from multiprocessing.pool import ThreadPool

import numpy as np

from .utils.dataclasses import Autocast 
from .utils import date as dateutils
from .utils.rinex import filepattern, filedate
from .meteo.meteomodel import VerticalProfileClient
from .meteo.openmeteo import OpenMeteoVerticalProfileClient
from .meteo.database import MeteoDatabase, MySQLMeteoDatabase
from .meteo.vertical import VerticalProfile

from .utils import config

SECOND = TimeDelta('1s')

@dataclass(frozen=True, kw_only=True)
class PWVModel:

    meteo_model: VerticalProfileClient | None = OpenMeteoVerticalProfileClient()
    meteo_db: MeteoDatabase | None = None
    period: TimeDelta('15min')

    def __call__(self, tab: Table):
    
        mean = dict(tab.groups.aggregate(np.mean)[0])
        mean = {(k + '_ppp' if k in ['ztd', 'zhd', 'zwd'] else k): v 
                        for k, v in mean.items()
                            if not k.startswith(('STDDEV', 'TG'))}
 
        t = Time(mean['date'], format='mjd')
        site = EarthLocation(mean['lon'], mean['lat'], mean['height'])

        profile = None
        measured = dict()
        if self.meteo_db is not None:
            T, P, h = self.meteo_db.mean_conditions(t, period)
            if all(not np.isnan(x) for x in [T, P, h]):
                measured = [T, P, h]
                profile = VerticalProfile(*measured, site=site, t=t)

        if self.meteo_model is not None:
            profile = self.meteo_model.query(
                site=site, t=t, measured_conditions=measured
            )

        if profile is None:
            return np.nan
        
        mean['ztd_weather'] = profile.zwd + profile.zhd
        mean['zhd_weather'] = profile.zhd
        mean['zwd_weather'] = profile.pwv
            
        mean['temperature'] = profile.temperature
        mean['pressure'] = profile.pressure
        mean['humidity'] = profile.humidity
        
        mean['zwd'] = mean['ztd_ppp'] - profile.zhd
        mean['pwv'] = profile.zwd_to_pwv(mean['zwd']) 
       
        return mean 

@dataclass(frozen=True, kw_only=True)
class PPPProductScanner(Autocast):
    """
    Retrieve PPP products files 

    Fields:
        marker:         RINEX marker name (site of observation)
        service:        Name of the service (e.g. csrs for CSRS-PPP)
        product:        Type of product (e.g. tro for zenithal delays)
        path:           path where RINEX files appear
        start_date:     start date
        frequency:      rate of GNSS measurements [s]
        period:         periodicity of file [s]
        gps_run_length: length of data in files [s]
        loop:           whether to constantly look for new files

    """
    marker: str
    product: str 
    service: str = 'csrs' 
    path: Path = "./products"
    start_date: Time = '2000-01-01'
    frequency: TimeDelta = '30s'
    period: TimeDelta = '15min'
    gps_run_length: TimeDelta = '8hr'
    loop: bool = False

    def __call__(self):

        file_start_date = self.start_date
        end = file_start_date

        nloop = 0
        loop_wait = self.period.sec / 100

        while (nloop := nloop + 1) <= 1 + self.loop * 10_000_000:
           
            #  print(f'loop: {nloop}') 
            obs = self.list_files(start=file_start_date - SECOND) # rounding

            for file_start_date, file in obs:

                # print(f"examine {file}")
                
                tab = Table.read(file, format='ascii.fixed_width_two_line')
                file_end_date = file_start_date + self.gps_run_length
               
                # Files have overlap, we don not want to give several
                # data for the same date and time range.

                if file_start_date >= end:
                    start = file_start_date
                else:
                    start = end
                    tab = tab[tab['date'] >= start.mjd]
                
                nperiods = self.gps_run_length / self.period
                while start < file_end_date - SECOND:

                    end = start + self.period 
                    mjd = tab['date']
                    eps = 1e-6
                    subtab = tab[(mjd >= start.mjd - eps) & (mjd < end.mjd - eps)]

                    if len(subtab):
                        yield subtab 

                    start = end
 
            if self.loop:
                # print(f'wait for {loop_wait} s before scanning again')
                sleep(loop_wait)

    def list_files(self, start: Time) -> Table:
        """List PPP products more recent than start"""

        # Filenames follow
        # ./<path>/<marker>/<YYYY-MM-DD>/<service>/MARKER/***.<product>
        #
        # 1. scan by date directory with date on or later than <start>
        # 2. scan in night directories, on the date, restrict to files
        #    later than <start>
        #
        # TODO: use watchdog

        path = self.path / self.marker 

        date_fmt = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]'
        date_dirs = sorted(path.glob(date_fmt))

        if start is not None:
            date_dirs = [d for d in date_dirs if start.iso[0:10] <= d.name]
            date_dirs = [d / self.service for d in date_dirs]

        if not date_dirs:
            return Table(rows=[], names=['date', 'file'])

        # find all files matching period, frequency, constellation, etc.

        prod_pattern = filepattern(
            marker=self.marker, constellation='M',
            period=self.gps_run_length, frequency=self.frequency,
            filetype=self.product
        )

        files = sorted(date_dirs[0].glob(prod_pattern))
        if start is not None:
            files = [f for f in files if filedate(f) >= start]

        for file in files:
            yield filedate(file), file

        for d in date_dirs[1:]:
            for file in sorted(d.glob(prod_pattern)):
                yield filedate(file), file

def pipeline():

    user_types = {
        "astropy.time.Time": Time,
        "astropy.time.TimeDelta": TimeDelta,
        "pathlib.Path": Path
    }
    description = "Scan the disk for Precise point positioning files and obtain the precipitable water vapour using meteo data"
    parser = config.ConfigParser(
        description=description, name='pwv_pipeline'
    )
    parser.add_options_from_config(user_types=user_types)
    args = parser.parse_args()

    scanner = PPPProductScanner(
        marker=args.marker, path=args.products_path,
        gps_run_length=args.gps_run_length, period=args.period,
        frequency=args.frequency, start_date=args.start_date,
        product='tro',
        loop=args.loop
    )
   
    meteo_model = OpenMeteoVerticalProfileClient()
 
    meteo_db = MySQLMeteoDatabase(
        server=args.meteodb_server,
        database=args.meteodb_database,
        user=args.meteodb_user,
        password=args.meteodb_password,
        table=args.meteodb_table,
        columns=args.meteodb_columns 
    )

    model = PWVModel(meteo_db=meteo_db, meteo_model=meteo_model)

    with ThreadPool(2) as pool:
        for r in pool.imap(model, scanner()):
            t = Time(r['date'], format='mjd')
            print(f"{t.isot[0:16]} {r['pwv'] * 1000:.1f}") 
            sleep(1)
