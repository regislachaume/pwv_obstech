from astropy.time import Time, TimeDelta
from astropy.table import Table
import logging

from pathlib import Path
from dataclasses import dataclass
from .dataclasses import Autocast

from time import sleep

from . import date as dateutils

from typing import Union

SECOND = TimeDelta(1, format='sec')

def _parse_frequency(frequency: TimeDelta) -> str:

    frequency = frequency.sec
    if frequency == 0:
        return ''
    if frequency <= 1 / 9950 or frequency > 99.5*86400:
        return 'OOU'
    if frequency < 0.01:
        return f"{round(0.01 / frequency):02d}C"
    if frequency < 1:
        return f"{round(1 / frequency):02d}H"
    if frequency < 55.5:
        return f"{round(frequency):02d}S"
    if frequency < 59.5 * 60:
        return f"{round(frequency / 60):02d}M"
    if frequency < 23.5 * 3600:
        return f"{round(frequency / 3600):02d}H"
        
    return f"{frequency / 86400:02d}D"

def _parse_period(period: TimeDelta) -> str:

    period = period.sec

    if period == 0:
        return ''
    if period >= 1 and period < 55.5:
        return f"{round(period):02d}S"
    if period < 59.5 * 60:
        return f"{round(period / 60):02d}M"
    if period < 23.5 * 3600:
        return f"{round(period / 3600):02d}H"
    if period < 99.5 * 86400:
        return f"{round(period / 86400):02d}D"
    if period > 365 * 86400 and period < 99.5 * 365.25 * 86400:
        return f"{round(period / (365.25 * 86400))}"
    
    return "00U"

def filedate(f: Union[str, Path]) -> Time:

    if isinstance(f, Path):
        f = f.name
    
    date = f.split('_')[2]

    return dateutils.date(date)

def filepattern(
    marker: str = '',
    date: Time = Time(0, format='mjd'),
    *,  
    source: str = '',
    period: TimeDelta = TimeDelta('0s'),
    frequency: TimeDelta=  TimeDelta('0s'),
    constellation: str = '',
    datatype: str = '',
    filetype: str = 'rnx',
    version: int = 3,
) -> str:
    return filename(
        marker, date, source=source, period=period, frequency=frequency, 
        constellation=constellation, datatype=datatype, filetype=filetype, 
        version=version
    )

def filename(
    marker: str, 
    date: Time = Time(0, format='mjd'), 
    *, 
    period: TimeDelta = TimeDelta('0s'), 
    frequency: TimeDelta = TimeDelta('0s'), 
    source: str = 'R', 
    datatype: str = 'O', 
    constellation: str = 'M', 
    filetype: str = 'rnx', 
    version: int = 3
) -> Path:

    period = _parse_period(period)
    frequency = _parse_frequency(frequency)

    rnxdate = dateutils.format(date, 'rnx') if date.mjd else None

    if version == 3:

        rnxdate = rnxdate if rnxdate else '[0-9]' * 11
        marker = marker if marker else '[A-Z]' * 4 + '[0-9]' * 2 + '[A-Z]' * 3
        period = period if period else '[0-9]' * 2 + '[YDHMU]' 
        frequency = frequency if frequency else '[0-9]' * 2 + '[DHMSZCU]'
        source = source if source else '[RSU]'
        constellation = constellation if constellation else '[GREJCISM]'
        datatype = datatype if datatype else '[ONM]'
 
        filename = f"{marker}_{source}_{rnxdate}_{period}_{frequency}_{constellation}{datatype}.{filetype}"
        return filename
    
    if version == 2:
    
        marker = marker[0:4].lower() if marker else '[a-z]' * 4
        datatype = datatype if datatype else '[omn]'

        yy = rnxdate[2:4] if rxndate else '[0-9]' * 2
        doy = rnxdate[4:7] if rnxdate else '[0-9]' * 3

        if not period or not date.mjd:
            hour = '[0a-x]'
        elif 'D' in period or 'Y' in period or period == '24H':
            hour = '0'
        else:
            hour = min(int(date.isot[11:12]), 23)
            hour = chr(ord('a') + hour)
 
        if filetype.startswith('crx'):
            datatype = 'd'
            filetype = filetype[0:4]
        elif filetype.startswith('rnx'):
            filetype = filetype[0:4] 
        
        filename = f"{marker}{doy}{hour}.{yy}o.{filetype}"

        return filename

    raise NotImplementedError('Only RINEX v2 o 3 filenames are implemented')

def file(
    marker: str, 
    date: Time, 
    *, 
    source: str = 'R', 
    period: TimeDelta, 
    frequency: TimeDelta, 
    datatype: str = 'O', 
    constellation: str = 'M', 
    filetype: str = 'rnx', 
    version: int = 3,
    path: Union[Path, str] = '.'
):

    night = date.iso[0:10]
    dir = Path(path).expanduser().absolute() / marker / night
    name = filename(marker, date, period=period,  frequency=frequency,
                constellation=constellation, datatype=datatype, 
                filetype=filetype, version=3) 

    return dir / name

def merge(files: list[Union[Path, str]]) -> str:

    for i, file in enumerate(files):

        with open(file, 'r') as stream:
            lines = stream.readlines()

        for l, line in enumerate(lines):
            if 'TIME OF FIRST OBS' in line:
                break

        time_end = lines[l+1]
        if i == 0:
            header = ''.join(lines[:l+1])
            content = ''.join(lines[l+2:])
        else:
            content += ''.join(lines[l+3:])

        merged = header + time_end + content

    return merged

@dataclass(
    frozen=False, 
    # kw_only=True # > python 3.10
)
class RinexObsScanner(Autocast):
    """
    Retrieve RINEX observation data from files and merge them to
    a given run duration.

    Fields:
        marker:         RINEX marker name (site of observation)
        path:           path where RINEX files appear
        start_date:     start date
        frequency:      rate of GNSS measurements [s]
        period:         file duration [s]
        gps_run_length: length of GPS runs to be submitted [s]
        loop:           whether to constantly look for new RINEX files

    """
    marker: str
    path: Path = "./obsdata"
    start_date: Time = '2000-01-01'
    frequency: TimeDelta = '30s'
    period: TimeDelta = '15min'
    gps_run_length: TimeDelta = '8h'
    loop: bool = False

    def __call__(self):

        logger = logging.getLogger()

        expected_gps_run_size = int(self.gps_run_length / self.period + .5)

        # it's a bit tricky here, we track the end of the gps run
        # so the start date is before the asked start date.  We'll
        # trim files earlier than start date later.

        end = self.start_date + self.period - SECOND
        nloop = 0

        while (nloop := nloop + 1):

            if nloop > 1:
                if not self.loop:
                    break
                dt = self.period / 10
                logger.info(f"Wait {dt:.0f}s before scanning again")
                sleep(dt)
 
            logger.info('Looking for new RINEX files')  

            obs = self.list_obs(start=end - self.gps_run_length)
            dates = obs['date']

            last_date = obs[-1][0]
            if last_date < end:
                # print(f"No new files scanned, try again later")
                continue

            for date, file in obs:

                # do not allow to reprocess files so we need end to
                # increase in each step. 
                if date < end:
                    continue

                end = date + self.period - SECOND # avoid num. rounding errors
                start = end - self.gps_run_length
                keep = (dates  >= max(start, self.start_date)) * (dates < end)
                gps_run = obs[keep]

                # first = gps_run['date'][0]
                if len(gps_run) < 0.5 * expected_gps_run_size:
                    continue

                try:
                    merged_data = merge(gps_run['file'])
                except Exception as e:
                    logger.error('could not merge RINEX files: {e}')
                    continue

                merged_name = filename(
                    self.marker, start + 2*SECOND, # again the rounding stuff
                    period=self.gps_run_length, frequency=self.frequency,
                )
                logger.info(f'New RINEX data batch {merged_name}')               
 
                yield (merged_name, merged_data)

    def list_obs(self, start: Time) -> Table:
        """List RINEX observations"""

        # list night directories under ./obsdata/MARKER

        path = Path(self.path) / self.marker

        date_fmt = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]'
        date_dirs = sorted(path.glob(date_fmt))

        if start is not None:
            date_dirs = [d for d in date_dirs if start.iso[0:10] <= d.name]

        if not date_dirs:
            return Table(rows=[], names=['date', 'file'])

        # find all rnx files matching period, frequency, constellation, etc.

        rnx_pattern = filepattern(
            marker=self.marker, constellation='M',
            period=self.period, frequency=self.frequency
        )
        files = sorted(date_dirs[0].glob(rnx_pattern))
        if start is not None:
            files = [f for f in files if filedate(f) >= start]

        for d in date_dirs[1:]:
            files += sorted(d.glob(rnx_pattern))


        obs = Table(
            rows=[(filedate(f), f) for f in files],
            names=['date', 'file']
        )

        return obs
