from astropy.time import Time, TimeDelta
from astropy.table import Table

from pathlib import Path
from dataclasses import dataclass
from time import sleep

from . import date as dateutils

SECOND = TimeDelta(1, format='sec')

def _parse_frequency(frequency: str | TimeDelta | float) -> str:

    if not isinstance(frequency, str):

        if isinstance(frequency, TimeDelta):
            frequency = frequency.sec

        if frequency <= 1 / 9950 or frequency > 99.5*86400:
            frequency = 'OOU'
        elif frequency < 0.01:
            frequency = f"{round(0.01 / frequency):02d}C"
        elif frequency < 1:
            frequency = f"{round(1 / frequency):02d}H"
        elif frequency < 55.5:
            frequency = f"{round(frequency):02d}S"
        elif frequency < 59.5 * 60:
            frequency = f"{round(frequency / 60):02d}M"
        elif frequency < 23.5 * 3600:
            frequency = f"{round(frequency / 3600):02d}H"
        else:
            frequency = f"{frequency / 86400:02d}D"

    return frequency

def _parse_period(period: str | TimeDelta | float) -> str:

    if not isinstance(period, str):

        if isinstance(period, TimeDelta):
            period = period.sec

        if period >= 1 and period < 55.5:
            period = f"{round(period):02d}S"
        elif period < 59.5 * 60:
            period = f"{round(period / 60):02d}M"
        elif period < 23.5 * 3600:
            period = f"{round(period / 3600):02d}H"
        elif period < 99.5 * 86400:
            period = f"{round(period / 86400):02d}D"
        elif period > 365 * 86400 and period < 99.5 * 365.25 * 86400:
            period = f"{round(period / (365.25 * 86400))}"
        else:
            period = "00U"

    return period

def filedate(f: Path | str) -> Time:

    if isinstance(f, Path):
        f = f.name
    date = f.split('_')[2]

    return dateutils.date(date)


def filepattern(
    marker: str | None = None,
    date: str | Time | None =None, *,  
    source: str | None = None,
    period: str | TimeDelta | float | None = None,
    frequency: str | TimeDelta | float | None = None,
    constellation: str | None = None,
    datatype: str | None = None,
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
    date: str | Time, *, 
    period: str | TimeDelta | float, 
    frequency: str | TimeDelta | float, 
    source: str = 'R', 
    datatype: str = 'O', 
    constellation: str = 'M', 
    filetype: str = 'rnx', 
    version: int = 3
) -> Path:

    if period:
        period = _parse_period(period)
    if frequency:
        frequency = _parse_frequency(frequency)

    rnxdate = dateutils.format(date, 'rnx') if date else None

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

        if not period or not date:
            hour = '[0a-x]'
        elif 'D' in period or 'Y' in period or period == '24H':
            hour = '0'
        else:
            hour = min(int(Time(date).isot[11:12]), 23)
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
    date: str | Time, *, 
    source: str = 'R', 
    period: str | TimeDelta | float, 
    frequency: str | TimeDelta | float, 
    datatype: str = 'O', 
    constellation: str = 'M', 
    filetype: str = 'rnx', 
    version: int = 3,
    path: Path | str = '.'
):

    night = Time(date).iso[0:10]
    dir = Path(path).expanduser().absolute() / marker / night
    name = filename(marker, date, period=period,  frequency=frequency,
                constellation=constellation, datatype=datatype, 
                filetype=filetype, version=3) 

    return dir / name

def merge(files: list[Path | str]) -> str:

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

@dataclass(frozen=True, kw_only=True)
class RinexObsScanner:
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
    path: str | Path = "./obsdata"
    start_date: str | Time = '2000-01-01'
    frequency: int = 30
    period: int = 900
    gps_run_length: int = 8 * 3600
    loop: bool = False

    def __call__(self):

        self.start_date = Time(self.start_date)

        period = TimeDelta(self.period, format='sec')
        gps_run_length = TimeDelta(self.gps_run_length, format='sec')
        expected_gps_run_size = gps_run_length / period

        # it's a bit tricky here, we track the end of the gps run
        # so the start date is before the asked start date.  We'll
        # trim files earlier than start date later.

        end = self.start_date + period - SECOND
        nloop = 0

        while (nloop := nloop + 1):
            
            if nloop > 1:
                if not self.loop:
                    break
                dt = self.period / 10
                # print(f"Wait {dt:.0f}s before scanning again")
                sleep(dt)

            obs = self.list_obs(start=end - gps_run_length)
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

                end = date + period - SECOND # to avoid num. rounding errors
                start = end - gps_run_length
                keep = (dates  >= max(start, self.start_date)) * (dates < end)
                gps_run = obs[keep]

                # first = gps_run['date'][0]
                if len(gps_run) < 0.5 * expected_gps_run_size:
                    continue

                try:
                    merged_data = merge(gps_run['file'])
                except Exception as e:
                    print('could not merge RINEX files: {e}')
                    continue

                merged_name = filename(
                    self.marker, start + 2*SECOND, # again the rounding stuff
                    period=gps_run_length, frequency=self.frequency,
                )
                
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
