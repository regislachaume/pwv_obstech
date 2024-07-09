from .utils import path as pathutils
from .utils import date as dateutils
from .utils import rinex as rnxutils
from .utils import config
from .utils.dataclasses import Autocast
from .comm.mqtt import obstech_mqtt_client, MQTTClient
from .utils.logging import daily_logger

from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.coordinates import EarthLocation

from multiprocessing.pool import ThreadPool
from threading import current_thread
from dataclasses import dataclass, field

import logging
from typing import Union

from random import random
import numpy as np
import time
import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder
from bs4 import BeautifulSoup
from pathlib import Path
from io import BytesIO
import os
import gzip
import json
import re

from zipfile import ZipFile

DOMAIN = 'https://webapp.csrs-scrs.nrcan-rncan.gc.ca'
POST_URL = f"{DOMAIN}/CSRS-PPP/service/submit"
RESULTS_URL = f"{DOMAIN}/CSRS-PPP/service/results"

def file(
    site: str, 
    date: Time, 
    /, *, 
    path: Path = Path('./products'), 
    product: str, 
    period: TimeDelta, 
    frequency: TimeDelta, 
    constellation: str
) -> Path:
    """
    Return a fully qualified file path for a given CSRS product
    correspondingt to a given start date and measurement characteristics
    """

    path = pathutils.product_directory('csrs', date, path=path, site=site)

    filetype = product
    filename = rnxutils.filename(
        site, date, 
        period=period, frequency=frequency, constellation=constellation, 
        datatype='O', filetype=filetype
    )

    return path / filename

def file_from_rinex(
    rnx_file: str, 
    *, 
    product: str, 
    path: Path = Path('./products')
) -> Path:
    """
    Return a fully qualified file path for a given CSRS product
    corresponding to a RINEX filename
    """

    date = rnxutils.filedate(rnx_file)
    site = rnx_file.split('_')[0]
    path = pathutils.product_directory('csrs', date, path=path, site=site)

    product_name = re.sub('.rnx(.gz)?', f'.{product}', rnx_file)

    return path / product_name

def process_epochs(epochs: list[str]) -> Time:

    yday = np.array([e[0:6] for e in epochs])
    sec = np.array([float(e[7:]) for e in epochs])
    mjd = np.empty_like(sec)

    for yd in np.unique(yday):
        ok = yd == yday
        m = Time(f"20{yd}", format='yday').mjd
        mjd[ok] = m + sec[ok] / 86400

    return Time(mjd, format='mjd')

def read_residuals_from_data(
    data: str, *, 
    path: Path = Path('./products'), 
    period: TimeDelta = TimeDelta('15min'), 
    frequency: TimeDelta = TimeDelta('30s'), 
    constellation: str = 'M'
) -> Table:

    data = json.loads(data)

    rows = []

    for epoch in data:

        date = epoch['epoch']
        doy = epoch['doy']
        sta  = epoch['sta']
        dir = epoch['dir']

        for sat in epoch['sat']:

            name = sat['prn']
            alt = float(sat['el'])
            az = float(sat['az'])
            row = [date, doy, sta, dir, name, alt, az]
            names = ['date', 'doy', 'station', 'dir', 
                            'sat_name', 'sat_alt', 'sat_az']

            for sig in sat['sig']:

                type = sig['type']
                res = float(sig['res'])
                for key in ['res', 'amb', 'ar', 'ini', 'samb']:
                    val = sig.get(key, None)
                    if val is not None:
                        val = (float if key != 'ini' else bool)(val)
                        names.append(f'{type}_{key}')
                        row.append(val)

            rows.append(row)

    return Table(rows=rows, names=names)

def read_pos_from_data(data: str) -> Table:

    lines = data.split('\n')[5:]
    for i, line in enumerate(lines):
        if line[0:4] not in ['HDR ', 'NOTE']:
            break

    tab = Table.read(lines[i:], format='ascii.csv', delimiter=' ')

    date_cols = ['YEAR-MM-DD','HR:MN:SS.SS']
    date = [f"{a} {b}" for a, b in tab[date_cols]]

    lon_cols = ['LONDD', 'LONMN', 'LONSS']
    a, b, c = [tab[c] for c in lon_cols]
    lon = np.sign(a) * (c / 3600 + b / 60) + a
    
    lat_cols = ['LATDD', 'LATMN', 'LATSS']
    a, b, c = [tab[c] for c in lat_cols]
    lat = np.sign(a) * (c / 3600 + b / 60) + a
   
    height  = tab['HGT(m)']
 
    tab.add_columns(
        [date, lat, lon, height], 
        names=['date', 'lat', 'lon', 'height'],
        indexes=[0,0,0,0]
    )
   
    tab['lon'].unit = 'deg' 
    tab['lat'].unit = 'deg' 
    tab.remove_columns([
        *date_cols, *lon_cols, *lat_cols, 'HGT(m)', 'DAYofYEAR'
    ])
    
    to_rename = [n for n in tab.colnames if n.endswith("(m)")]
    renamed = [n[:-3] for n in to_rename]
    tab.rename_columns(to_rename, renamed)
    for name in renamed:
        tab[name].unit = 'm'

    return tab

def read_zd_from_data(data: str) -> Table:
    
    lines = data.split("\n")

    # retrieve station coordinates
    for k, line in enumerate(lines):
        if line[0:14] == '+TROP/STA_COOR':
            break

    tab = Table.read(lines[k+1:k+3], format='ascii.csv', delimiter=' ')
    xyz = [tab[f'__STA_{c}_____'] for c in 'XYZ']
    loc = EarthLocation.from_geocentric(*xyz, unit='m')
    lon = loc.lon.value[0]
    lat = loc.lat.value[0]
    height = loc.height.value[0]

    # retrieve ZTD solution
    lines = lines[k+4:]
    for i, line in enumerate(lines):
        if line[0:14] == '+TROP/SOLUTION':
            break
    for j, line in enumerate(lines[i+1:], start=i+1):
        if line[0:14] == "-TROP/SOLUTION":
            break

    lines = lines[i+1:j]

    tab = Table.read(lines, format='ascii.csv', delimiter=' ')
    t = process_epochs(tab['____EPOCH___'])

    tab.remove_columns(['*SITE', '____EPOCH___'])
    tab.add_column(t, name='date', index=0)
    for name in tab.colnames[1:]:
        tab[name] /= 1000
        tab[name].unit = 'm'

    tab['TRODRY'].name = 'zhd'
    tab['TROWET'].name = 'zwd'

    ztd = tab['zhd'] + tab['zwd']
    tab.add_column(ztd, name='ztd', index=1)
    tab['ztd'].unit = 'm'

    tab.add_column(height, index=1, name='height')
    tab.add_column(lat, index=1, name='lat')
    tab.add_column(lon, index=1, name='lon')
    
    return tab

def wait_for_download_link(
    reqid: str, 
    max_time: float = 1800, 
    sleeptime: float = 120,
    sleeptime_variation: float = 0.3 # +/- 30% useful when threading
) -> str:

    logger = logging.getLogger()

    results_url = f"{RESULTS_URL}?id={reqid}"
       
    if max_time < 300 or max_time > 7200:
        msg = 'CSRS. Max time for request processing should be in [300, 7200] s'
        raise RuntimeError(msg)
    
    totaltime = 0 
    while totaltime <= max_time:

        dt = int(sleeptime * (1 + 2*sleeptime_variation*(random() - 0.5)))
        logger.info(f"sleep for {dt}s")
        time.sleep(dt)
        totaltime += dt
 
        logger.info(f'check results at {results_url}')
        r = requests.get(results_url)
        try:
            status = r.content.decode(encoding='utf-8', errors='strict')
        except UnicodeError:
            raise UnicodeError('ERROR: Problem with status! Try again!')

        html = BeautifulSoup(status, features='html5lib')

        # find a possible error message 

        tab = html.find_all('table')
        if len(tab) < 3:
            continue 

        tab = tab[2]
        rows = tab.find_all('tr')
        if not rows:
            continue 

        msg = rows[0].find_all('td')[1].text.lower()
        if 'error' in msg:
            raise RuntimeError(msg)

        for a in html.find_all('a'):
            if a.text == 'full_output.zip':
                return a['href']

    msg = f"CSRS: maximum job time exceeded ({max_time} s)"
    raise TimeoutError(msg)

def download_products(
    link: str, *, 
    products: list[str], 
    sleeptime: float = 10,
    sleeptime_variation: float = 0.3 # +/- 30%. useful when threading
) -> dict[str, str]:

    logger = logging.getLogger()

    max_tries = 30
    logger.info(f'download link: {link}')
 
    for n in range(max_tries):

        logger.info(f"atempt to download results #{n + 1} of {max_tries}")
        try:
            
            r = requests.get(link, timeout=5)
            content = r.content

            with ZipFile(BytesIO(content), 'r') as zip:

                members = zip.namelist()
                members = [m for m in members if m[-3:] in products]

                if not members:
                    raise RuntimeError('no data') 

                unzipped = {m: zip.read(m) for m in members}
                msg = 'the following products were downloaded and unzipped:'
                for m in members:
                    msg += f' {m}'
                logger.info(msg)

                return unzipped                       

        except Exception as e:
            message(e)
            dt = int(sleeptime * (1 + 2*sleeptime_variation*(random() - 0.5)))
            message(f"sleep for {dt} s")
            time.sleep(dt)
                
    msg = f"CSRS: maximum number of download attempts reached ({max_tries} s)"
    raise TimeoutError(msg)

def results_formatter(name: str, data: bytes) -> str:

        product = name[-3:]
        if isinstance(data, bytes):
            data = data.decode()

        if product == 'pos':
            data = read_pos_from_data(data)
        elif product == 'tro':
            data = read_zd_from_data(data)
        elif product == 'res':
            data = read_residuals_from_data(data)

        return data

def single_submit_attempt(
    zipped_name: str, 
    zipped_data: bytes, *, 
    username: str, 
    max_time: float = 1800,
    products: list[str]
) -> dict[str, str]:

    logger = logging.getLogger()
    zipped_stream = BytesIO(zipped_data)
     
    content = {
        'return_email': 'dummy_email',
        'cmd_process_type': 'std',
        'ppp_access': 'nobrowser_status',
        'language': 'en',
        'user_name': username,
        'process_type': 'Static',
        'sysref': 'ITRF',
        'nad83_epoch': 'CURR',
        'v_datum': 'cgvd2013',
        'rfile_upload': (zipped_name, zipped_stream, 'application/gzip'),
        'output_pdf': 'full',
    }
    mtp_data = MultipartEncoder(fields=content)
    header = {
        'User-Agent': 'python script',
        'Content-Type': mtp_data.content_type,
        'Accept': 'text/plain'
    }
    req = requests.post(POST_URL, data=mtp_data, headers=header)

    reqid = req.text

    if not reqid:

        logger.info("request ID does not exist. Resubmit {rnx_file}.")
        return {}

    if 'DOCTYPE' in reqid:
        logger.info(f"request has a weird value. Resubmit {rnx_file}.")
        return {}

    if reqid == 'ERROR [002]':
        msg =  "temporarily blocked from using CSRS-PPP.\n"
        msg += "Please contact CGS for further information."
        raise RuntimeError(msg)

    logger.info(f'request ID: {reqid}')

    residuals = 'res' in products
    products = [p for p in products if p != 'res']

    try:
        
        download_url = wait_for_download_link(reqid, max_time=max_time)
        csrs_products = download_products(download_url, products=products)

    except Exception as e:
        print(e)
        return {} 

    if residuals:

        res_url = download_url[:-17] + 'fid=001&type=res'

        try:
            res_products = download_products(res_url, products=['res'])
        except Exception as e:
            message(e)
            return {}

        csrs_products['res'] = res_products['res']

    return csrs_products

def submit_data(
    rnx_name: str, 
    rnx_data: Union[str, bytes],
    /, *, 
    username: str, 
    products: str = 'tro,pos,pdf,clk,sum,csv', 
    overwrite: bool = False, 
    path: Path = Path('./products'),
) -> list[Path]:
    """
    Submit RINEX data rnx_data with associated filename rnx_name to the
    Canadian CSRS PPP service.  

    Arguments:
        rnx_name:   file name
        rnx_data:   file content
        username:   e-mail of the registered user
        products:   list of products to download (file extension)
        overwrite:  whether to resubmit of result files exist
        path:       local directory where the results are stored

    Returned values: 
        files:      list of result files
    """ 
    logger = logging.getLogger()

    max_requests=5
    max_time=7200

    # If all product files exist return them, otherwise create the
    # directory where they shall go

    product_files = []

    products = [p.strip() for p in products.split(',')]

    for p in products:

        if p not in ['txt', 'tro', 'clk', 'sum', 'csv', 'pos', 'pdf', 'res']:
            raise RuntimeError(f"Unknown CSRS product: {p}")

        product_file = file_from_rinex(rnx_name, product=p, path=path)
        if not product_file.exists():
            product_file.parent.mkdir(parents=True, exist_ok=True)
            break

        product_files.append(product_file)

    else:

        if not overwrite:
            logger.info(f'Already been run for {rnx_name}. Skip')
            return product_files
        else:
            logger.info('Already been run for {rnx_name} but overwrite')

    # encode & zip data if necessary 

    if isinstance(rnx_data, str):
        rnx_data = rnx_data.encode()

    if rnx_name[-3:] != '.gz':
        rnx_name += '.gz'
        rnx_data = gzip.compress(rnx_data)
 
    # send requests, download and extract

    for n in range(1, 1 + max_requests):

        logger.info(f"sending request for {rnx_name} ({n} of {max_requests})")

        csrs_products = single_submit_attempt(
            rnx_name, rnx_data,
            username=username, products=products, max_time=max_time
        )

        if csrs_products:
            break

    else: 

        raise TimeoutError('CSRS: number of request attemps exceeded')

    # Rewrite some to easily processable astropy table 

    csrs_products = {
        name: results_formatter(name, data)
            for name, data in csrs_products.items()
    }

    # Write files to disk

    for name, data in csrs_products.items():

        product = os.path.splitext(name)[-1]
        file = file_from_rinex(name, product=product, path=path)

        format = 'ascii.fixed_width_two_line'
        if isinstance(data, Table):
            data.write(file, overwrite=True, format=format)
        else:
            with open(file, 'w') as out:
                out.write(data)

    return product_files
  
def submit_file(
    rinex_file: Path,
    username: str = 'regis.lachaume@gmail.com',
    path: Path = Path('./products'),
    products: str = 'tro,pos,pdf,clk,sum,csv', 
    overwrite: bool = False 
):
    """
    Submit a RINEX file to the Canadian CSRS PPP service.  

    Arguments:
        rnx_file:   file 
        username:   e-mail of the registered user
        products:   list of products to download (file extension)
        overwrite:  whether to resubmit of result files exist
        path:       local directory where the results are stored

    Returned values: 
        files:      list of result files
    """ 

    if not rnx_file.exists():
        raise FileNotFoundError(f'CSRS: inexistent file: {file}')
 
    # read from file
        
    with open(rnx_file, 'rb') as in_:
        rnx_data = in_.read()

    # submit data
    
    rnx_name = rnx_file.name
    product_files = submit_data(
        rnx_name, rnx_data, path=path, 
        username=username, products=products, overwrite=overwrite
    )

    return product_files

def submit(
    site: str, 
    date: Time, 
    /, *, 
    username: str, 
    path: Path = Path('./products'), 
    obsdata_path: Path = Path('./obsdata'), 
    period: TimeDelta = TimeDelta('1d'), 
    frequency: TimeDelta = TimeDelta('30s'), 
    constellation: str = 'M',
    overwrite: bool = False, 
    products: str = 'tro,pos,pdf,clk,sum,csv', 
) -> list[Path]:
    """
    Submit RINEX data rnx_data to the Canadian CSRS PPP service.  

    Arguments:
        site:           RINEX marker name
        date:           date of first data
        username:       e-mail of the registered user
        products:       list of products to download (file extension)
        overwrite:      whether to resubmit of result files exist
        path:           local directory where the results are stored
        obsdata_path:   local directory where data are stored
        period:         duration of a single file
        frequency:      rate at which GNSS are aquired 
        constellation:  GNSS constellation(s) used

    Returned values: 
        files:      list of result files
    """ 

    rnx_file = rnxutils.file(
        site, date, path=obsdata_path,
        period=period, frequency=frequency, datatype='O', 
        constellation=constellation
    )
    result_files = submit_file(
        rnx_file, path=path, 
        username=username, overwrite=overwrite, products=products,
    )

    return result_files

@dataclass(
    frozen=True, 
    # kw_only=True
)
class RinexObsSubmitter(Autocast):

    username: str 
    path: Path = Path('./products')
    products: str = 'tro,pdf'
    overwrite: bool = False 
    mqtt_client: MQTTClient = obstech_mqtt_client('CSRS') 

    def __call__(self, item: tuple[str, str]) -> list[Path]:
      
        logger = logging.getLogger()
 
        rinex_name, rinex_data = item
        logger.info(f'CSRS-PPP service: submit analysis for {rinex_name}') 

        product_files = submit_data(
            rinex_name, rinex_data, username=self.username, 
            path=self.path, products=self.products, overwrite=self.overwrite
        )

        msg = f'CSRS-PPP service: results available for {rinex_name}'
        logger.info(msg) 
        mqtt_payload = dict(
            event=msg, 
            data=[str(f) for f in product_files]
        )
        self.mqtt_client.publish(mqtt_payload)

        return product_files

def pipeline(args: Union[list[str], None] = None) -> None:

    user_types = {
        "astropy.time.Time": Time,
        "astropy.time.TimeDelta": TimeDelta,
        "pathlib.Path": Path
    }

    description = "Scan the disk for RINEX obs files and submit them to the precision point positioning service of the Canadian Spatial Reference System"

    name = 'csrs_pipeline'
    parser = config.ConfigParser(
        name=name, description=description
    )
    parser.add_options_from_config(user_types=user_types)
    args = parser.parse_args(args=args)

    # logging 
    log_path = args.products_path / args.marker / 'logs'
    logger = daily_logger(path=log_path, basename=name)
    logger.setLevel(logging.DEBUG)

    username = 'regis.lachaume@gmail.com'
    scanner = rnxutils.RinexObsScanner(
        marker=args.marker, start_date=args.start_date, path=args.data_path,
        period=args.period, frequency=args.frequency,
        gps_run_length=args.gps_run_length,
        loop=args.loop
    )
    submitter = RinexObsSubmitter(
        username=args.username, path=args.products_path,
        overwrite=args.overwrite
    )

    with ThreadPool(4) as pool:
        for item in pool.imap(submitter, scanner()):
            for file in item:
                logger.info(f"product file: {file}")
