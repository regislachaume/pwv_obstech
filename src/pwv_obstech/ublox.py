from threading import Thread
from pathlib import Path
from queue import Queue
from serial import Serial, SerialException
from serial.tools.list_ports import comports as list_serial_ports
from time import sleep
import os

from dataclasses import dataclass
# import pyudev 

from typing import BinaryIO
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation 

from pygnssutils.ubxload import UBXLoader

from pyubx2.exceptions import (
    UBXMessageError, 
    UBXParseError,
    UBXStreamError, 
    UBXTypeError
)
from pyubx2 import (
    GET, 
    ERR_LOG,
    UBX_PROTOCOL, 
    UBXMessage, 
    UBXReader
)

from . import get_resource
from .utils import rinex as rnxutils
from .utils import config 

def find_ublox_device(*, vid: int = 5446, pid: int = 0) -> Path | None:

    # try instead with platform-independant list_ports.comports()
    #
    # context: pyudev.Context = pyudev.Context()
    #
    # for dev in context.list_devices(subsystem='tty'):
    #    
    #    vendor = dev.get('ID_VENDOR_FROM_DATABASE', '')
    #
    #    if not vendor.startswith('U-Blox'):
    #        continue
    #    
    #    if model_id and model_id != dev.get('ID_USB_MODEL_ID'):
    #        continue
    #            
    #    return dev.device_node

    for port in list_serial_ports():
    
        if vid and self.vid != vid:
            continue

        if pid and self.pid != pid:
            continue

        return port.device

    return None
    
def start_mjd_time(mjd: float, per: float) -> tuple[float, float]:


    mjd, sod = mjd // 1, (mjd % 1) * 86400
    try:
        sod_start = (sod // per.sec) * per.sec
    except:
        print(type(sod), type(per), sod, per)
        print("-----")
        raise
       
    return (mjd, sod_start) 


@dataclass
class RawUBXReader:

    model_id: int 
    queue: Queue
    receiver: str | None = None 
    baudrate: int = 256000 
    timeout: float = 3
    max_msg: int = 10**12 

    def load_config(self, port) -> None:

        if not self.receiver:
            return
 
        config = get_resource(f'config/{self.receiver}.conf')
        print(f"Load {config}...")
        with open(config, "rb") as in_:
            with Serial(port, self.baudrate, timeout=self.timeout) as out:
                ubl = UBXLoader(in_, out, verbosity=0, waittime=self.timeout)
                ubl.run()
        print(f"Load {config} done")

    def keep_message(self, parsed) -> bool:
        
        identity = parsed.identity
        is_ubx_message = isinstance(parsed, UBXMessage)
        
        return is_ubx_message and identity in ['RXM-RAW', 'RXM-RAWX']
            
    def parse_stream(self, stream: BinaryIO) -> None:
    
        reader = UBXReader(
            stream,
            quitonerror=ERR_LOG, protfilter=UBX_PROTOCOL,
            validate=1, msgmode=GET, parsebitfield=1
        )
           
        n_msg = 0 
        while n_msg < self.max_msg:
            
            try:
                raw = None
                (raw, parsed) = reader.read()
            except (UBXMessageError, UBXParseError,
                    UBXStreamError, UBXTypeError) as e:
                print(f"UBX parsing issue: {type(e)} {e}")
                continue
            except (EOFError, KeyboardInterrupt, SerialException):
                raise
            except Exception as e:
                print(f"unexpected exception {type(e)} {e}")
                continue
 
            if raw is None:  # EOF or timeout
                print("No input: will stop.")
                raise EOFError

            if not self.keep_message(parsed):
                continue
  
            n_msg += 1
 
            rcv_gps_time = parsed.rcvTow + parsed.week * 604800
            print(f'Got message at {rcv_gps_time:.1f}')
            self.queue.put((rcv_gps_time, raw, parsed))
    
    def start(self) -> None:

        while(True):

            try:

                # find u-blox device

                if port := find_ublox_device(pid=self.model_id):
                    print(f"U-blox device found at {port}")
                else:
                    msg = "No U-Blox device connected... will retry soon"
                    raise SerialException(msg)
               
                # write configureation to device

                self.load_config(port)

                # open serial port for reading
                
                rate = self.baudrate
                timeout = self.timeout
                with Serial(port, rate, timeout=timeout) as in_:
                    self.parse_stream(in_)

            except SerialException as e:
                print(e)
                sleep(30)

            except EOFError:
                break 

            except Exception as e:
                print(f'unforeseen termination: {e}')
                break

        print('End of messages')
        self.queue.put(None)        


@dataclass
class RawUBXConverter:

    queue: Queue
    marker: str
    period: int = 900
    frequency: int = 30
    path: Path | str = './obsdata'
    position: EarthLocation | None = None
    receiver: str = 'UNKNOWN'
    antenna: str = 'UNKNOWN'
    observer: str = 'UNKNOWN'
    institution: str = 'UNKNOWN'
    clean: bool = False

    def convert_to_rinex(self, ubx_file: Path) -> None:

        convbin_options = "-od -os -v 3 -hm {self.marker} -ht GEODETIC"
        convbin_options += " -ho {self.observer}/{self.institution}"
        convbin_options += f" -ha UNKNOWN/{self.antenna}"
        convbin_options += f" -hr UNKNOWN/{self.receiver}"

        if self.position is not None:
            pos = [p.value for p in self.position.geocentric]
            pos = '/'.join(format(p, '.3f') for p in pos)
            convbin_options += f" -hp {pos}"

        rnx_file = str(ubx_file)[:-3] + 'rnx'

        cmd = f'convbin {convbin_options} -o {rnx_file} {ubx_file}'
        if os.system(cmd):
            print(cmd)
        elif self.clean:
            ubx_file.unlink()

    def new_file(
        self, 
        file_start_time: tuple[float, float], 
        ubx_stream: BinaryIO = None
    ) -> BinaryIO:

        print(f'new file at {file_start_time[1]} s')        
        if ubx_stream is not None:
            ubx_file = ubx_stream.name
            ubx_stream.close()
            self.convert_to_rinex(ubx_file)

        mjd = file_start_time[0] + file_start_time[1] / 86400
        t = Time(mjd, format='mjd')

        ubx_file = rnxutils.file(
            marker=self.marker,
            date=t,
            period=self.period,
            frequency=self.frequency,
            constellation='M',
            filetype='ubx',
            path=self.path
        )
        ubx_file.parent.mkdir(parents=True, exist_ok=True)
        ubx_stream = open(ubx_file, 'ab')

        return ubx_stream

    def run(self) -> None:

        # in RINEX parlance, file duration is the file period and spacing
        # between measurements, the date frequency 
        # (RINEX 3.05, Table A1, pp. 47-48)

        prev_freq_start_time = (0, 0)
        prev_period_start_time = (0, 0)
        ubx_stream = None

        while (item := self.queue.get()) is not None:
    
            (rcv_gps_time, raw, parsed) = item
            mjd = Time(rcv_gps_time, format='gps').tai.mjd - 19/86400

            freq_start_time = start_mjd_time(mjd, self.frequency)

            if freq_start_time != prev_freq_start_time:
                
                prev_freq_start_time = freq_start_time 
                period_start_time = start_mjd_time(mjd, self.period)

                if period_start_time != prev_period_start_time:
                    ubx_stream = self.new_file(period_start_time, ubx_stream)
                    prev_period_start_time = period_start_time

                print(f' * kept message at {rcv_gps_time:.1f}')
                ubx_stream.write(raw)

            self.queue.task_done()

        if item is None:
            self.queue.put(item) 
    
        print('Processing done')

def record(args: list[str] | None = None) -> None: 
        
    user_types={
        'astropy.time.Time': Time, 
        'astropy.time.TimeDelta': TimeDelta,
        'pathlib.Path': Path,
    }

    parser = config.ConfigParser(
        name='ublox_record',
        description='Read raw GPS measurements from ublox receiver every FREQUENCY and write into RINEX obs files of length PERIOD'
    )
    parser.add_options_from_config(user_types=user_types)
    args = parser.parse_args(args=args)

    queue = Queue()

    if args.geodetic:
        position = EarthLocation.from_geodetic(
                *args.position, ellipsoid=args.ellipsoid
        )
    else:
        position = EarthLocation.from_geocentric(*args.position, unit="m")

    converter = RawUBXConverter(
        queue=queue, path=args.path, 
        period=args.period, frequency=args.frequency,
        marker=args.marker, antenna=args.antenna, receiver=args.receiver,
        position=position, clean=args.clean
    )
    reader = RawUBXReader(
        model_id=args.model_id, receiver=args.receiver, queue=queue
    )
        
    converter_thread = Thread(target=converter.run)
    try: 
        converter_thread.start()
        reader.start() 
    except Exception as e:
        print(e)
    else:
        converter_thread.join()
