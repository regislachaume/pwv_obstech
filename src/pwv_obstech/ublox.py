from threading import Thread
from pathlib import Path
from queue import Queue
from serial import Serial, SerialException
from serial.tools.list_ports import comports as list_serial_ports
from time import sleep
import os
import logging

from dataclasses import dataclass
from .utils.dataclasses import Autocast
from .utils.logging import daily_logger
# import pyudev 

from typing import BinaryIO, Union
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
from .comm.mqtt import obstech_mqtt_client, MQTTClient

MQTT_TOPIC = '/ElSauce/Weather/GNSS'

def find_ublox_device(*, vid: int = 5446, pid: int = 0) -> Path:

    for port in list_serial_ports():
    
        if vid and port.vid != vid:
            continue

        if pid and port.pid != pid:
            continue

        return port.device
               
    msg = "No U-Blox device connected..."
    raise SerialException(msg)
    
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
class RawUBXReader(Autocast):

    queue: Queue = Queue()
    model_id: int = 0 
    receiver: str = ''
    baudrate: int = 256000 
    timeout: float = 3
    max_msg: int = 10**12 
    mqtt_client: MQTTClient = obstech_mqtt_client(topic=MQTT_TOPIC)

    def load_config(self, port) -> None:
        
        if not self.receiver:
            return

        logger = logging.getLogger()
 
        config = get_resource(f'config/{self.receiver}.conf')
        logger.info('Load UBX {config}')
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
    
        logger = logging.getLogger()

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
                logger.info(f"UBX parsing issue: {type(e)} {e}")
                continue
            except (EOFError, KeyboardInterrupt, SerialException):
                raise
            except Exception as e:
                msg = f"unexpected UBX reading error {type(e).__name__} {e}"
                logger.warning(msg)
                continue
 
            if raw is None:  # EOF or timeout
                print("No input: will stop.")
                logger.error('End of UBX messages')
                raise EOFError

            if not self.keep_message(parsed):
                continue
  
            n_msg += 1
            
            # logging and communication

            rcv_gps_time = parsed.rcvTow + parsed.week * 604800
            t = Time(rcv_gps_time, format='gps').isot

            if n_msg % 20 == 0:
                logger.debug(f'Received 20 RAW UBX messages at GPS time {t}')
            
            if n_msg % 120 == 0:
                mqtt_payload = dict(
                    date=t, 
                    event='RAW GNSS messages being received',
                    filename=None
                )
                self.mqtt_client.publish(mqtt_payload)

            self.queue.put((rcv_gps_time, raw, parsed))
    
    def start(self) -> None:

        logger = logging.getLogger()

        while(True):

            try:

                # find u-blox device

                port = find_ublox_device(pid=self.model_id)
                print(f"U-blox device found at {port}")
                logger.info(f'U-blox device found at {port}')
               
                # write configureation to device

                self.load_config(port)

                # open serial port for reading
                
                rate = self.baudrate
                timeout = self.timeout
                with Serial(port, rate, timeout=timeout) as in_:
                    self.parse_stream(in_)

            except SerialException as e:
                logger.warning(f'No serial port connection: e')
                logger.warning(f'will retry in 30 s')
                print(e)
                sleep(30)

            except EOFError:
                break 

            except Exception as e:
                print(f'unforeseen termination: {e}')
                logger.error(f'unforeseen termination: {e}')
                break

        print('End of messages')
        self.queue.put(None)        


@dataclass
class RawUBXConverter(Autocast):

    position: EarthLocation 
    marker: str
    queue: Queue = Queue()
    period: TimeDelta = TimeDelta('15min')
    frequency: TimeDelta = TimeDelta('30s')
    path: Path = Path('./obsdata')
    receiver: str = 'UNKNOWN'
    antenna: str = 'UNKNOWN'
    observer: str = 'UNKNOWN'
    institution: str = 'UNKNOWN'
    clean: bool = False
    mqtt_client: MQTTClient = obstech_mqtt_client(topic=MQTT_TOPIC)

    def convert_to_rinex(self, ubx_file: Path) -> None:

        logger = logging.getLogger()

        convbin_options = f"-od -os -v 3 -hm {self.marker} -ht GEODETIC"
        convbin_options += f" -ho {self.observer}/{self.institution}"
        convbin_options += f" -ha UNKNOWN/{self.antenna}"
        convbin_options += f" -hr UNKNOWN/{self.receiver}"

        pos = [p.value for p in self.position.geocentric]
        pos = '/'.join(format(p, '.3f') for p in pos)
        convbin_options += f" -hp {pos}"

        rnx_file = str(ubx_file)[:-3] + 'rnx'

        cmd = f'convbin {convbin_options} -o {rnx_file} {ubx_file}'
        if os.system(cmd):
            logger.info(f'New RINEX file: {rnx_file}')
            logger.info(f'Created with: {cmd}')
            mqtt_payload = dict(
                date=Time.now().isot,
                event='RINEX file generated', 
                filename=rnx_file,
            )
            self.mqtt_client.publish(mqtt_payload)
        elif self.clean:
            ubx_file.unlink()

    def new_file(
        self, 
        file_start_time: tuple[float, float], 
        ubx_stream: BinaryIO = None
    ) -> BinaryIO:

        logger = logging.getLogger()

        t = Time(file_start_time[1], format='gps').isot
        msg = f'Starting new UBX file at GPS time {file_start_time[1]} s'
        logger.info(msg)
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

        logger = logging.getLogger()

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

                t = Time(rcv_gps_time, format='gps').isot
                logger.info(f'Kept RAW UBX message at GPS time {t}')
                ubx_stream.write(raw)

            self.queue.task_done()

        if item is None:
            self.queue.put(item) 
   
        logger.info('Processing done') 
        print('Processing done')

def record(args: Union[list[str], None] = None) -> None: 
        
    user_types={
        'astropy.time.Time': Time, 
        'astropy.time.TimeDelta': TimeDelta,
        'pathlib.Path': Path,
    }

    name = 'ublox_record'
    parser = config.ConfigParser(
        name=name,
        description='Read raw GPS measurements from ublox receiver every FREQUENCY and write into RINEX obs files of length PERIOD'
    )
    parser.add_options_from_config(user_types=user_types)
    args = parser.parse_args(args=args)
    
    args.position = EarthLocation(*args.position)

    queue = Queue()

    converter = RawUBXConverter(
        queue=queue, path=args.path, 
        period=args.period, frequency=args.frequency,
        marker=args.marker, antenna=args.antenna, receiver=args.receiver,     
        position=args.position, clean=args.clean,
    )
    reader = RawUBXReader(
        model_id=args.model_id, receiver=args.receiver, queue=queue
    )
   
    # logging 
    log_path = args.path / args.marker / 'logs'

    logger = daily_logger(path=log_path, basename=name)
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s - %(message)s", 
        datefmt="%Y-%m%dT%I:%M:%S",
        encoding='utf8',
    )
    
    converter_thread = Thread(target=converter.run)
    try: 
        converter_thread.start()
        reader.start() 
    except Exception as e:
        print(e)
    else:
        converter_thread.join()
