import logging
from logging import FileHandler, StreamHandler, LogRecord
from pathlib import Path

from datetime import datetime

class PeriodicFileHandler(FileHandler):
    
    def __init__(
        self, 
        pattern: str, 
        mode: str='a', 
        encoding: str = None, 
        delay: bool = False
    ) -> None:
       
        self.pattern = pattern
        filename = datetime.now().strftime(self.pattern)

        super().__init__(filename, mode, encoding, delay)

    def emit(
        self, 
        record: LogRecord
    ) -> None:

        date = datetime.fromtimestamp(record.created)
        new_filename = datetime.now().strftime(self.pattern)

        if self.stream is None:
            self.baseFileName = new_filename
        elif self.baseFilename != new_filename:
            self.close()
            self.baseFileName = new_filename

        super().emit(record)

def daily_logger(
    basename: str, 
    path: Path = Path(), 
    ext: str = 'log'
):

    """Change logfile daily"""
 
    formatter = logging.Formatter(
        "%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s"
    )
    
    path.mkdir(parents=True, exist_ok=True)
    pattern = f"{path}/{basename}-%Y%m%d.log"
    filehandler = PeriodicFileHandler(pattern)
    filehandler.setFormatter(formatter)

    consolehandler = StreamHandler()
    consolehandler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.addHandler(filehandler)
    logger.addHandler(consolehandler)

    return logger
