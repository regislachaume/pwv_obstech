from . import date as dateutils
from pathlib import Path
from importlib import resources
from astropy.time import Time

def product_directory(
    product: str, 
    date: Time | str, 
    path: Path | str = '.', 
    site: str | None = None
):

    path = Path(path)
    if site is not None:
        path = path / site
    path = path / dateutils.format(date, 'iso') / product

    return path

