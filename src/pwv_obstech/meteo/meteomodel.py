from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.table import Table
from abc import ABC, abstractmethod

class VerticalProfileClient(ABC):
    
    @abstractmethod
    def query(
        self,
        site: EarthLocation,
        t: Time | None = None,
        *,
        measured_conditions: dict[str, float] = {},
    ) -> Table:
        ...
    
