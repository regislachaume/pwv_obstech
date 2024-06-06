from dataclasses import dataclass, field
from mysql import connector
from astropy.time import Time, TimeDelta
import numpy as np
from abc import ABC, abstractmethod

@dataclass(frozen=True,kw_only=True)
class MeteoDatabase(ABC):

    @abstractmethod
    def mean_conditions(self, t: Time, period: TimeDelta): ...

@dataclass(frozen=True,kw_only=True)
class MySQLMeteoDatabase(MeteoDatabase):

    server: str
    database: str
    user: str
    password: str
    table: str
    columns: tuple[str, str, str, str]
    _cursors: list = field(default_factory=[], repr=False, compare=False) 

    def __post_init__(self):

        conn = connector.connect(
            host=server,
            database=database,
            user=user,
            password=password
        )
        self._cursors.append(conn.cursor())

    def mean_conditions(self, t: Time, period: TimeDelta):

        c1, c2, c3, c4 = self.columns
        dmin = (t - period / 2).isot
        dmax = (t + period / 2).isot

        query = (
            f"SELECT {c1}, {c2}, {c3}, {c4} "
            "FROM {self.table} "
            "WHERE {c1} BETWEEN '{dmin}' AND '{dmax}' "
            "ORDER by fecha" 
        )

        self._cursor.execute(query)
        t, T, P, h = list(zip(*self._cursor))

        if not len(t):
            return np.nan, np.nan, np.nan

        return np.mean(T), np.mean(P), np.mean(h)

