from astropy.time import Time, TimeDelta

DAY = TimeDelta(1, format='jd')

# Accept more formats than astropy.time.Time

def date(s: str) -> Time:
        
    d = s

    if ':' not in s and '-' not in s: # compact date format

        if len(s) in [7, 9, 11, 13]: # RINEX3 yyyydddHHMM format
            d = f"{s[0:4]}:{s[4:7]}"
            s = ":" + s[7:]
        elif len(s) in [8, 10, 12, 14]: # yyyymmddHHMMSS
            d = f"{s[0:4]}-{s[4:6]}-{s[6:8]}"
            s = "T" + s[8:]
            
        if len(s) == 3:
            s += "0000"
        elif len(s) == 5:
            s += "00"

        d += f"{s[0:3]}:{s[3:5]}:{s[5:7]}"

    return Time(d)

def format(date: str | Time, fmt: str) -> str:
    """format_date(date, fmt)

Arguments

    date [string or astropy.Time]
        Date in any format readable by astropy.Time

    fmt [string]

        Desired format:
            rnx - RINEX 3 format "YYYYMMDDhhmm"
            pride - PRIDE PPP-AR format "YYYYJJJhhmm" using day-of-year
            iso - "YYYY-MM-DD"
            pdp3 - "YYYY/MM/DD"
            ymd - "YYYYMMDD"
        
    """
    if not isinstance(date, Time):
        date = Time(date)

    if fmt == 'rnx':
        return date.strftime("%Y%j%H%M")

    if fmt == 'pride':
        return date.strftime("%Y%j")

    if fmt == 'iso':
        return date.strftime("%Y-%m-%d")

    if fmt == 'pdp3':
        return date.strftime("%Y/%m/%d")

    if fmt == "ymd":
        return date.strftime("%Y%m%d")

    if fmt == "jyear":
        return date.jyear
        
    if fmt == "mjd":
        return date.mjd
    raise NotImplementedError(f"Unknown date format: {fmt}")

