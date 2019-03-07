import pandas as pd
from urllib.request import urlopen

def get_data(jyear, jday, jperiod, lonmin, lonmax, latmin, latmax):

    cols = [
        "UNIX_Time", "Sat", "Year", "Mo", "Dy", "Hr", "Mn",
        "Longitude", "Latitude", "B21", "B22", "B6", "B31", "B32",
        "SatZen", "SatAzi", "SunZen", "SunAzi", "Line", "Samp",
        "Ratio", "Glint", "Excess", "Temp", "Err"]

    url = (f"http://modis.higp.hawaii.edu/cgi-bin/mergeimage?maptype=alerts&"
            f"jyear={jyear}&"
            f"jday={jday}&"
            f"jperiod={jperiod}&"
            f"lonmin={lonmin}&"
            f"latmin={latmin}&"
            f"lonmax={lonmax}&"
            f"latmax={latmax}")

    df = pd.read_csv(urlopen(url),
            delimiter=' ', skipinitialspace=True, header=None)

    df.columns = cols

    return df
