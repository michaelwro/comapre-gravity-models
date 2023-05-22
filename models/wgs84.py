"""

WGS84 Gravity Model

Code By: Michael Wrona
Created: 22 May 2023

"""

import numpy as np


def wgs84(lat: float, lon: float) -> np.array:
    """Compute gravitational acceleration with the WGS84 gravity equation.

    Args
    ----
        lat (float): WGS84 geodetic latitude [deg].
        lon (float): WGS84 longitude [deg].

    Returns
    -------
        (float): Gravitational acceleration in [m/s/s]

    Resources
    ---------
        https://en.wikipedia.org/wiki/Theoretical_gravity#Somigliana_equation
    """
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    frac_num = 1 + (0.00193185265245827352087 * np.sin(lat) ** 2)
    frac_den = np.sqrt(1 - (0.006694379990141316996137 * np.sin(lat) ** 2))

    return 9.780325335903891718546 * frac_num / frac_den
