"""

Convert WGS84 lat/lon/alt to ECEF position.

Code By: Michael Wrona
Created: 21 May 2023.

"""

import numpy as np
from . import constants


def lla_to_ecef(lat: float, lon: float, alt: float) -> np.array:
    """Convert WGS84 lat/lon/alt to ECEF position.

    Args
    ----
        lat (float): WGS84 geodetic latitude [rad].
        lon (float): WGS84 longitude [rad].
        alt (float): WGS84 ellipsoidal altitude [m].

    Returns
    -------
        (numpy.array): ECEF position [x, y, z] in [m].
    """

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    c_term = constants.WGS84_A_M / np.sqrt(1 - (constants.WGS84_ECC2 * sin_lat**2))
    s_term = c_term * (1 - constants.WGS84_ECC2)

    return np.array(
        [
            (c_term + alt) * cos_lat * np.cos(lon),
            (c_term + alt) * cos_lat * np.sin(lon),
            (s_term + alt) * sin_lat,
        ]
    )
