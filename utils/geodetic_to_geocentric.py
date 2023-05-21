from typing import Tuple

import math as m
from . import constants


def geodetic_to_geocentric(lat_gd_rad: float, alt_m: float) -> Tuple[float, float]:
    """Convert geodetic latitude to geocentric latitude.

    Args
    ----
        lat_gd_rad (float): WGS84 geodetic latitude [rad].
        alt_m (float): WGS84 ellipsoidal altitude [m].

    Returns
    -------
        (float, float): ({geocentric latitude [rad]}, {geocentric radius [m]})
    """

    sin_lat_gd = m.sin(lat_gd_rad)
    cos_lat_gd = m.cos(lat_gd_rad)

    N = constants.WGS84_A_M / m.sqrt(1 - (constants.WGS84_ECC2 * sin_lat_gd**2))

    rho = (N + alt_m) * cos_lat_gd
    z = (alt_m + (N * (1 - constants.WGS84_ECC2))) * sin_lat_gd

    return m.atan2(z, rho), m.sqrt(z**2 + rho**2)
