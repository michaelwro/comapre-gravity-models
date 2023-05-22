"""

Zonal Harminic Gravity Model

Code By: Michael Wrona
Created: 21 May 2023

"""

from typing import Tuple
import numpy as np
import utils
from utils import constants

# coefs from * https://en.wikipedia.org/wiki/Geopotential_model
J2 = -0.1082635854e-2  # J2 zonal harmonic (BM&W pg. 422)
J3 = 0.2532435346e-5  # J3 zonal harmonic (BM&W pg. 422)
J4 = 0.1619331205e-5  # J4 zonal harmonic (BM&W pg. 422)


def zonal(lat: float, lon: float, order=2) -> float:
    """Compute gravitational acceleration from zonal harmonics.

    Args
    ----
        lat (float): WGS84 geodetic latitude [deg].
        lon (float): WGS84 longitude [deg].
        order (int): Order of model. J2=2, J3=3, etc.

    Returns
    -------
        (float): Gravitational acceleration in [m/s/s]

    Resources
    ---------
        "Fundamentals of Astrodynamics and Applications" (Vallado) section 8.7.1.
    """
    max_order = 4

    if order < 2:
        raise ValueError(f"Order must be between 2 and {max_order}.")

    if order > max_order:
        raise ValueError(f"Zonal harmonics supported up to J{max_order}.")

    harmonic_functions = [_zonal_j2, _zonal_j3]

    # cvt. to [rad]
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # convert LLA to ECEF
    r_vec = utils.lla_to_ecef(lat, lon, 0)
    x, y, z = r_vec
    r = np.linalg.norm(r_vec)

    # compute grav. accel.
    x_accel = -constants.WGS84_GM_M3PS2 * x / r**3
    y_accel = -constants.WGS84_GM_M3PS2 * y / r**3
    z_accel = -constants.WGS84_GM_M3PS2 * z / r**3
    for harmonic_funct in harmonic_functions:
        x_term, y_term, z_term = harmonic_funct(x, y, z)
        x_accel -= x_term
        y_accel -= y_term
        z_accel -= z_term

    grav_ecef = np.array([x_accel, y_accel, z_accel])
    grav_ned = utils.dcm_ecef_to_ned(lat, lon) @ grav_ecef

    return grav_ned[2]


def _zonal_j2(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Compute J2 gravitational acceleration contribution."""
    re = constants.WGS84_A_M
    mu = constants.WGS84_GM_M3PS2
    r = np.sqrt(x**2 + y**2 + z**2)

    factor = -3 * J2 * mu * re**2 / (2 * r**5)
    x_accel = factor * x * (1 - (5 * z**2 / r**2))
    y_accel = factor * y * (1 - (5 * z**2 / r**2))
    z_accel = factor * z * (3 - (5 * z**2 / r**2))

    return x_accel, y_accel, z_accel


def _zonal_j3(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Compute J3 gravitational acceleration contribution."""
    re = constants.WGS84_A_M
    mu = constants.WGS84_GM_M3PS2
    r = np.sqrt(x**2 + y**2 + z**2)

    factor = -5 * J3 * mu * re**3 / (2 * r**7)
    x_accel = factor * x * ((3 * z) - (7 * z**3 / r**2))
    y_accel = factor * y * ((3 * z) - (7 * z**3 / r**2))
    z_accel = factor * ((6 * z**2) - (7 * z**4 / r**2) - (3 * r**2 / 5))

    return x_accel, y_accel, z_accel


def _zonal_j4(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Compute J4 gravitational acceleration contribution."""
    re = constants.WGS84_A_M
    mu = constants.WGS84_GM_M3PS2
    r = np.sqrt(x**2 + y**2 + z**2)

    factor = 15 * J4 * mu * re**4 / (8 * r**7)
    x_accel = factor * x * (1 - (14 * z**2 / r**2) + (21 * z**4 / r**4))
    y_accel = factor * y * (1 - (14 * z**2 / r**2) + (21 * z**4 / r**4))
    z_accel = factor * z * (5 - (70 * z**2 / (3 * r**2)) + (21 * z**4 / r**4))

    return x_accel, y_accel, z_accel
