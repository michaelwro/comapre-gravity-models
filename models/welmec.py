"""

WELMEC gravity model.

Code By: Michael Wrona
Created: 22 May 2023s

"""

import numpy as np


def welmec(lat: float, lon: float) -> float:
    """Compute gravitational acceleration with the WELMEC gravity equation.

    Args
    ----
        lat (float): WGS84 geodetic latitude [deg].
        lon (float): WGS84 longitude [deg].

    Returns
    -------
        (float): Gravitational acceleration in [m/s/s]

    Resources
    ---------
        https://en.wikipedia.org/wiki/Theoretical_gravity#WELMEC_formula
    """
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    # ignore height part
    return 9.780318 * (
        1 + (0.0053024 * np.sin(lat) ** 2) - (0.0000058 * np.sin(2 * lat) ** 2)
    )
