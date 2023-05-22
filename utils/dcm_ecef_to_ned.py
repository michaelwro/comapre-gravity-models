"""

Compute ECEF to NED DCM.

Code By: Michael Wrona
Created: 22 May 2023

"""

import numpy as np


def dcm_ecef_to_ned(lat: float, lon: float) -> np.array:
    """Compute DCM to rotate a vector from ECEF to NED.

    Args
    ----
        lat (float): WGS84 geodetic latitude in [rad].
        lon (float): WGS84 geodetic longitude in [rad].

    Returns
    -------
        (np.array): ECEF to NED DCM.

    Resources
    ---------
        https://www.mathworks.com/help/aeroblks/directioncosinematrixeceftoned.html
    """

    slat = np.sin(lat)
    clat = np.cos(lat)
    slon = np.sin(lon)
    clon = np.cos(lon)

    return np.array(
        [
            [-slat * clon, -slat * slon, clat],
            [-slon, clon, 0],
            [-clat * clon, -clat * slon, -slat],
        ]
    )
