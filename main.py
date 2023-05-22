"""

Compute and plot gravity for various gravity models.

Code By: Michael Wrona
Created: 20 May 2023

"""

import models
import numpy as np
import matplotlib.pyplot as plt

latitudes = list(np.linspace(-90, 90, 100))

zonal_j2 = []
zonal_j3 = []
zonal_j4 = []
wgs84 = []

for lat in latitudes:
    zonal_j2_grav = models.zonal(lat, 0, order=2)
    zonal_j2.append(np.linalg.norm(zonal_j2_grav))

    zonal_j3_grav = models.zonal(lat, 0, order=3)
    zonal_j3.append(zonal_j3_grav)

    zonal_j4_grav = models.zonal(lat, 0, order=4)
    zonal_j4.append(zonal_j4_grav)

    wgs84_grav = models.wgs84(lat, 0)
    wgs84.append(wgs84_grav)

plt.figure()
plt.plot(latitudes, zonal_j2, label="J2")
plt.plot(latitudes, zonal_j3, label="J3")
plt.plot(latitudes, zonal_j4, label="J4")
plt.plot(latitudes, wgs84, label="WGS84")

plt.title("Comparing Gravity Models")
plt.xlabel("Latitude [deg]")
plt.ylabel(r"Gravitational Accel. [$m/s^2$]")
plt.legend()
plt.grid()

plt.savefig("gravity.png")
# plt.show()
