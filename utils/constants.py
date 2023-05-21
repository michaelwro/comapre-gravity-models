"""
Physical constants

Code By: Michael Wrona
Created: 21 May 2023

Resources
* https://en.wikipedia.org/wiki/Geopotential_model

"""

WGS84_A_M = 6378137  # semi-major axis [m]
WGS84_ECC = 0.0818191908426215  # eccentricity
WGS84_ECC2 = WGS84_ECC**2  # eccentricity squared
WGS84_GM_M3PS2 = 3.986004418e14  # Gravitational parameter [m^3 / s^2].

J2 = -0.1082635854e-2  # J2 zonal harmonic
J3 = 0.2532435346e-5  # J3 zonal harmonic
J4 = 0.1619331205e-5  # J4 zonal harmonic
