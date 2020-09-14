import numpy as np

SOLAR_CONST = 1360.8

# Data for year 2000 (from GISS)
# https://data.giss.nasa.gov/ar5/srorbpar.html
# rae - ratio of atmosphere to earth radius, from ECHAM6
EARTH = {'ecc': 0.016704,
         'obl': np.radians(23.4398),
         'lph': np.radians(282.895),
         'rae': 0.1227e-2}


def SOLAR_CONST_adjusted(day_of_year):
    # Fixed version from wikipedia
    # https://en.wikipedia.org/wiki/Solar_irradiance
    return SOLAR_CONST * (1 + 0.034 * np.cos(2 * np.pi * day_of_year / 365.25))


def solar_zenith_angle_angles(observer_latitude,
                              subsolar_latitude,
                              hour_angle):
    sza = np.arccos(
        (np.sin(observer_latitude) *
         np.sin(subsolar_latitude)) +
        (np.cos(observer_latitude) *
         np.cos(subsolar_latitude) *
         np.cos(hour_angle)))
    return sza


def solar_zenith_angle_time(observer_latitude,
                            day_of_year,
                            local_solar_time):
    return solar_zenith_angle_angles(
        observer_latitude,
        declination_doy(day_of_year),
        (local_solar_time / 12 - 1) * np.pi)


def declination_doy(day_of_year):
    equinox_adj = (EARTH['lph']-np.radians(270))
    return (-1 * EARTH['obl'] *
            np.cos((2 * np.pi * (day_of_year-1) / 365.25) +
                   equinox_adj))


def sunrise_hour_angle(observer_latitude,
                       subsolar_latitude):
    cosh0 = - np.tan(observer_latitude) * np.tan(subsolar_latitude)
    cosh0 = np.clip(cosh0, -1, 1)
    return np.arccos(cosh0)


def sunrise_time(observer_latitude,
                 day_of_year):
    hour_angle = sunrise_hour_angle(
        observer_latitude,
        declination_doy(day_of_year))
    return 12 - (12 * np.abs(hour_angle) / np.pi)


def insolation(observer_latitude,
               day_of_year,
               local_solar_time):
    return SOLAR_CONST_adjusted(day_of_year) * np.cos(
        solar_zenith_angle_time(observer_latitude,
                                day_of_year,
                                local_solar_time)).clip(0, None)


def daily_mean_insolation(observer_latitude,
                          day_of_year):
    subsolar_latitude = declination_doy(day_of_year)
    h0 = sunrise_hour_angle(observer_latitude,
                            subsolar_latitude)
    qmean = (SOLAR_CONST_adjusted(day_of_year) / np.pi *
             ((h0 * np.sin(observer_latitude) *
               np.sin(subsolar_latitude)) +
              (np.cos(observer_latitude) *
               np.cos(subsolar_latitude) *
               np.sin(h0))))
    return qmean
