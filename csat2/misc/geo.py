import math
import numpy as np

R_E = 6372.8  # see http://rosettacode.org/wiki/Haversine_formula


def haversine(lon1, lat1, lon2, lat2):
    '''Computes the Haversine distance between two points in km'''
    dlat = math.pi*(lat2-lat1)/180.
    dlon = math.pi*(lon2-lon1)/180.

    lat1 = math.pi*(lat1)/180.
    lat2 = math.pi*(lat2)/180.

    arc = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*(np.sin(dlon/2)**2)
    c = 2*np.arcsin(np.sqrt(arc))
    return R_E*c


def bearing(lon1, lat1, lon2, lat2):
    ''' Calculates the bearing between two points
    NOTE: switched lon/lat order to match TASIC code

    All arguments in degrees

    Code from: https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/'''
    dlon = np.deg2rad(lon2-lon1)
    X = np.cos(np.deg2rad(lat2)) * np.sin(dlon)
    Y = np.cos(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2)) - \
        np.sin(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.cos(dlon)
    return np.rad2deg(np.arctan2(X, Y))


def coordinate_rotate(lon, lat, centre_lon, centre_lat):
    '''Rotates a set of spherical coordinates so that the centre location
    is on the equator'''
    if abs(centre_lon) < 10:
        icentre_lon = np.deg2rad(centre_lon-90)
    else:
        icentre_lon = np.deg2rad(centre_lon)
    icentre_lat = np.deg2rad(centre_lat)
    ilon = np.deg2rad(lon)
    ilat = np.deg2rad(lat)
    alpha = ilon - icentre_lon
    x1 = np.cos(ilat) * np.cos(alpha)
    x2 = np.cos(ilat) * np.sin(alpha)
    x3 = np.sin(ilat)

    x1p = np.cos(icentre_lat) * x1 + np.sin(icentre_lat) * x3
    x2p = x2
    x3p = -np.sin(icentre_lat) * x1 + np.cos(icentre_lat) * x3

    if abs(centre_lon) > 10:
        return (np.rad2deg(np.arctan2(x2p, x1p)),
                np.rad2deg(np.arcsin(x3p)))
    else:
        return (np.rad2deg(np.arcsin(x3p)),
                -1*np.rad2deg(np.arctan2(x2p, x1p)))
