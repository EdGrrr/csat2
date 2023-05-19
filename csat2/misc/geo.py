import math
import numpy as np

R_E = 6372.8  # see http://rosettacode.org/wiki/Haversine_formula


def haversine(lon1, lat1, lon2, lat2):
    """Computes the Haversine distance between two points in km"""
    dlat = math.pi * (lat2 - lat1) / 180.0
    dlon = math.pi * (lon2 - lon1) / 180.0

    lat1 = math.pi * (lat1) / 180.0
    lat2 = math.pi * (lat2) / 180.0

    arc = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon / 2) ** 2)
    c = 2 * np.arcsin(np.sqrt(arc))
    return R_E * c


def bearing(lon1, lat1, lon2, lat2):
    """Calculates the bearing between two points
    NOTE: switched lon/lat order to match TASIC code

    All arguments in degrees

    Code from: https://www.igismap.com/formula-to-find-bearing-or-heading-angle-between-two-points-latitude-longitude/"""
    dlon = np.deg2rad(lon2 - lon1)
    X = np.cos(np.deg2rad(lat2)) * np.sin(dlon)
    Y = np.cos(np.deg2rad(lat1)) * np.sin(np.deg2rad(lat2)) - np.sin(
        np.deg2rad(lat1)
    ) * np.cos(np.deg2rad(lat2)) * np.cos(dlon)
    return np.rad2deg(np.arctan2(X, Y))


def coordinate_rotate(lon, lat, centre_lon, centre_lat):
    """Rotates a set of spherical coordinates so that the centre location
    is on the equator"""
    if abs(centre_lon) < 10:
        icentre_lon = np.deg2rad(centre_lon - 90)
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
        return (np.rad2deg(np.arctan2(x2p, x1p)), np.rad2deg(np.arcsin(x3p)))
    else:
        return (np.rad2deg(np.arcsin(x3p)), -1 * np.rad2deg(np.arctan2(x2p, x1p)))


def great_circle_path(
    dep_pos, arr_pos, n_points=None, max_seg_dist=None, evenly_spaced=False
):
    """Returns a list of points along a great circle
    - Either n_points, evenly spaced (excluding ends)
    - Or a set of points spaced by no more than max_seg_dist

    evenly_spaced - points are evenly spaced (max_seg_dist method)

    Formula from http://www.edwilliams.org/avform.htm#Intermediate"""

    dist = haversine(dep_pos[0], dep_pos[1], arr_pos[0], arr_pos[1])
    if n_points and not max_seg_dist:
        f = np.linspace(0, 1, n_points + 2)
    elif max_seg_dist and not n_points:
        if dist < max_seg_dist:
            # if distance is already short, return!
            return [dep_pos, arr_pos]
        if evenly_spaced:
            f = np.linspace(0, 1, int(dist / max_seg_dist) + 2)
        else:
            f = np.arange(0, 1, max_seg_dist / dist)

    dp = np.deg2rad(dep_pos)
    ap = np.deg2rad(arr_pos)

    delta = dist / R_E
    a = np.sin((1 - f) * delta) / np.sin(delta)
    b = np.sin(f * delta) / np.sin(delta)
    x = a * np.cos(dp[1]) * np.cos(dp[0]) + b * np.cos(ap[1]) * np.cos(ap[0])
    y = a * np.cos(dp[1]) * np.sin(dp[0]) + b * np.cos(ap[1]) * np.sin(ap[0])
    z = a * np.sin(dp[1]) + b * np.sin(ap[1])

    op_lon = np.rad2deg(np.arctan2(y, x))
    op_lat = np.rad2deg(np.arctan2(z, np.sqrt(x**2 + y**2)))

    return list(zip(op_lon, op_lat))


def cross_track_error(dep_pos, arr_pos, err_pos):
    """Cross track position error from a great circle
    GC defined by dep_pos and arr_pos, err_pos is the off-track location"""

    if len(np.array(err_pos).shape) > 1:
        ep = np.array(err_pos)
        angdist_err = haversine(dep_pos[0], dep_pos[1], ep[..., 0], ep[..., 1]) / R_E
        bearing_err = np.deg2rad(
            bearing(dep_pos[0], dep_pos[1], ep[..., 0], ep[..., 1])
        )
    else:
        angdist_err = haversine(*dep_pos, *err_pos) / R_E
        bearing_err = np.deg2rad(bearing(*dep_pos, *err_pos))
    bearing_gc = np.deg2rad(bearing(*dep_pos, *arr_pos))

    return np.arcsin(np.sin(angdist_err) * np.sin(bearing_err - bearing_gc)) * R_E


def great_circle_destination(dep_pos, bearing, dist):
    dp = np.deg2rad(dep_pos)
    br = np.deg2rad(bearing)
    delta = dist / R_E
    arr_lat = np.arcsin(
        np.sin(dp[1]) * np.cos(delta) + np.cos(dp[1]) * np.sin(delta) * np.cos(br)
    )
    arr_lon = dp[0] + np.arctan2(
        np.sin(br) * np.sin(delta) * np.cos(dp[1]),
        np.cos(delta) - np.sin(dp[1]) * np.sin(arr_lat),
    )
    return list(zip(np.rad2deg(arr_lon), np.rad2deg(arr_lat)))
