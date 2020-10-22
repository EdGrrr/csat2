'''Stats functions'''

import numpy as np
import scipy.stats


def nanlinregress(data_x, data_y):
    ''' A version of the scipy linear regression function that can cope
    with missing data - by ignoring it'''
    mask = np.isfinite(data_x+data_y)
    return scipy.stats.linregress(data_x[mask], data_y[mask])


def nanmean(data, *args, **kwargs):
    '''Returns of mean of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data, mask=np.bitwise_not(np.isfinite(data))).mean(*args, **kwargs),
                        np.nan)


def nanstd(data, *args, **kwargs):
    '''Returns of std of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data,
        mask=np.bitwise_not(np.isfinite(data))).std(*args, **kwargs), np.nan)


def nanmax(data, *args, **kwargs):
    '''Returns of max of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data,
        mask=np.isnan(data)).max(*args, **kwargs), np.nan)


def nanmin(data, *args, **kwargs):
    '''Returns of min of a numpy array ignoring missing data'''
    return np.ma.filled(np.ma.array(
        data,
        mask=np.isnan(data)).min(*args, **kwargs), np.nan)


def nanmask(data):
    '''Returns a masked array with the nans masked'''
    return np.ma.array(data, mask=np.isnan(data))


def zero_nans(data):
    '''Returns an array where nans have been replaced by zeros'''
    return np.where(np.isnan(data), 0, data)


def lin_av(data):
    '''Returns an array that is the average of each neighbouring pair'''
    data = np.array(data)
    return (data[:-1]+data[1:])/2.


def lat_weighted_av(data, latrange, lat_axis=0):
    '''Returns a latitude weighted mean of the data.  Kind of approximate, but
    perhaps still good-ish.  Nans are treated as missing and
    ignored from the average (since 21/9/15).

    data : array_like, shape(N1,N2)
        An array (lat-lon) with regularly spaced latitudes, to be averaged
    latrange : array_like, shape(2)
        The min and max latitudes, used to generate the weights
    lat_axis : Int - default = 0
        !!!!!!The latitude axis of the array - not yet working!!!!!'''
    nlats = data.shape[lat_axis]
    lats = np.linspace(latrange[0], latrange[1], nlats+1)
    lats = lin_av(lats*np.pi/180.)
    lats = np.cos(lats)
    lats = lats/np.sum(lats)
    lats = lats[:, None].repeat(data.shape[1], axis=1)
    lats[np.isnan(data)] = np.nan
    return np.nansum(data*lats)/np.nansum(lats)


def reduce_res(data, res, axis=(0, 1), func=np.ma.mean):
    '''Reduces the resolution along axis. Currently only the first two dimensions'''
    dshape = data.shape
    newdshape = [dshape[0]//res, res, dshape[1]//res, res]
    return func(data.reshape(newdshape), axis=(1, 3))


def reduce_res_old(data, res, func=np.ma.mean):
    dshape = data.shape
    outdata = np.zeros((dshape[0]//res, dshape[1]//res))
    for i in range(outdata.shape[0]):
        for j in range(outdata.shape[1]):
            outdata[i, j] = func(np.ma.array(
                data[res*i:res*(i+1), res*j:res*(j+1)],
                mask=np.isnan(data[res*i:res*(i+1), res*j:res*(j+1)])))
    return outdata


def reduce_res_axis(data, res, axis=0, func=nanmean):
    '''Reduce the resolution of an array along a specified axis (axis) by
    a factor (res), using a function (func).'''
    if not isinstance(res, int):
        print('res converted to int')
        res = int(res)
    dshape = data.shape
    if axis < 0:
        axis = len(dshape) + axis
    dshape = np.concatenate((dshape[:axis], np.array(
        [dshape[axis]//res, res]), dshape[(axis+1):]))
    return func(data.reshape(dshape.astype(int)), axis=axis+1)


def normalise(data, axis=0):
    '''Normalises a histogram along a specified axis'''
    return data/np.expand_dims(data.sum(axis=axis), axis=axis)


def regression_gridded(data_x, data_y, weights=None, axis=0, min_no=3):
    '''Returns an array of the slopes, errors on the slopes and pearson's pmcc
    for a set of gridded data

    Parameters
    ---------------
    data_x : array_like, shape(N1,N2,N3)
       An array of data to be considered the independant variable,
    with missing data signified by np.nan

    data_y : array_like, shape(N1,N2,N3)
       An array of data to be considered the dependant variable,
    with missing data signified by np.nan

    weights : array_like, shape(N1,N2,N3)
       An array of data to be considered the weights,
    with missing data signified by np.nan.  These are treated as multiplying
    factors in the simple least squares formula - it is equivalent to
    multiplying the number of points by some multiple of the weights
    (at least for a simple version),  This needs to be checked.

    axis : integer
       The axis along which the regression should be performed

    min_no : integer
       The minimum number of valid points necessary for a regression to be
    performed. Minimum value is 3
    '''
    if weights is None:
        weights = np.ones(data_x.shape)

    mask = np.where(np.logical_not(np.isfinite(data_x+data_y+weights)))
    data_x[mask] = np.nan
    data_y[mask] = np.nan
    weights[mask] = np.nan
    n = np.isfinite(data_x+data_y).sum(axis=axis)
    nmask = (n <= min_no)
    n[nmask] = 3
    sum_x = np.nansum(weights*data_x, axis=axis)
    sum_y = np.nansum(weights*data_y, axis=axis)
    sum_y2 = np.nansum(weights*(data_y**2), axis=axis)
    sum_x2 = np.nansum(weights*(data_x**2), axis=axis)
    sum_xy = np.nansum(weights*data_x*data_y, axis=axis)
    sxy = (sum_xy - (sum_x*sum_y)/n)
    sxy[nmask] = 1
    sxx = (sum_x2 - (sum_x**2)/n)
    sxx[nmask] = 1
    syy = (sum_y2 - (sum_y**2)/n)
    syy[nmask] = 2
    return np.where((n > min_no),
                    (np.array((sxy/sxx)),
                     np.array(np.sqrt((1./(n-2))*((syy/sxx)-((sxy/sxx)**2)))),
                     np.array((sxy/np.sqrt(sxx*syy)))),
                    np.nan)


def regression_gridded_reed(data_x, data_y, err_x, err_y, axis=0, min_no=0):
    '''Returns an array of the slopes, errors on the slopes and pearson's pmcc
    for a set of gridded data

    Parameters
    ---------------
    data_x : array_like, shape(N1,N2,N3)
       An array of data to be considered the independant variable,
    with missing data signified by np.nan

    data_y : array_like, shape(N1,N2,N3)
       An array of data to be considered the dependant variable,
    with missing data signified by np.nan

    axis : integer
       The axis along which the regression should be performed

    min_no : integer
       The minimum number of valid points necessary for a
    regression to be performed
    '''

    m = regression_gridded(data_x, data_y, axis=axis, min_no=min_no)

    mask = np.where(np.logical_not(np.isfinite(data_x+data_y)))
    data_x[mask] = np.nan
    data_y[mask] = np.nan
    n = np.isfinite(data_x+data_y).sum(axis=axis)
    sum_x = np.nansum(data_x, axis=axis)
    sum_y = np.nansum(data_y, axis=axis)
    sum_y2 = np.nansum(data_y**2, axis=axis)
    sum_x2 = np.nansum(data_x**2, axis=axis)
    sum_xy = np.nansum(data_x*data_y, axis=axis)
    sxy = (sum_xy - (sum_x*sum_y)/n)
    sxx = (sum_x2 - (sum_x**2)/n)
    syy = (sum_y2 - (sum_y**2)/n)
    wx = 1./(err_x**2)
    wy = 1./(err_y**2)

    W = wy/(1+m**2 * (wy/wx))
    A = (W**2 / wx)*sxy
    B = (W**2) * (sxx/wy - syy/wx)
    C = -(wx/wy)*A
    discrm = np.sqrt(B**2 - 4*A*C)
    roots = np.concatenate((((-B + discrm)/(2*A))[None, :],
                            ((-B - discrm)/(2*A))[None, :]), axis=0)
    return roots


def gridbox_regression(data_x, data_y, lat, lon, bins):
    '''Returns an array of pearson's correlation coefficients from 1d arrays
    of each data type. Use for ungridded data, regression_gridded is
    for gridded data

    Parameters
    ---------------
    data_x : array_like, shape(N,)
       A sequence of values to be considered as the independant coordinate
    for calculating r

    data_y : array_like, shape(N,)
       A sequence of values to be considered as the dependant coordinate for
    calculating r

    lat : array_like, shape(N,)
       The latitude values for each of the points

    lon : array_like, shape(N,)
       The longitude values for each of the points

    bins : int or [int, int] or array-like or [array, array], optional
      The bin specification:

        * the number of bins for the two dimensions (nx=ny=bins),
        * the number of bins in each dimension (nx, ny = bins),
        * the bin edges for the two dimensions (x_edges=y_edges=bins),
        * the bin edges in each dimension (x_edges, y_edges = bins).'''

    mask = np.where(np.isfinite(data_x+data_y))
    total = np.histogram2d(lat[mask], lon[mask], bins=bins)[0]
    sum_x = np.histogram2d(lat[mask], lon[mask],
                           bins=bins, weights=data_x[mask])[0]
    sum_y = np.histogram2d(lat[mask], lon[mask],
                           bins=bins, weights=data_y[mask])[0]
    sum_xx = np.histogram2d(
        lat[mask], lon[mask], bins=bins, weights=data_x[mask]*data_x[mask])[0]
    sum_xy = np.histogram2d(
        lat[mask], lon[mask], bins=bins, weights=data_x[mask]*data_y[mask])[0]
    return ((total*sum_xy - sum_x*sum_y).astype('float32') /
            (total*sum_xx-sum_x*sum_x).astype('float32'))


def partial_corr_gridded(data_x, data_y, data_covar, *args, **kwargs):
    '''Calculates the partial correlation between data_x and data_y assuming
    data_covar as the controlling variable
    Uses regression_gridded, so only for gridded data'''
    r_ab = regression_gridded(data_x, data_y, *args, **kwargs)[2]
    r_ac = regression_gridded(data_x, data_covar, *args, **kwargs)[2]
    r_bc = regression_gridded(data_covar, data_y, *args, **kwargs)[2]
    number = np.isfinite(data_x+data_y+data_covar).sum(axis=0)
    return ((r_ab - r_ac*r_bc) / (np.sqrt(1 - r_ac**2) * np.sqrt(1 - r_bc**2)),
            number)
