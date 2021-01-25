''' A library of functions for converting between different time
    specifications (UTC, TAI and LST)
    and different time formats (year, month, day and year, day of year - doy)'''

import datetime
import numpy as np
from copy import copy

#################
# DOY functions #
#################


def doy_to_date(year, doy):
    '''doy_to_date(year,doy)
    Converts a date from DOY representation to day.month.year
    returns tuple(year,month,day)

    Raises ValueError if the doy is not valid for that year'''
    dat = datetime.date(year, 1, 1)
    dat += datetime.timedelta(int(doy) - 1)
    if dat.year != year:
        raise ValueError('Day not within year')
    return (dat.year, dat.month, dat.day)


def date_to_doy(year, month, day):
    '''Converts date format from y,m,d to a tuple (year,doy)'''
    return (year, datetime.date(year, month, day).timetuple().tm_yday)


def doy_step(year, doy, step):
    '''Adds "step" number of days to the date specified by {year,doy},
    taking into account lengths of yaers etc.'''
    dat = datetime.date(year, 1, 1)
    dat += datetime.timedelta(int(doy) - 1 + step)
    return date_to_doy(dat.year, dat.month, dat.day)


def ydh_to_datetime(year, doy, hour):
    dat = (datetime.datetime(*doy_to_date(year, doy)) +
           datetime.timedelta(
               hours=int(hour//1),
               minutes=int((hour % 1)*60)))
    return dat


def datetime_to_ydh(dtime):
    year = dtime.year
    doy = dtime.timetuple().tm_yday
    dtime_diff = datetime.datetime(*doy_to_date(year, doy))
    hour = (dtime-dtime_diff).total_seconds()/3600
    return year, doy, hour


def datetime_to_lst(dtime, lon):
    '''Returns the local solar time for a given datetime and longitude'''
    return utc_to_lst(dtime.hour+dtime.minute/60, lon)


def doy_exists(year, doy):
    "Returns True if the doy is a valid date, False otherwise"
    try:
        y, m, d = doy_to_date(year, doy)
        return True
    except ValueError:
        # doy_to_date raises ValueError for an invalid doy
        return False


def get_season(doy, year_length=364):
    '''Returns the season for a given year and doy
    0 - DJF
    1 - MAM
    2 - JJA
    3 - SON

    Note, this assumes a 364 day year for simplicity. It should not 
    matter for the majority of use cases (seasons are still within
    the range above)'''
    return ((doy+30) % year_length)//91


#####################
# TAI/UTC functions #
#####################
def tai_to_utc(time_tai):
    '''Returns the UTC given a decimal TAI time (seconds since Jan 1st 1993)'''
    time_utc = datetime.datetime(1993, 1, 1) + \
        datetime.timedelta(seconds=time_tai)
    return (time_utc.year,
            time_utc.timetuple().tm_yday,
            time_utc.hour + time_utc.minute / 60. + time_utc.second / 3600.)


def utc_to_tai(year, doy, time_utc):
    '''Returns the TAI given the year, day of year and decimal utc time'''
    td = ((datetime.datetime(year, 1, 1) + datetime.timedelta(doy - 1) +
           datetime.timedelta(hours=time_utc)) - datetime.datetime(1993, 1, 1))
    return (td.seconds + td.days * 24 * 3600)


def lst_to_utc(time_lst, longitude):
    '''Returns the UTC time given a decimal LST and longitude'''
    longitude = np.mod(longitude, 360.)
    return np.mod(time_lst - longitude * 24. / 360, 24)


def utc_to_lst(time_utc, longitude):
    '''Returns the LST time given a decimal UTC and longitude'''
    longitude = np.mod(longitude, 360.)
    return np.mod(time_utc + longitude * 24. / 360, 24)


def utc_to_sat_offset(utc, lon, sattime, col='5'):
    '''Utc [0,24), lon [0,360)
    Returns hours to satellite track on day doy days ago

    Sattime (and utc) is a decimal hour - e.g. 13.5 for Aqua

    The satellite is assumed to have a time like MODIS L3 (C5)
    (e.g. DOY is defined as the UTC day)

    For collection 6, use flag 'col='6'', see example below


    For example, using Aqua (sattime=13.5)

    At longitude 0, utc 14, the offset is 0.5
               (half and hour behind the satellite)

    At longitude -170, utc 1, the offset is 0.16 (approximately 10 min behind)

    At longitude -170, utc 23, the offset is 22.16 (22 hours behind)

    import matplotlib.pyplot as plt
    t,lon = np.meshgrid(np.arange(0, 24, 0.25), np.arange(0, 360, 1))
    offset = utc_to_sat_offset(t, lon, 13.5, col='6')
    plt.imshow(offset)
    plt.show()

'''
    lon[lon < 0] += 360
    lst = utc + lon / 15.
    lst = lst - sattime
    if col == '5':
        lst[(lon / 15.) > sattime] -= 24
    else:
        # Collection 6 modifies the dateline so that the effective
        # switch time is 12, even though the satellties are earlier and later
        # This is to avoid data gaps close to the dateline
        lst[(lon / 15.) > 12] -= 24
    return lst


def toLocalSolarTime(lst_time, gmt_times, longitudes, data,
                     interpolation='linear', DEBUG=False):
    ''' Converts a set of data at GMT times (gmt_times) and 
    longitudes (longitudes) to local solar time (lst_time),

    Input
    lst_time - number defining local solar time required
                 (between 0 and 24, fractions for minutes)
    gmt_times[times] - array indicating the GMT times of the data slices
    longitudes[lons] - array containing the longitudes of the grid
    data[lats,lons,times] - containing the data to be retimed
    interpolation (linear or nearest) - specifies the method
                 to use when retiming

    Output - array[lats,lons] of the data retimed to lst_time
    '''
    if len(set(gmt_times)) != len(gmt_times):
        raise ValueError(
            'Data should not contain duplicate times' +
            '(gmt_times,toLocalSolarTime))')
    outdata = np.zeros(data[:, :, 0].shape)

    # Deal with non-sorted gmt_times
    times_sort = copy(gmt_times)
    times_sort.sort()
    times_index = np.zeros(gmt_times.shape)
    for i in range(len(times_index)):
        times_index[i] = list(gmt_times).index(times_sort[i])

    # Do the re-timing
    for lon_ind in range(len(longitudes)):
        local_gmt = lst_to_utc(lst_time, longitudes[lon_ind])

        if local_gmt > times_sort[-1]:
            raise ValueError(
                'Cannot calculate LST, require extra data timeslice')
        # Calculate the weights for each timeslice
        weights = np.zeros(gmt_times.shape)
        for i in range(len(times_sort)):
            if times_sort[i] >= local_gmt:
                break
        high_ind = i

        for i in range(high_ind, -1, -1):
            if times_sort[i] <= local_gmt:
                break
        low_ind = i

        if high_ind == low_ind:
            weights[low_ind] = 1
        else:
            if interpolation == 'linear':
                fract = ((local_gmt - times_sort[high_ind]) /
                         (times_sort[low_ind] - times_sort[high_ind]))
                weights[low_ind] = fract
                weights[high_ind] = 1 - fract
            elif interpolation == 'nearest':
                if ((abs(longitudes[lon_ind] - longitudes[low_ind]) <
                     abs(longitudes[lon_ind] - longitudes[low_ind]))):
                    weights[low_ind] = 1
                else:
                    weights[high_ind] = 1

        outdata[:, lon_ind] = (data[:, lon_ind, :] *
                               weights[np.newaxis, :]).sum(axis=1)
        if interpolation == 'nearest':
            outdata[:, lon_ind] = np.nansum(
                (data[:, lon_ind, :] * weights[np.newaxis, :]), axis=1)
    return outdata
