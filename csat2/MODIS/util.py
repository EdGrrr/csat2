import csat2
import csat2.misc.geo
import csat2.misc.fileops
from .readfiles import readin, field_interpolate

import numpy as np
from netCDF4 import Dataset

_bands_cent = {
    # Reflective bands
    '1':  0.645,
    '2':  0.865,
    '3':  0.470,
    '4':  0.555,
    '5':  1.240,
    '6':  1.640,
    '7':  2.130,
    '8':  0.415,
    '9':  0.443,
    '10':  0.490,
    '11':  0.531,
    '12':  0.565,
    '13':  0.653,
    '14':  0.681,
    '15':  0.750,
    '16':  0.865,
    '17':  0.905,
    '18':  0.936,
    '19':  0.940,
    '26':  1.375,
    # Emissive bands
    '20':  3.750,
    '21':  3.959,
    '22':  3.959,
    '23':  4.050,
    '24':  4.465,
    '25':  4.515,
    '27':  6.715,
    '28':  7.325,
    '29':  8.550,
    '30':  9.730,
    '31': 11.030,
    '32': 12.020,
    '33': 13.335,
    '34': 13.635,
    '35': 13.935,
    '36': 14.235, }


def make_composite(gran, filename, bands=None, scale=None, offset=None, gamma=None, imscale=1, dpi=200):
    import matplotlib.pyplot as plt
    bt = crosstrack_rectify(gran.get_band_bt(31, col=61))
    btd = crosstrack_rectify(gran.get_band_bt(32, col=61)) - bt
    refl_blu = crosstrack_rectify(gran.get_band_radiance(1, col=61, refl=True))

    con_scale = np.array([1, -15, 220])[None, None, :]
    con_offset = np.array([0, 5, 190])[None, None, :]
    con_gamma = np.array([1, 1, 1])[None, None, :]

    toplot = np.concatenate((refl_blu[:, :, None],
                             btd[:, :, None],
                             bt[:, :, None]), axis=2)
    toplot = np.clip(((toplot-con_offset)/con_scale)**con_gamma, 0, 1)

    ax = plt.axes([0, 0, 1, 1])
    ax.imshow(np.rot90(toplot), aspect='auto')
    fig = plt.gcf()
    fig.set_size_inches(
        imscale*toplot.shape[0]/dpi, imscale*toplot.shape[1]/dpi)
    fig.savefig(filename, dpi=dpi)
    fig.clf()


def band_centres(band):
    '''Central wavelength of a band in um'''
    return _bands_cent[str(band)]


def bowtie_correct(field_1km):
    '''Corrects a standard 1km MODIS granule filed for the bowtie effect
by re-ordering pixels. Note this is really a costmetic change as the swath
edge pixels are still large'''
    try:
        file_correct = csat2.locator.search(
            'MODIS', 'bowtie',
            res='1km', length=field_1km.shape[0])[0]
    except IndexError:  # No files returned
        raise FileNotFoundError(
            '''Bowtie correct file not found. You need 
            to create this file first with the functions
            in csat2.MODIS.util''')

    with Dataset(file_correct) as ncdf:
        along_track_index = ncdf.variables['at_ind'][:]
        cross_track_index = ncdf.variables['ct_ind'][:]
        return field_1km[along_track_index.astype('int'),
                         cross_track_index.astype('int')]


def crosstrack_rectify(field_1km, nearest=False, missing=False):
    '''Expands out the crosstrack dimension to maintain a constant distance
between pixels. Note this is cosmetic, it doesn't account for the changes
in pixel size'''
    try:
        file_correct = csat2.locator.search(
            'MODIS', 'ctrect', res='1km')[0]
    except IndexError:  # No files returned
        raise FileNotFoundError(
            '''Crosstrack rectify file not found. You need
            to create this file first with the functions
            in csat2.MODIS.util''')
    with Dataset(file_correct) as ncdf:
        if nearest:
            remapped = ncdf.variables['remapped'][:]
            new_data = np.zeros((field_1km.shape[0], len(remapped)))
            for i in range(new_data.shape[0]):
                new_data[i] = field_1km[i][remapped.astype('int')]
        else:
            remapped_ratio = ncdf.variables['remapped_ratio'][:]
            new_data = np.zeros((field_1km.shape[0], len(remapped_ratio)))
            for i in range(new_data.shape[0]):
                new_data[i] = (
                    ((1-remapped_ratio[:, 1]) *
                     field_1km[i][remapped_ratio[:, 0].astype('int')]) +
                    ((remapped_ratio[:, 1]) *
                     field_1km[i][remapped_ratio[:, 0].astype('int')+1]))
            if missing:
                remapped = ncdf.variables['remapped'][:]
                new_mdata = np.zeros((field_1km.shape[0], len(remapped)))
                for i in range(new_data.shape[0]):
                    new_mdata[i] = field_1km[i][remapped.astype('int')]
                new_data = np.where(np.isfinite(new_data), new_data, new_mdata)
    return new_data


def crosstrack_rectcoords():
    '''Returns the coordinates in the rectified frame'''
    try:
        file_correct = csat2.locator.search(
            'MODIS', 'ctrect', res='1km')[0]
    except IndexError:  # No files returned
        raise FileNotFoundError(
            '''Crosstrack rectify file not found. You need
            to create this file first with the functions
            in csat2.MODIS.util''')
    with Dataset(file_correct) as ncdf:
        return ncdf.variables['remap_locate'][:]


def _create_bowtie_correct(mod_data, output_file):
    lon = np.zeros(mod_data['Cloud_Effective_Radius'][:, :-4].shape)
    lat = np.zeros(lon.shape)

    for i in range(lon.shape[0]//10):
        if i % 20 == 0:
            print(i)
        centre_lon = np.mean(mod_data['Longitude'][(2*i):((2*i)+2), :], axis=0)
        centre_lat = np.mean(mod_data['Latitude'][(2*i):((2*i)+2), :], axis=0)

        inc_lon = (centre_lon - mod_data['Longitude'][(2*i), :])/5
        inc_lat = (centre_lat - mod_data['Latitude'][(2*i), :])/5

        centre_lon = MODIS.field_interpolate(centre_lon[:, None])[:, 2]
        centre_lat = MODIS.field_interpolate(centre_lat[:, None])[:, 2]
        inc_lon = MODIS.field_interpolate(inc_lon[:, None])[:, 2]
        inc_lat = MODIS.field_interpolate(inc_lat[:, None])[:, 2]

        lon[(10*i):(10*(i+1))] = centre_lon + \
            np.arange(-9, 10, 2)[:, None]*inc_lon[None, :]
        lat[(10*i):(10*(i+1))] = centre_lat + \
            np.arange(-9, 10, 2)[:, None]*inc_lat[None, :]

    distance_from_baseline = csat2.misc.geo.haversine(
        lon[[0], :], lat[[0], :], lon, lat)
    along_track_index = np.zeros(distance_from_baseline.shape)
    for j in range(distance_from_baseline.shape[1]):
        sdata = list(zip(distance_from_baseline[:, j],
                         range(distance_from_baseline.shape[0])))
        sdata.sort()
        along_track_index[:, j] = list(zip(*sdata))[1]
        if j % 50 == 0:
            print(j)

    cross_track_index = np.fromfunction(
        lambda x, y: y, along_track_index.shape)

    csat2.misc.fileops.nc4_dump_mv(output_file,
                                   {'dist_bl': distance_from_baseline,
                                    'at_ind': along_track_index,
                                    'ct_ind': cross_track_index})


def _create_bowtie_correct_files():
    mod_data = MODIS.readin('MYD06_L2', 2015, 24, [
                            'Cloud_Effective_Radius'], times=['2220'])
    _create_bowtie_correct(mod_data, 'bowtie_correction_1km_2030.nc')
    mod_data = MODIS.readin('MYD06_L2', 2015, 4, [
                            'Cloud_Effective_Radius'], times=['2105'])
    _create_bowtie_correct(mod_data, 'bowtie_correction_1km_2040.nc')


def _create_crosstrack_rectify(mod_data, output_file):
    lon = MODIS.field_interpolate(mod_data['Longitude'])
    lat = MODIS.field_interpolate(mod_data['Latitude'])

    dist = csat2.misc.geo.haversine(lon[0], lat[0], lon[0, 0], lat[0, 0])
    # The nearest point in the stardard MODIS across-track grid
    remapped = np.zeros(int(dist.max()+1))+len(dist)-1
    # [0 - the point before the location,
    #  1 - the fraction of data to come from the point after]
    remapped_ratio = np.zeros((int(dist.max()+1), 2)) + \
        np.array((len(dist)-2, 1))
    # Locate MODIS pixels in the remapped image
    remap_locate = np.zeros(len(dist))

    j = 0
    k = 0
    for i in range(len(remapped)):
        try:
            if abs(i-dist[j]) > (abs(i-dist[j+1])):
                j += 1
            remapped[i] = j
            if (i > dist[k+1]):
                k += 1
            remapped_ratio[i, 0] = k
            remapped_ratio[i, 1] = (i-dist[k])/(dist[k+1]-dist[k])
        except:
            pass

    for i in range(len(remap_locate)):
        avloc = np.arange(len(remapped))[remapped == i]
        remap_locate[i] = np.mean(avloc)

    with Dataset(output_file, 'w') as ncdf:
        ncdf.createDimension('remap', len(remapped))
        ncdf.createDimension('ratio', 2)
        ncdf.createDimension('standard_modis', len(remap_locate))

        Var = ncdf.createVariable('remapped', 'f', ('remap',))
        Var[:] = remapped
        Var = ncdf.createVariable('remapped_ratio', 'f', ('remap', 'ratio'))
        Var[:] = remapped_ratio
        Var = ncdf.createVariable('remap_locate', 'f', ('standard_modis',))
        Var[:] = remap_locate


def _create_crosstrack_rectify_files():
    mod_data = MODIS.readin('MYD06_L2', 2015, 24, [
                            'Cloud_Effective_Radius'], times=['2220'], col='61')
    _create_crosstrack_rectify(mod_data, 'crosstrack_rectify_1km.nc')
