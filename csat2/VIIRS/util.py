from .. import locator
from .readfiles import readin
import numpy as np
import csat2.misc.geo
import csat2.misc.fileops
from netCDF4 import Dataset
import pkg_resources

_bands_cent = {
    # Imager bands (375m)
    "I01": 0.640,
    "I02": 0.865,
    "I03": 1.61,
    "I04": 3.74,
    "I05": 11.45,
    # Dual gain mod (750x250m)
    "M01": 0.412,
    "M02": 0.445,
    "M03": 0.488,
    "M04": 0.555,
    "M05": 0.672,
    "M07": 0.865,
    "M13": 4.05,
    # Moderate resolution (750m)
    "M06": 0.746,
    "M08": 1.24,
    "M09": 1.38,
    "M10": 1.61,
    "M11": 2.25,
    "M12": 3.70,
    "M14": 8.55,
    "M15": 10.76,
    "M16": 12.01,
    "DNB": 0.7,
}

_bands_res = {
    # Imager bands (375m)
    "I01": "375m",
    "I02": "375m",
    "I03": "375m",
    "I04": "375m",
    "I05": "375m",
    # Dual gain mod (750x250m)
    "M01": "750m250m",
    "M02": "750m250m",
    "M03": "750m250m",
    "M04": "750m250m",
    "M05": "750m250m",
    "M07": "750m250m",
    "M13": "750m250m",
    # Moderate resolution (750m)
    "M06": "750m",
    "M08": "750m",
    "M09": "750m",
    "M10": "750m",
    "M11": "750m",
    "M12": "750m",
    "M14": "750m",
    "M15": "750m",
    "M16": "750m",
    "DNB": "750m",
}


def band_centres(band):
    return _bands_cent[str(band)]


def band_res(band):
    return _bands_res[str(band)]


def _remap_by_file(datafield, filetype, res):
    file_correct = pkg_resources.resource_filename(
            "csat2", f"data/viirs_remapping/{filetype}_correction_{res[:4]}.nc")
    with Dataset(file_correct) as ncdf:
        along_track_index = ncdf.variables["at_ind"][:]
        cross_track_index = ncdf.variables["ct_ind"][:]
    scan_width = along_track_index.shape[0]
    atind = np.fromfunction(
        lambda x, y: along_track_index[x % scan_width, y] + x - x % scan_width,
        datafield.shape,
        dtype=int,
    ).astype("int")
    ctind = np.fromfunction(
        lambda x, y: cross_track_index[x % scan_width, y], datafield.shape, dtype=int
    ).astype("int")
    return datafield[atind, ctind]


def bowtie_correct(datafield, res="375m"):
    return _remap_by_file(datafield, "bowtie", res=res)


def mesh_correct(datafield, res="375m"):
    return _remap_by_file(datafield, "mesh", res=res)


def _create_mesh_correct(img_data, mod_data, output_file, res=375):
    scan_width = int(12000 // res)

    img_vals = np.isfinite(img_data[:scan_width, :])
    img_min = np.argmax(img_vals, axis=0)

    along_track_index = np.fromfunction(lambda x, y: x, img_vals.shape, dtype="i")
    for i in range(1, scan_width):
        rep_inds = np.where(img_min == i)[0]
        along_track_index[:i, rep_inds] = i
        along_track_index[(-i - 1) :, rep_inds] = scan_width - i - 1

    cross_track_index = np.fromfunction(
        lambda x, y: y, along_track_index.shape, dtype=np.int
    )

    csat2.misc.fileops.nc4_dump_mv(
        output_file, {"at_ind": along_track_index, "ct_ind": cross_track_index}
    )


def _create_bowtie_correct(img_data, mod_data, output_file, res=375):
    scan_width = int(12000 // res)

    centre_scan_row = int(np.floor(scan_width // 2))
    inc_lon = (
        mod_data["longitude"][centre_scan_row + scan_width]
        - mod_data["longitude"][centre_scan_row]
    ) / scan_width
    inc_lat = (
        mod_data["latitude"][centre_scan_row + scan_width]
        - mod_data["latitude"][centre_scan_row]
    ) / scan_width

    lon = (
        np.arange(-np.floor(scan_width // 2), -np.floor(scan_width // 2) + scan_width)[
            :, None
        ]
        * inc_lon[None, :]
    ) + mod_data["longitude"][centre_scan_row][None, :]
    lat = (
        np.arange(-np.floor(scan_width // 2), -np.floor(scan_width // 2) + scan_width)[
            :, None
        ]
        * inc_lat[None, :]
    ) + mod_data["latitude"][centre_scan_row][None, :]

    old_lon = mod_data["longitude"][:scan_width, :]
    old_lon[np.isnan(img_data[:scan_width, :])] = -90
    old_lat = mod_data["latitude"][:scan_width, :]
    old_lat[np.isnan(img_data[:scan_width, :])] = -90

    along_track_index = np.zeros(lon.shape)
    for j in range(scan_width):
        along_track_index[j, :] = np.argmin(
            csat2.misc.geo.haversine(lat[j, :], lon[j, :], old_lat, old_lon), axis=0
        )

    cross_track_index = np.fromfunction(
        lambda x, y: y, along_track_index.shape, dtype=np.int
    )

    csat2.misc.fileops.nc4_dump_mv(
        output_file, {"at_ind": along_track_index, "ct_ind": cross_track_index}
    )


def _create_bowtie_correct_files():
    mod_data = readin("VNP03IMG", 2015, 80, ["longitude", "latitude"], time="0800")
    mod_data["I01"] = readin("VNP02IMG", 2015, 80, ["I01"], time="0800")["I01"]
    _create_bowtie_correct(
        mod_data["I01"], mod_data, "bowtie_correction_375m.nc", res=375
    )
    _create_mesh_correct(mod_data["I01"], mod_data, "mesh_correction_375m.nc", res=375)

    mod_data = readin("VNP03MOD", 2015, 80, ["longitude", "latitude"], time="0800")
    mod_data["M01"] = readin("VNP02MOD", 2015, 80, ["M01"], time="0800")["M01"]
    _create_bowtie_correct(
        mod_data["M01"], mod_data, "bowtie_correction_750m.nc", res=750
    )
    _create_mesh_correct(mod_data["M01"], mod_data, "mesh_correction_750m.nc", res=750)
