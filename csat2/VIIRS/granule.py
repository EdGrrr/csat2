from datetime import datetime, timedelta
from csat2 import misc
from csat2 import locator

# , field_interpolate, bowtie_correct
from .readfiles import readin, readin_metadata, DEFAULT_COLLECTION
from .util import band_res, bowtie_correct, mesh_correct
from .download import download, download_geometa
import numpy as np

granule_length_minutes = 6


class Granule(object):
    """An object for plotting data for a MODIS granule"""

    def __init__(self, year, doy, time, sat, col=DEFAULT_COLLECTION):
        """year, doy, time - all int
        sat = 'N' or 'J' (Suomi or JPSS-1)"""
        super(Granule, self).__init__()
        self.year = year
        self.doy = doy
        self.time = int(time)
        if sat not in ["N", "J"]:
            raise ValueError("Sat must be N or J")
        self.sat = sat
        self.locator = None
        self.lonlat = None
        self._orbit = None
        self.col = col

    @classmethod
    def fromtext(cls, gran_text, col=DEFAULT_COLLECTION):
        return cls(
            int(gran_text[:4]),
            int(gran_text[4:7]),
            int(gran_text[8:12]),
            gran_text[-1],
            col,
        )

    @classmethod
    def fromfilename(cls, filename):
        basename = filename.split("/")[-1]
        sat = {"N": "N", "J": "J"}[basename[1]]
        year = int(basename.split(".")[1][1:5])
        doy = int(basename.split(".")[1][5:8])
        time = int(basename.split(".")[2])
        return cls(year, doy, time, sat)

    def orbit(self):
        if not self._orbit:
            dat = readin("GEOMETA", self.year, self.doy, self.sat)
            self._orbit = dat["Orbit"][
                np.where(
                    dat["StartTime"] == self.datetime().strftime("%Y-%m-%d %H:%M")
                )[0]
            ][0]
        return self._orbit

    def timestr(self):
        return "{:0>4}".format(self.time)

    def astext(self):
        return "{}{:0>3}.{:0>4}{}".format(self.year, self.doy, self.time, self.sat)

    def datetime(self):
        _, mon, day = misc.time.doy_to_date(self.year, self.doy)
        return datetime(
            self.year, mon, day, int(self.time) // 100, int(self.time) % 100
        )

    def product_expand(self, product):
        if not product.startswith("V"):
            # Add on the correct satellite prefix
            product = "V{}{}".format({"N": "NP", "J": "J1"}[self.sat], product)
        return product

    def get_lonlat(self, product=None, col=None, dateline=False):
        """Get lon lat data - can specify the product or the collection
        used to the geolocation data from"""
        if not self.lonlat:
            self._read_lonlat(product=product, col=col, dateline=dateline)
        return self.lonlat

    def _read_lonlat(self, product=None, col=None, dateline=False):
        if not col:
            col = self.col
        if not product:
            product = "03MOD"
        data = readin(
            self.product_expand(product),
            self.year,
            self.doy,
            sds=["longitude", "latitude"],
            time=self.timestr(),
            col=col,
        )
        if dateline:
            # Ignores the poles
            if (data["longitude"].max() > 160) and (data["longitude"].min() < -160):
                data["longitude"] = np.where(
                    data["longitude"] < 0, data["longitude"] + 360, data["Longitude"]
                )
        self.lonlat = [data["longitude"], data["latitude"]]

    def get_band_radiance(
        self, band, col=None, refl=False, bowtie_corr=False, mesh_corr=False, nrt=False
    ):
        """
        --Bowtie Corrected-- => bowtie_corr=True

        Return the radiance for a specific band from a specific granule.
        Set 'refl' to True to get the reflectance instead (if supported)"""
        if col is None:
            col = self.col

        metadata = self.get_metadata_band(band, col=col, nrt=nrt)
        vir_data = readin(
            self.get_radiance_product(band, nrt=nrt),
            self.year,
            self.doy,
            [band],
            time=self.timestr(),
            col=col,
        )[band]
        if "radiance" in metadata["long_name"]:  # TIR bands
            if refl:
                raise ValueError("Reflectance cannot be computed for {}".format(band))
            else:
                # Already radiance, so continue
                pass
        else:  # A VIS/NIR band
            if refl:
                # Already reflectance, so continue
                pass
            else:
                vir_data /= metadata["radiance_scale_factor"]

        if bowtie_corr:
            return bowtie_correct(vir_data, res=band_res(band))
        elif mesh_corr:
            return mesh_correct(vir_data, res=band_res(band))
        else:
            return vir_data

    def get_band_bt(
        self, band, col=None, bowtie_corr=False, mesh_corr=False, nrt=False
    ):
        if col is None:
            col = self.col

        vir_data = readin(
            self.get_radiance_product(band, nrt=nrt),
            self.year,
            self.doy,
            [band, band + "_brightness_temperature_lut"],
            time=self.timestr(),
            col=col,
            scale=False,
        )
        if bowtie_corr:
            return bowtie_correct(
                vir_data[band + "_brightness_temperature_lut"][
                    vir_data[band].astype("int")
                ],
                res=band_res(band),
            )
        elif mesh_corr:
            return mesh_correct(
                vir_data[band + "_brightness_temperature_lut"][
                    vir_data[band].astype("int")
                ],
                res=band_res(band),
            )
        else:
            return vir_data[band + "_brightness_temperature_lut"][
                vir_data[band].astype("int")
            ]

    def download_product(self, product, col=None, nrt=False):
        if col is None:
            col = self.col

        if not self.check_product(product, col):
            download(
                self.product_expand(product),
                self.year,
                self.doy,
                times=[self.timestr()],
                col=col,
            )

    def download_geometa(self, col=None, force_redownload=False):
        if not col:
            col = self.col
        download_geometa(
            self.year,
            self.doy,
            {"N": "NPP", "J": "JPSS1"}[self.sat],
            col=col,
            force_redownload=force_redownload,
        )

    def get_filename(self, product, col=None, return_primary=True):
        if not col:
            col = self.col
        if product == "GEOMETA":
            fnames = locator.search(
                "VIIRS",
                "GEOMETA",
                year=self.year,
                doy=self.doy,
                sat={"N": "NPP", "J": "JPSS1"}[self.sat],
                collection=col,
            )
        else:
            fnames = locator.search(
                "VIIRS",
                self.product_expand(product),
                year=self.year,
                doy=self.doy,
                time=self.timestr(),
                collection=col,
            )
        if return_primary:
            return fnames[0]
        else:
            return fnames

    def check_product(self, product, col=None):
        if col is None:
            col = self.col

        filename = locator.search(
            "VIIRS",
            self.product_expand(product),
            year=self.year,
            doy=self.doy,
            collection=col,
            time=self.timestr(),
        )
        if len(filename) == 1:
            return 1
        else:
            return 0

    def get_radiance_product(self, band, nrt=False):
        prefix = {"N": "VNP", "J": "VJ1"}[self.sat]
        if nrt:
            suffix = "_NRT"
        else:
            suffix = ""
        if band.startswith("I"):
            return prefix + "02IMG" + suffix
        elif band.startswith("M"):
            return prefix + "02MOD" + suffix
        elif band.startswith("DNB"):
            return prefix + "02DNB" + suffix
        else:
            raise KeyError("Not a valid band")

    def get_metadata_band(self, band, col=None, nrt=False):
        if col is None:
            col = self.col

        metadata = readin_metadata(
            self.get_radiance_product(band, nrt=nrt),
            self.year,
            self.doy,
            self.timestr(),
            col,
            band,
        )
        return metadata

    def __repr__(self):
        return "Granule: " + self.astext()

    def increment(self, number=1):
        dt = self.datetime()
        dt += timedelta(minutes=granule_length_minutes * number)
        _, doy = misc.time.date_to_doy(dt.year, dt.month, dt.day)
        ntime = int("{:0>2}{:0>2}".format(dt.hour, dt.minute))
        return Granule(dt.year, doy, ntime, self.sat)
