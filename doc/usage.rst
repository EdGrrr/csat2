Usage
=====

Getting started
---------------

One you have setup a configuration file for you local computer, you can get started using csat2.

MODIS files
-----------

The main file for reading MODIS data is ``csat2.MODIS.readin``, which will read a variety of different MODIS data formats, from the L3 daily data to swath L2 data and the subset and cdnc_best files at Imperial. This function will return an xarray with the requested data.

.. code-block:: python

   import csat2.MODIS
   data = csat2.MODIS.readin('cdnc_best', year=2018, doy=100, sds=['CDNC_best', 'LWP_all'])


This provides location-independent access to MODIS files, assuming machine files have been setup on the required computers. However, the csat2 library also provides the ``Granule`` interface, for dealing primarily with L2 (swath) data.

MODIS Granule interface
-----------------------

A MODIS Granule is a single 5 minute section of data collected by the instrument. It is the sub-division used in many of the NASA derived products, particularly those that focus on the atmosphere. The ``Granule`` class is designed to make using these files as simple as possible.

A ``Granule`` object contains information about the granule, allows a step forwards in time and links to a ``MODISLocator`` object, that aids in geolocation of MODIS pixels. A granule can be constructed explicitly, from a filename or a granule shortname.

.. code-block:: python

   >>> import csat2.MODIS
   >>> gran = csat2.MODIS.Granule(year=2015, doy=11, time='2110', sat='A')
   >>> gran = csat2.MODIS.Granule.fromtext('2015011.2110A')
   >>> print(gran)
   Granule: 2015011.2110A

You can then return information about the granule - note that this requires the GEOMETA data from LAADS to be downloaded.

.. code-block:: python

   >>> gran.datetime()
   datetime.datetime(2015, 1, 11, 21, 10)
   >>> gran.orbit()
   67507
   >>> gran.daynight()
   'D'


The granule object can also provide the approximate nearest pixel to a list of lon-lat points (points are always given in a ``[lon, lat]`` format. Sun and observation zenith angles can also be calculated, along with the interpolated lon-lat locations for each pixel.

.. code-block:: python

   >>> gran.locate([[-130, 31]])
   array([[ 774, 1322]])
   >>> sunz, suna, satz, sata, azidiff = gran.get_angles(mst=False)
   >>> lon, lat = gran.get_lonlat()


To get data for a granule, you can read a variable direct from a file, using the ``get_variable`` method, or get brightness temperature/reflectance data from the radiance product ``021KM`` (use ``refl=True``) to get reflectance. Set ``bowtie_corr=True`` to re-order the edge of swath pixels to account for the bowtie effect.

.. code-block:: python

   # Get the effective radius data from the cloud '06_L2' product
   >>> cer = gran.get_variable('06_L2', ['Cloud_Effective_Radius'])
   # Reflectivity from band 2 (0.865um)
   >>> refl2 = gran.get_band_radiance(2, refl=True, bowtie_corr=False)


Finally, you can step forward a specified number of granules

.. code-block:: python

   >>> newgran = gran.increment(number=1)

   
Downloading MODIS data
^^^^^^^^^^^^^^^^^^^^^^

The granule object can be used to download MODIS files if you place a NASA Earthdata username and password in the csat2 configuration directory (``${HOME}/.csat2/earthdata_auth.json``). This should be in json format, such that the contents of the file read something like

.. code-block:: json

   {
       "username": "<YOUR_USERNAME>",
       "password": "<YOUR_PASSWORD"
   }
   

This replaces the old LAADS token API. You can then download granule data as follows. This will also work for VIIRS data.

.. code-block:: python

   >>> gran.download('06_L2')  # Download cloud data
   >>> gran.download('03')  # Download the geolocation data


   

ECMWF files
-----------

The ECMWF files are read from a set of pre-processed files, stored in a one day per file, one level per file format. Code to create these files will be included in csat2 shortly (it is currently on seldon).

There are two ways to access ECMWF data. The ``ECMWF.readin_ERA`` function (accessed through ``ECMWF.readin``) returns an xarray with the data for the requested data and time. It is also able to calculate the LTS and EIS, assuming appropriate input data

.. code-block:: python

   >>> import csat2.ECMWF
   >>> eis = ECMWF.readin('ERA5', 2015, 11, 'EIS')
   >>> t1000_1330LST = ECMWF.readin('ERA5', 2015, 11, 'Temperature', level='1000hPa', time='LST')[2]


The second method is through the ``ECMWF.ERA5Data`` object. This stores the netcdf data, allowing faster access to data that is already in use. It is designed to provide access to a single variable and level. ``ECMWF.ERA5WindData`` provides access to the U and V wind components together.

.. code-block:: python

   >>> temp_data = ECMWF.ERA5Data('Tenperature', level='1000hPa', res='1grid')
   >>> t1000 = temp_data.get_data([100, 101, 102], [10, 9 ,8], datetime.datetime(2015, 1, 1, 10))


ERA5 data can be downloaded using the ``ECMWF.download.download`` function. This will place the ERA5 data in the location specified in the machine file, as well as calculating the local solar time files (if required). To keep the requests manageable, it will only request one month and one level at a time, but multiple variables can be requested on the same level. This requires a Copernicus data store key (`see this guide from ECMWF <https://cds.climate.copernicus.eu/api-how-to>`_)

.. code-block:: python

   >>> ECMWF.download(2020, 1, ['Temperature', 'Relative_humidity', 'U-wind-component'], level='1000hPa', resolution='1grid')

Note that this requires `CDO <https://code.mpimet.mpg.de/projects/cdo/>`_ to be installed (which you can do through anaconda), as it uses it for splitting up the netcdf files (and regridding where required).

The variable names here are the local names, which mostly (but not always) match the Copernicus names. Windspeed and SST are the main exceptions.



GOES data
---------

The GOES data uses a similar Granule structure to MODIS. You create a granule either by defining all the relevant time properties, from a filename or from a text granule name.

.. code-block:: python

   >>> from csat2 import GOES
   >>> gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
   >>> gran = GOES.Granule.fromfilename('OR_ABI-L1b-RadC-M3C01_G16_s20180020002199_e20180020004572_c20180020005016.nc')

As with MODIS, you can then read in scientific data from the granule object. The scan mode can be specified if desired, but as there is only one scan mode at any time, the default is to ignore it.

.. code-block:: python

   >>> gran.get_band_bt(channel=13)
   <xarray.DataArray (y: 1500, x: 2500)>
   array([[       nan,        nan,        nan, ..., 237.863262, 238.101511, 238.810748],
          [       nan,        nan,        nan, ..., 236.718719, 236.53601 , 238.160925],
          [       nan,        nan,        nan, ..., 235.8613  , 235.675966, 238.042034],
          ...,
          [292.08001 , 292.229971, 292.199996, ..., 294.306864, 294.336236, 294.336236],
          [292.199996, 292.469458, 292.439552, ..., 294.189298, 294.277484, 294.306864],
          [292.319846, 292.559118, 292.588989, ..., 293.953766, 294.101032, 294.218703]])
   Coordinates:
       t        datetime64[ns] 2018-01-02T00:03:39.182033984
       y_image  float32 0.08624
       x_image  float32 -0.03136
     * y        (y) float32 0.128212 0.128156 ... 0.044324003 0.044268005
     * x        (x) float32 -0.101332 -0.101276 ... 0.038556002 0.038612

You can get information about the granule, or step forwards in time. If you have the files on disk, you can also use the granule object to locate the relevant file (for a given channel). Note that the exact timing of the granule depends on the GOES scan pattern, so the search here is done only for a granule within the current granule increment time.

.. code-block:: python

   >>> gran.datetime()
   datetime.datetime(2018, 1, 2, 0, 0)
   >>> gran.next()
   G16.2018002.0005.RadC

The Granule object also allows you to geolocate pixels, or to locate a pixel given a lon/lat array. A channel is currently required, due to the varying resolution of the instrument. You should *not* switch resolutions for the same granule object and the locator is cached.

The locate function returns integers (for use an indices), unless you ask for a float using the option ``interp=True``. If the requested pixel is outside the current granule grid, a large negative is returned (or a ``np.nan`` if using ``interp=True``).

.. code-block:: python

   >>> gran.geolocate(np.array([[ 190, 840]]), channel=13)
   array([[-113.11428559,   28.98774504]])

   >>> gran.locate(np.array([[-113.114, 28.9877]]))
   array([[190, 840]])


Downloading GOES data
^^^^^^^^^^^^^^^^^^^^^

Setting this up is more complicated that for MODIS. The granule object is currently using Google cloud storage, which although public requires an API-key and project setup to use.

Start by logging into the `Google API console <https://console.developers.google.com/>`_. Create a new project, not linked to any organisation.

Once you have created a project, create a service account. You will then need to add a key to this account. When the dialogue box opens, pick ``json``. Download this key to your csat2 configuration folder as ``goes-service-key.json``. You should set the permissions on this file so that only you can read it.

This should then allow you to use the granule download functions for GOES-16 and GOES-17 data from the Google Cloud.

.. code-block:: python

   >>> gran.download(channel=13)


Troubleshooting
---------------

The module writes diagnostic output into the default logging location. You can change the logging level and print it to screen to aid with debugging.

.. code-block:: python

   >>> import logging
   >>> logging.basicConfig(level=logging.DEBUG)
