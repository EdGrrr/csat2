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
   

The granule object can be used to download MODIS files if you have the location defined in your machine file and an API-key for the LAADS DAAC.

.. code-block:: python
                
   >>> gran.download('06_L2')  # Download cloud data
   >>> gran.download('03')  # Download the geolocation data


To get data for a granule, you can read a variable direct from a file, using the ``get_variable`` method, or get brightness temperature/reflectance data from the radiance product ``021KM`` (use ``refl=True``) to get reflectance. Set ``bowtie_corr=True`` to re-order the edge of swath pixels to account for the bowtie effect.

.. code-block:: python

   # Get the effective radius data from the cloud '06_L2' product
   >>> cer = gran.variable('06_L2', ['Cloud_Effective_Radius'])
   # Reflectivity from band 2 (0.865um)
   >>> refl2 = gran.get_band_radiance(2, refl=True, bowtie_corr=False)


Finally, you can step forward a specified number of granules

.. code-block:: python
                
   >>> newgran = gran.increment(number=1)


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


GOES data
---------

Coming soon, but the overall structure is very similar to the MODIS Granule structure. The download method uses google-cloud storage, so you will need to set that up to use it.
   
   
Troubleshooting
---------------

The module writes diagnostic output into the default logging location. You can change the logging level and print it to screen to aid with debugging.

.. code-block:: python

   >>> import logging
   >>> logging.basicConfig(level=logging.DEBUG)
