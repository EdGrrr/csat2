Usage
=====

Getting started
---------------

To ensure csat2 can find files on your system, setup a configuration and machine file for your system (see the introduction).


Data model
----------

The primary class of csat2 for each instrument is the ``Granule`` class - e.g. ``csat2.MODIS.Granule``. A granule represents a discrete 'chunk' of data. This is usually a chunk of data recorded at around the same time. The Granule object then store information about the time and location of the observations, as well as providing methods to access data and to geolocate observations. A Granule object will also provide a ``next()`` method, to support the processing of a series of observations.

A Granule could represent a single file, but it does not have to. For each granule, you can then have multiple *products*. These all are related to the same data 'chunk'.

As an example, MODIS data is organised into Granules that are each 5 minutes long. This is the subdivision used in many NASA products, particularly those that focus on the atmosphere. The different Granule products are then the different NASA products (e.g. MOD06_L2 - the cloud retrieval).

Downloading data
----------------

csat2 simplifies the accessing of data on various computers, but it also provides methods for downloading data.

For a granule object, a specified product is downloaded using the ``download()`` method. The data is then stored in a location according to the machine file, such that it can then be found and loaded in future.

The ``download()`` function has a counterpart ``check()`` function, that checks whether data is already available. By default, the ``download()`` function returns without downloading any files if the data already exists within the file structure defined in the relevant machine file. This behaviour can be overridden using the ``force_redownload`` argument for the ``download()`` method.

This requires some setup for each dataset used, such as applying for data access accounts. The details of this are provides in the instrument specific guides.


Troubleshooting
---------------

The module writes diagnostic output into the default logging location. You can change the logging level and print it to screen to aid with debugging.

.. code-block:: python

   >>> import logging
   >>> logging.basicConfig(level=logging.DEBUG)
