Installation
============

From conda
----------

Coming soon!


From PyPI
---------

Coming soon!


..
   Using pip
   ---------

   For the most up-to-date version of the code, you should use ``pip`` pointed at the git repository.


   .. code-block:: bash

       $ pip install git+git://github.com/edgrrr/csat2.git


From source code
----------------

If you are using and working on the code (or might want to update it from git), you can install an editable version

.. code-block:: bash

    $ cd path/to/csat2
    $ python setup.py develop


Note that the netCDF4 library should be built to read HDF files. This is the case if you install it from conda.

Configuration
-------------

csat2 comes with a simple configuration for the JASMIN service, run by CEDA/NERC. To use csat2 for a different machine, you will need to create a configuration file for your local machine.

When you import csat2, it searches the default configuration directory ``${HOME}/.csat2``. The configuration file ``config.cfg``, stored in this directory is used to determine which machine-file to use. An example config file will define a glob-able pattern which then links to a machine-file. Many patterns can point to the same machine file and they are searched in order until a match is found.

.. code-block:: yaml

   machines:
        "*.jasmin.ac.uk": jasmin.txt
        "*.ceda.ac.uk": jasmin.txt
        "*": local-machine.txt


Note that although only one config file can be specified, the config can fall-back to the default config, so only extra configuration needs to be supplied.


Machine files
-------------

Machine files contain the patterns used to locate files by the csat2 locator module and are designed such that code can be machine independent (e.g. code written on your laptop will run on JASMIN with no changes of file paths). These patterns are python string formatters. As such, they allow the relevant formatting commands when building a filename (e.g. zero-padding).

An example configuration file (only for MODIS files on JASMIN) looks like

.. code-block:: python

   # Comment
   M neodc_folder=/neodc

   [MODIS]
   -[MYD03|MOD03|MYD04_L2|MOD04_L2|MYD06_L2|MOD06_L2|MYD35_L2|MOD35_L2|MYD021KM|MOD021KM]
    {neodc_folder}/modis/data/{product}/collection{collection}/{year}/{mon:0>2}/{day:0>2}/*.{time}.*.hdf


Macros are strings that are immediately substituted. This is useful if you might change the root location of the data but keep the overall filetree structure the same. Note that ``product``, the second level of detail (e.g. MYD03, MOD03) is automatically a macro and so immediately substituted.

Each instrument is defined in a line starting with a ``[``. Following this are lines starting with ``-[``. These define individual products, with the file pattern defined on the following line. You can have several products for each instrument, either defined in multiple lines or in a single line (as above, separated by ``|``.

Path syntax is the same as the python string formatting syntax (what a coincidence...). Use ``{variable}`` to make a substitution - can be a macro or a variable from the readin function. ``{doy:0>3}`` left pads the variable doy with 0s to a width of three characters, ``*`` can be used as a wildcard (e.g. for MODIS files that include the processing date.

If ``(year, mon, day)`` or ``(year, doy)`` are supplied to the locator function, the corresponding one is also calculated.

The simplest machine file will store all files in a local folder. Here is an example I use for my laptop.

.. code-block:: python

   M data_folder=/home/edward/LocalData

   [MODIS]
   -[MYD06_L2|MOD06_L2|MYD021KM|MOD021KM]
    {data_folder}/MODIS/{product}.A{year}{doy:0>3}.{time}.*.hdf
   -[MOD08_D3|MYD08_D3]
    {data_folder}/MODIS/{product}.A{year}{doy:0>3}*
   -[bowtie]
    {data_folder}/MODIS/bowtie_correction_{res}_{length}.nc

   [ECMWF]
   -[ERA5]
    {data_folder}/ECMWF/{year}{time}/{variable}_{level}_{doy:0>3}.nc


