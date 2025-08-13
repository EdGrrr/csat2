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

    git clone https://github.com/EdGrrr/csat2.git
    pip install -e csat2

    
Setup
-----

One of the key aims of csat2 is to make code using observational data machine-independent. For that, you need to tell it where to find your files (and where to save anything you download). csat2 does this using two types of configuration file.

   1. Config (``.cfg``) files. These are used by csat2 to determine which system it is running on and hence which *machine file* to use.
   2. Machine files. These store the file location patterns for a specific computer (or set of computers). 

Example files are supplied in the ``csat2/config`` folder and are described below. These defaults are overriden by any local configuration supplied in the ``${HOME}/.csat2`` directory.


Configuration files
...................
     
csat2 comes with a simple configuration for the JASMIN service, run by CEDA/NERC. Along with example configurations for our local computers and an example one that includes all the datasets that csat2 is currently built to use. To use csat2 for a different machine, you will need to create a configuration file for your local machine.

The default directory for user configuration is ``${HOME}/.csat2``. The configuration file ``config.cfg``, stored in this directory is used to determine which system csat2 is running on (based on the fully-qualified domain name) and so which *machine file* (which contains the actual filepaths) to use.

The config file is a `YAML <https://yaml.org/>`_ file that contains a list of glob-able patterns. The machine fully-qualified domain name is compared to these in turn until a match is found. If no match is found in the user config, csat2 will fallback by default to the example config supplied in the ``csat2/config`` folder.

The matched name then links to a *machine file*. Many patterns can point to the same machine file.  For example, this config file snippet uses two domain name patterns to identify a computer that should use the 'jasmin.txt' machine file. If these matches fail, but the computer name is *my-local-computer*, then the local-machine.txt file will be used. Otherwise the supplied config list will be checked.

.. code-block:: yaml

   machines:
        "*.jasmin.ac.uk": jasmin.txt
        "*.ceda.ac.uk": jasmin.txt
        "my-local-computer": local-machine.txt

Note that although only one config file can be specified, the config can fallback to the default config (in the 'csat2/config' folder), so only extra configuration needs to be supplied. If not matches are found, csat2 raises an IOError.

If you want to use csat2 without any of the data loading functions, you can supply a minimal ``config.cfg`` file in your ``${HOME}/.csat2`` folder.

.. code-block:: yaml

   machines:
        "*": local-machine.txt


Which will trigger on any machine, but is unlikely to contain the correct file paths for your system.


Machine files
-------------

Machine files contain the patterns used to locate files by the csat2 locator module. These patterns are python string formatters. As such, they allow the relevant formatting commands when building a filename (e.g. zero-padding in date strings).

An example configuration file (only for MODIS files on JASMIN) looks like

.. code-block:: python

   # Comment
   M neodc_folder=/neodc

   [MODIS]
   -[MYD03|MOD03|MYD04_L2|MOD04_L2|MYD06_L2|MOD06_L2|MYD35_L2|MOD35_L2|MYD021KM|MOD021KM]
    {neodc_folder}/modis/data/{product}/collection{collection}/{year}/{mon:0>2}/{day:0>2}/*.{time}.*.hdf

Lines starting with ``M`` are macros. These are strings that are immediately substituted when the file is loaded. This is useful if you might change the root location of the data but keep the overall filetree structure the same. Note that ``product``, the second level of detail (e.g. MYD03, MOD03) is automatically a macro and so immediately substituted.

A line starting with a ``[`` defines an instrument. Each instrument contains one of more ``products``.

Lines starting with ``-[`` define individual products. You can have several products for each instrument, either defined in multiple lines or in a single line (as above, separated by ``|``).

Following each product definition line are one or more lines that contain file paths. Path syntax is the same as the python string formatting syntax (what a coincidence...). Use ``{variable}`` to make a substitution - this can be a macro (supplied in this file) or a variable from the readin function.  For example, ``{doy:0>3}`` left pads the variable doy with 0s to a width of three characters, ``*`` can be used as a wildcard (e.g. for MODIS files that include the processing date).

If ``(year, mon, day)`` or ``(year, doy)`` are supplied to the locator function, the corresponding one is also calculated. csat2 favours ``(year, doy)`` combination, with the January 1st being DOY 1.

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


