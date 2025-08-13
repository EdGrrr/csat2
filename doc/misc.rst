Misc
====

csat2 contains a large number of miscellaneous or helper functions that are useful for reading or working with satellite data.

astro
-----

These functions deal with calculating the solar zenith angle, along with related tasks, such as the calculation of sunrise/sunset times and insolation.


bitops
------

Functions for dealing with bitmask variables (such as quality flags).


dlist
-----

Reading in satellite data often involves dealing with dictionaries of lists. The functions support that usage pattern. See the examples for sample usage patterns.


fileops
-------

A set of functions for dumping numpy arrays to netcdf files, along with reading some rarer data foormats (such as ICARTT).


geo
---

Functions for working on a sphere, including distance and bearing calculations, coordinate rotation and great circle path calculation.


hdf
---

If the netcdf library you are using has been build correct, it can read hdf4 files (such as MODIS data) directly. This is the case if you install csat2 using conda, but not if you install it via pip. The hdf library provides a class that hides this distinction, providing a consistent interface to hdf files whether you are using the netCDF4 package or the pyhdf one.


plotting
--------

Functions relating to plotting. Includes some subplot labelling and colourbar functions.


stats
-----

Statistical functions, largely for dealing with dataset that including missing data. This includes nanmean, nanmax (which are now included in numpy) as well as linear regressions that can account for missing data.

stats.py also includes functions (``csat2.misc.stats.regression_gridded``) for calculating linear regressions and Pearson correlation coefficients along specific axes between two datasets.

The ``csat2.misc.stats.DataAccumulator`` class is also defined here. This class is designed to support the processing of very large datasets, storing only the summary variables. This class can easily calculate means across large datasets, but can also compute other statistics and a space-efficient manner, including linear regressions and histograms.


time
----

Functions to convert between different time systems commonly used in satellite datasets, including 'year, month, day', 'year, doy' and TAI/International atomic time.

This module also includes a set of functions for converting between UTC and local solar time, as well as for calculating time offsets to the previous/next satellite overpass for a given location.
