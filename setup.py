#! /usr/bin/env python3

from setuptools import setup
import re

# Version code from https://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
VERSIONFILE = "csat2/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setup(name = 'csat2',
      version = verstr,
      author='Edward Gryspeerdt',
      author_email='e.gryspeerdt@imperial.ac.uk',
      maintainer='Edward Gryspeerdt',
      maintainer_email='e.gryspeerdt@imperial.ac.uk',
      description='Python library for satellite and meteorological data ',
      license='Not currently open source',
      install_requires=['xarray',
                        'numpy',
                        'scipy',
                        'netCDF4',
                        'matplotlib',
                        'scikit-learn',
                        'google-cloud-storage'],
      packages = ['csat2'],
      package_dir = {'csat2': 'csat2'},
      package_data = {'csat2': ['data/*']}
)
