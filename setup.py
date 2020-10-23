#! /usr/bin/env python3

from setuptools import setup

setup(name = 'csat2',
      version = '0.1',
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
