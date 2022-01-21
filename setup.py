#!/usr/bin/env python
# Joao Couto 

import os
from os.path import join as pjoin
from setuptools import setup
from setuptools.command.install import install


longdescription = ''' Random tools for the neuropixels paris neuro course workshop. '''

setup(
    name = 'pnc_spks',
    version = '0.0b',
    author = 'Joao Couto',
    author_email = 'jpcouto@gmail.com',
    description = ('Neuropixels workshop.'),
    long_description = longdescription,
    license = 'GPL',
    packages = ['pnc_spks'],
    entry_points = {
      'console_scripts': [
        'plot-bin = pnc_spks.scripts.plot_binary:main',
        'concat-bin = pnc_spks.scripts.concatenate_files:main',
      ]
    }
)

