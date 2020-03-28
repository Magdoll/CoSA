from setuptools import setup, Extension, find_packages
import sys

__author__ = "etseng@pacb.com"
version = "1.0.0"


setup(
    name = 'CoSA',
    version=version,
    description='Coronavirus Sequence Analysis',
    author='Elizabeth Tseng',
    author_email='etseng@pacb.com',
    zip_safe=False,
    packages = ['cosa'],
    package_dir = {'cosa': 'cosa'},
    install_requires=[
        'biopython',
        ],
    scripts = ['cosa/clean_up_metadata.py',
               'cosa/filter_gappedshort.py',
               ],
    )
