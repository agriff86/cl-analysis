#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="clanalysis",
    version="0.1.0",
    description="Detect PBL in Vaisala NetCDF files and generate plots",
    author="Alan Griffiths",
    author_email="alan.griffiths@ansto.gov.au",
    license="MIT",
    py_modules=["clanalysis"],
    entry_points={
        "console_scripts": ["clanalysis=clanalysis:main"],
    },
    install_requires=[
        "netCDF4>=1.2.9",
        "xarray>=0.12.1",
        "scipy>=1.3.1",
        "matplotlib>=3.1",
    ],
    keywords=["ceilometer", "cl31", "cl51", "boundary-layer", "lidar"],
    url="https://github.com/agriff86/clanalysis",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
)
