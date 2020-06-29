#!/usr/bin/env python

"""
Add the missing vars required to run CABLE

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (11.03.2020)"
__email__ = "mdekauwe@gmail.com"


import numpy as np
import xarray as xr
import os
import sys
import netCDF4
import shutil

def main(fname, out_fname):

    f = netCDF4.Dataset(fname, 'r+')

    f.setncattr("Qair", 'units': u"kg/kg")


    f.close()



if __name__ == "__main__":

    fname = "raw/01_preston_metforcing_1993_2004_UTC_v1.nc"
    out_fname = "01_preston_metforcing_1993_2004_UTC_v1.nc"
    main(fname, out_fname)
