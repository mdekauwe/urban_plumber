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
import xarray as xr

def main(fname):

    f = netCDF4.Dataset(fname, 'r+')

    nc_attrs = f.ncattrs()
    nc_dims = [dim for dim in f.dimensions]
    nc_vars = [var for var in f.variables]

    ndim = 1
    f.createDimension('z', ndim)
    f.createDimension('y', ndim)
    f.createDimension('x', ndim)

    latitude = f.createVariable('latitude', 'f8', ('y', 'x',))
    latitude.units = "degrees_north"
    latitude.missing_value = -9999.
    latitude.long_name = "Latitude"

    longitude = f.createVariable('longitude', 'f8', ('y', 'x',))
    longitude.units = "degrees_east"
    longitude.missing_value = -9999.
    longitude.long_name = "Longitude"


    z = f.createVariable('z', 'f4', ('z',))
    z.units = " "

    y = f.createVariable('y', 'f4', ('y',))
    y.units = " "

    x = f.createVariable('x', 'f4', ('x',))
    x.units = " "

    x[:] = ndim
    y[:] = ndim
    z[:] = ndim
    latitude[:] = -37.73  # Coutts et al 2007a
    longitude[:] = 145.01 # Coutts et al 2007a

    f.close()



if __name__ == "__main__":

    in_fname = "raw/01_preston_metforcing_1993_2004_UTC_v1.nc"
    out_fname = "01_preston_metforcing_1993_2004_UTC_v1.nc"
    shutil.copyfile(in_fname, out_fname)

    main(out_fname)
