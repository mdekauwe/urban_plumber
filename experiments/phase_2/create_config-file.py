#!/usr/bin/env python


''' This file forms part of the Urban-PLUMBER benchmarking evaluation project
for urban areas. Copyright 2020. This program is distributed freely in the hope
that it will be useful, but with NO WARRANTY, implied or otherwise. Users are
responsible for ensuring it is fit for their purpose. '''

__title__   = 'Urban-PLUMBER model configuration script'
__version__ = 'v1.0 (2020-05-22)'
__author__  = 'Mathew Lipson'
__email__   = 'm.lipson@unsw.edu.au'

'''
This python3 script is an example template which converts site parameters into
a model configuration file, e.g. a fortran namelist. This example will produce
a configuration based on the requirements of model NOAH-LSM, participants need
to update the script to fit the requirements of their model. The values and
assumptions are indicative of the method only.

Run the script "as-is" to see output produced in ./XX-Test/output

Instructions:
    1. MAIN: Update number of configuration files required (c1, c2 etc)
    2. INPUT INFORMATION: Update paths to your files
    3. CONFIG FILES: Update the create_config function(s) to your requirements
    4. run script

'''

import numpy as np
import pandas as pd
import netCDF4 as nc

# add sites to list as required
#sitelist = ['AU-Preston','GR-HECKOR','PL-Lipowa','US-Minneapolis1',\
#            'AU-SurreyHills','JP-Yoyogi','PL-Narutowicza','US-Minneapolis2',\
#            'CA-Sunset','KR-Jungnang','SG-TelokKurau','US-WestPhoenix',\
#            'FI-Kumpula','KR-Ochang','UK-KingsCollege',\
#            'FI-Torni','MX-Escandon','UK-Swindon',\
#            'FR-Capitole','NL-Amsterdam','US-Baltimore']
sitelist = ['AU-Preston','PL-Lipowa','US-Minneapolis1',\
            'AU-SurreyHills','JP-Yoyogi','PL-Narutowicza','US-Minneapolis2',\
            'CA-Sunset','KR-Jungnang','SG-TelokKurau','US-WestPhoenix',\
            'FI-Kumpula','KR-Ochang','UK-KingsCollege',\
            'FI-Torni','MX-Escandon','UK-Swindon',\
            'FR-Capitole','US-Baltimore']

###############################################################################
##### 1. MAIN
#####    Update number of configuration files as required (c1, c2 etc)
###############################################################################

def main(sitename):
    print('configuring %s' %sitename)

    print('loading site data and forcing information')
    info = set_info(sitename)
    sitedata,info = set_more_info(info)

    print('creating model configuration files')

    # config 1
    c1 = create_config1(sitedata,info)
    # write config file
    f1  = "all_sites/%s" %(info['fname_config1'])
    with open(file=f1, mode='w') as ofile:
        ofile.write(c1)


    print('done! see %s for new files' %info['outpath'])

    return

###############################################################################
##### 2. INPUT INFORMATION:
#####    Input path information.
###############################################################################

def set_info(sitename):
    '''This function sets basic path input information.

    Parameters
    ----------
    sitename (string) : the sitename (set in __main__ scope at end of file)

    Outputs
    -------
    info (dictionary) : information for use by other functions
    '''

    info = {}

    # defines sitename from sitelist
    info['sitename'] = sitename

    # path to site information
    info['sitepath'] = 'all_sites'

    # path to model configuration (output)
    info['outpath'] = './%s/namelists' %sitename

    # name of input forcing file

    if sitename != "AU-Preston":
        info['fname_forcing'] = 'met/%s_metforcing_v1.nc' %sitename
    else:
        info['fname_forcing'] = 'met/%s_metforcing_v3.nc' %sitename

    #info['fname_forcing'] = 'met/%s_metforcing_v1.nc' %sitename

    # name of configuration output files (add or remove files as needed)
    info['fname_config1'] = 'cable_%s_up.nml' %sitename

    return info

###############################################################################
##### 3. CONFIG FILES
#####    Update the create_config function(s)
#####    a) Replace "template" string with your model's default configuration.
#####    b) Define variables that will change from the default configuration.
#####    c) Use those variables to replace values from in the template string.
#####        Avoid hardcoding value changes in the default configuration
#####        template string, define new variables if required so changes can
#####        be tracked. Here python3.6 'f-string' functionality is utilised,
#####        see: https://www.python.org/dev/peps/pep-0498/
###############################################################################

def create_config1(sitedata,info):
    '''create string for model configuration file

    Parameters
    ----------
    sitedata (dataframe):  site configuration information
    info (dictionary):     path and timing information

    Output
    ------
    template (str): full string of namelist configuration'''

    ###############################################################################
    #### define any parameters that change the default configuration

    startdate   = pd.to_datetime(info['time_coverage_start']).strftime('%Y%m%d%H%M')
    enddate     = pd.to_datetime(info['time_coverage_end']).strftime('%Y%m%d%H%M')
    loops       = 0  # no looping of spinup data
    outdir      = info['outpath']
    latitude    = sitedata['latitude']
    longitude   = sitedata['longitude']
    timestep    = info['timestep_interval_seconds']
    meas_height = sitedata['measurement_height_above_ground']
    usemonalb   = '.TRUE.'
    alb         = sitedata['average_albedo_at_midday']

    ###############################################################################
    ##### CONFIGURATION TEMPLATE
    ##### REPLACE THE "template" STRING WITH YOUR MODEL'S CONFIGURATION
    ##### Default source:
    ###############################################################################

    template= f'''\
&METADATA_NAMELIST
 startdate             = {startdate}
 enddate               = {enddate}
 loop_for_a_while      = {loops}
 output_dir            = "{outdir}"
 Latitude              = {latitude}
 Longitude             = {longitude}
 Forcing_Timestep      = {timestep}
 CABLE_Timestep      = 900
 Sea_ice_point         = .FALSE.
 Soil_layer_thickness  =    0.1            0.3            0.6            1.0
 Soil_Temperature      =  266.0995       274.0445       276.8954       279.9152
 Soil_Moisture         =  0.2981597      0.2940254      0.2713114      0.3070948
 Soil_Liquid           =  0.1611681      0.2633106      0.2713114      0.3070948
 Skin_Temperature      =  263.6909
 Canopy_water          =  3.9353027E-04
 Snow_depth            =  1.0600531E-03
 Snow_equivalent       =  2.0956997E-04
 Deep_Soil_Temperature = 285
 Landuse_dataset       = "USGS"
 Soil_type_index       = 8
 Vegetation_type_index = 7
 Urban_veg_category    = 1
 glacial_veg_category  = 24
 Slope_type_index      = 1
 Max_snow_albedo       = 0.75
 Air_temperature_level = {meas_height}
 Wind_level            = {meas_height}
 Green_Vegetation_Min  = 0.01
 Green_Vegetation_Max  = 0.96
 Usemonalb             = {usemonalb}
 Rdlai2d               = .FALSE.
 sfcdif_option         = 1
 iz0tlnd               = 0
 Albedo_monthly        = {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb},  {alb}
 Shdfac_monthly        = 0.01,  0.02,   0.07,  0.17,  0.27,  0.58,  0.93,  0.96,  0.65,  0.24,  0.11,  0.02
 lai_monthly           = 4.00,  4.00,   4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00,  4.00
 Z0brd_monthly         = 0.020, 0.020, 0.025, 0.030, 0.035, 0.036, 0.035, 0.030, 0.027, 0.025, 0.020, 0.020
! Z0brd_monthly         = 0.100, 0.140, 0.180, 0.220, 0.250, 0.280, 0.250, 0.220, 0.180, 0.140, 0.120, 0.100
 Use_urban_module      = .TRUE.
/


'''
    return template


###############################################################################
##### 4. OTHER FUNCTIONS
###############################################################################

def set_more_info(info):
    '''sets timing and other info from forcing file and places into dictionary.

    Inputs
    ------
    info (dictionary): script information

    Outputs
    -------
    info (dictionary): updated script information
    sitedata (dataframe): site information
    '''

    # load site data table
    if sitename != "AU-Preston":
        path_sitedata = '%s/site_info/%s_sitedata_v1.csv' %(info['sitepath'], info['sitename'])
    else:
        path_sitedata = '%s/site_info/%s_sitedata_v3.csv' %(info['sitepath'], info['sitename'])
    sitedata_full = pd.read_csv(path_sitedata, index_col=1, delimiter=',')
    sitedata      = pd.to_numeric(sitedata_full['value'])

    # load forcing information
    fpath = '%s/%s' %(info['sitepath'],info['fname_forcing'])

    print(fpath)
    print(info['sitepath'])
    print(info['fname_forcing'])
    with nc.Dataset(filename=fpath, mode='r', format='NETCDF4') as f:

        info['time_coverage_start']       = f.time_coverage_start
        info['time_coverage_end']         = f.time_coverage_end
        info['time_analysis_start']       = f.time_analysis_start
        info['local_utc_offset_hours']    = f.local_utc_offset_hours
        info['timestep_interval_seconds'] = f.timestep_interval_seconds
        info['timestep_number_spinup']    = f.timestep_number_spinup
        info['timestep_number_analysis']  = f.timestep_number_analysis

    return sitedata, info

###############################################################################
##### __main__  scope
###############################################################################

if __name__ == "__main__":

    for sitename in sitelist:

        main(sitename)

###############################################################################
