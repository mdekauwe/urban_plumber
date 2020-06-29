#!/usr/bin/env python


''' This file forms part of the Urban-PLUMBER benchmarking evaluation project
for urban areas. Copyright 2020. This program is distributed freely in the hope
that it will be useful, but with NO WARRANTY, implied or otherwise. Users are
responsible for ensuring it is fit for their purpose. '''

__title__   = 'Urban-PLUMBER netCDF creation script'
__version__ = 'v1.04 (2020-06-19)'
__author__  = 'Mathew Lipson'
__email__   = 'm.lipson@unsw.edu.au'

'''
This python3 script constructs a netCDF file which is formatted to comply
with the Urban-PLUMBER protocol and fills it with output model data.
Run the script "as-is" to create a netcdf from example model output.

Instructions:
    1. INPUT INFORMATION: update paths to your files and soil layer number
    2. GET OUTPUT: import your model output and change variables as required.
    3. SET OUTPUT: set model output into netcdf file.
    4. run script

Notes:
    standard_name attributes are based on CF conventions

Acknowledgments:
    With thanks Aristofanis Tsiringakis

Changelog:
    v1.01: Correct RootMoist attributes
    v1.02: add x,y dimensions to all variables to align with ALMA
    v1.03: change time datatype from i8 (64bit) to i4 (32bit), add history
    v1.04: added x,y dimension variable definitions
'''

import numpy as np
import pandas as pd
import netCDF4 as nc
import sys

missing_float = -9999.

# add sites to list as required
sitelist = ['AU-Preston']

###############################################################################
##### MAIN: no action required
###############################################################################

def main(sitename):

    print('loading site data and forcing information')
    info = set_info(sitename)
    sitedata,info = set_more_info(info)

    print('building empty netcdf in complying form')
    create_empty_netcdf(info)

    data = get_model_data(info)

    print('setting netcdf with output data')
    set_netcdf_data(data,info)

    print('done! see %s' %info['post_fname'])

    return

###############################################################################
##### 1. INPUT INFORMATION: update paths to your files and soil layer number
###############################################################################

def set_info(sitename):
    '''sets user inputs'''
    print('processing %s' %sitename)

    info = {}

    # site name
    info['sitename'] = '%s' %sitename

    # path to site data
    info['sitepath'] = './%s' %sitename

    # path to netcdf output
    info['outpath']  = './%s/outputs' %sitename

    # name of netcdf forcing file for site
    info['fname_forcing'] = 'met/%s_metforcing_v1.nc' %sitename

    # name of netcdf output file:
    ofn = "%s_out.nc" % (sitename)
    info['fname_output'] = '%s/%s' % (info['outpath'], ofn)

    # list number of soil layers in model
    info['num_soil_layers']  = 6

    # name of post process file
    info['post_fname'] = '%s/%s' %(info['outpath'],"%s_post_processed.nc" % (info['sitename']))

    return info

###############################################################################
##### 2. GET OUTPUT: import your model output and change variables as required.
###############################################################################

def get_model_data(info):
    '''
    Reads model output (in this example text) into a pandas dataframe.
    Update the function to match your model output data types.

    Inputs
    ------
    info (dictionary): script information

    Outputs
    -------
    data (dataframe): model output
    '''

    f = nc.Dataset(info['fname_output'])
    time = nc.num2date(f.variables['time'][:],
                        f.variables['time'].units)

    df = pd.DataFrame(f.variables['SWnet'][:,0,0], columns=['SWnet'])
    df['LWnet'] = f.variables['LWnet'][:,0,0]
    df['Qle'] = f.variables['Qle'][:,0,0]
    df['Qh'] = f.variables['Qh'][:,0,0]
    df['Qg'] = f.variables['Qg'][:,0,0]
    df['Snowf'] = f.variables['Snowf'][:,0,0]
    df['Rainf'] = f.variables['Rainf'][:,0,0]
    df['Evap'] = f.variables['Evap'][:,0,0]
    df['Qs'] = f.variables['Qs'][:,0,0]
    df['Qsb'] = f.variables['Qsb'][:,0,0]
    df['Qsm'] = f.variables['SnowMelt'][:,0,0]
    df['SnowT'] = f.variables['SnowT'][:,0,0]
    df['VegT'] = f.variables['CanT'][:,0,0]
    df['BaresoilT'] = f.variables['BaresoilT'][:,0,0]
    df['RadT'] = f.variables['RadT'][:,0,0]
    df['Albedo'] = f.variables['Albedo'][:,0,0]
    df['SWE'] = f.variables['SWE'][:,0,0]
    df['LAI'] = f.variables['LAI'][:,0,0]

    df['SoilMoist1'] = f.variables['SoilMoist'][:,0,0,0]
    df['SoilMoist2'] = f.variables['SoilMoist'][:,1,0,0]
    df['SoilMoist3'] = f.variables['SoilMoist'][:,2,0,0]
    df['SoilMoist4'] = f.variables['SoilMoist'][:,3,0,0]
    df['SoilMoist5'] = f.variables['SoilMoist'][:,4,0,0]
    df['SoilMoist6'] = f.variables['SoilMoist'][:,5,0,0]

    df['SoilTemp1'] = f.variables['SoilTemp'][:,0,0,0]
    df['SoilTemp2'] = f.variables['SoilTemp'][:,1,0,0]
    df['SoilTemp3'] = f.variables['SoilTemp'][:,2,0,0]
    df['SoilTemp4'] = f.variables['SoilTemp'][:,3,0,0]
    df['SoilTemp5'] = f.variables['SoilTemp'][:,4,0,0]
    df['SoilTemp6'] = f.variables['SoilTemp'][:,5,0,0]

    df['TVeg'] = f.variables['TVeg'][:,0,0]
    df['ESoil'] = f.variables['ESoil'][:,0,0]

    zse = [.022, .058, .154, .409, 1.085, 2.872]
    frac1 = zse[0] / (zse[0] + zse[1])
    frac2 = zse[1] / (zse[0] + zse[1])
    frac3 = zse[2] / (zse[2] + zse[3])
    frac4 = zse[3] / (zse[2] + zse[3])
    frac5 = zse[4] / (zse[4] + zse[4])
    frac6 = zse[5] / (zse[5] + zse[5])
    df['RootMoist'] = (f.variables['SoilMoist'][:,0,0,0] * frac1) +  \
                      (f.variables['SoilMoist'][:,1,0,0] * frac2) + \
                      (f.variables['SoilMoist'][:,2,0,0] * frac3) + \
                      (f.variables['SoilMoist'][:,3,0,0] * frac4) + \
                      (f.variables['SoilMoist'][:,4,0,0] * frac4) + \
                      (f.variables['SoilMoist'][:,5,0,0] * frac5)


    df['SWdown'] = f.variables['SWdown'][:,0,0]
    df['LWdown'] = f.variables['LWdown'][:,0,0]
    df['Tair'] = f.variables['Tair'][:,0,0]
    df['Qair'] = f.variables['Qair'][:,0,0]
    df['PSurf'] = f.variables['PSurf'][:,0,0]
    df['Wind'] = f.variables['Wind'][:,0,0]

    #df['dates'] = time
    #df = df.set_index('dates')

    # set timesteps using forcing information
    times = pd.date_range( start = info['time_coverage_start'],
                           end   = info['time_coverage_end'],
                           freq  = '%sS' %(info['timestep_interval_seconds']))
    df.index = times  # define time index in UTC

    return df

###############################################################################
##### 3. SET OUTPUT: set model output into netcdf file.
###############################################################################

def set_netcdf_data(data,info):
    '''
    Write model data values into existing netCDF file.
    This function sets a small number of variables as an example.
    Please include all available model output.
    Note subsurface state variables are two dimensional {time, soil_layer}.

    Inputs
    ------
    data (dataframe): model data
    info (dictionary): script information
    '''

    timesteps = info['timestep_number_spinup'] + info['timestep_number_analysis']

    # define empty numpy arrays for missing data
    no_data_1D = np.full([ timesteps ], missing_float)
    no_data_2D = np.full([ timesteps, info['num_soil_layers'] ], missing_float)

    # open netcdf files (r = read only, r+ = append existing)

    with nc.Dataset(filename=info['post_fname'], mode='r+', format='NETCDF4') as o:

        # set metadata
        o.title             = 'CABLE output for the Urban-PLUMBER project'
        o.site              = '%s' % (info['sitename'])
        o.experiment        = 'Baseline'
        o.institution       = 'UNSW'
        o.primary_contact   = 'Martin De Kauwe (mdekauwe@gmail.com)'
        o.secondary_contact = 'Martin De Kauwe (mdekauwe@gmail.com)'
        o.model             = 'CABLE'
        o.source            = 'The Community Atmosphere–Biosphere Land Exchange model, revision: 7287'
        o.references        = 'Kowalczyk EA, Wang YP, Wang P, Law RH, Davies HL. 2006. The csiro atmosphere biosphere land exchange (CABLE) model for use in climate models and as an offline model. CSIRO.; Wang YP, Kowalczyk E, Leuning R, Abramowitz G, Raupach MR, Pak B, Gorsel E van, Luhar A. 2011. Diagnosing errors in a land surface model (CABLE) in the time and frequency domains. Journal of Geophysical Research: Biogeosciences (2005–2012) 116.'
        o.repository        = 'https://trac.nci.org.au/trac/cable/browser/branches/Users/mgk576/trunk_plumber'
        o.site_experience   = 'No'
        o.additional_data   = 'No'
        o.comment           = 'No'
        o.history           = 'Created with %s at %s' %(__file__,pd.Timestamp.now())

        # Critical energy balance components
        o.variables['SWnet'][:]        = data['SWnet'].values       # Net shortwave radiation (downward)
        o.variables['LWnet'][:]        = data['LWnet'].values       # Net longwave radiation (downward)
        o.variables['Qle'][:]          = data['Qle'].values        # Latent heat flux (upward)
        o.variables['Qh'][:]           = data['Qh'].values        # Sensible heat flux (upward)
        o.variables['Qanth'][:]        = no_data_1D       # Anthropogenic heat flux (upward)
        o.variables['Qstor'][:]        = no_data_1D         # Net storage heat flux in all materials (increase)
        # Additional energy balance components
        o.variables['Qg'][:]           = data['Qg'].values   # Ground heat flux (downward)
        o.variables['Qanth_Qh'][:]     = no_data_1D   # Anthropogenic sensible heat flux (upward)
        o.variables['Qanth_Qle'][:]    = no_data_1D   # Anthropogenic latent heat flux (upward)
        o.variables['Qtau'][:]         = no_data_1D   # Momentum flux (downward)
        # General water balance components
        o.variables['Snowf'][:]        = data['Snowf'].values   # Snowfall rate (downward)
        o.variables['Rainf'][:]        = data['Rainf'].values   # Rainfall rate (downward)
        o.variables['Evap'][:]         = data['Evap'].values   # Total evapotranspiration (upward)
        o.variables['Qs'][:]           = data['Qs'].values   # Surface runoff (out of gridcell)
        o.variables['Qsb'][:]          = data['Qsb'].values   # Subsurface runoff (out of gridcell)
        o.variables['Qsm'][:]          = data['Qsm'].values   # Snowmelt (solid to liquid)
        o.variables['Qfz'][:]          = no_data_1D   # Re-freezing of water in the snow (liquid to solid)
        o.variables['DelSoilMoist'][:] = no_data_1D   # Change in soil moisture (increase)
        o.variables['DelSWE'][:]       = no_data_1D   # Change in snow water equivalent (increase)
        o.variables['DelIntercept'][:] = no_data_1D   # Change in interception storage (increase)
        o.variables['Qirrig'][:]       = no_data_1D   # Anthropogenic water flux from irrigation (increase)
        # Surface state variables
        o.variables['SnowT'][:]        = data['SnowT'].values   # Snow surface temperature
        o.variables['VegT'][:]         = data['VegT'].values   # Vegetation canopy temperature
        o.variables['BaresoilT'][:]    = data['BaresoilT'].values   # Temperature of bare soil (skin)
        o.variables['AvgSurfT'][:]     = no_data_1D   # Average surface temperature (skin)
        o.variables['RadT'][:]         = data['RadT'].values   # Surface radiative temperature
        o.variables['Albedo'][:]       = data['Albedo'].values   # Surface albedo
        o.variables['SWE'][:]          = data['SWE'].values   # Snow water equivalent
        o.variables['SurfStor'][:]     = no_data_1D   # Surface water storage
        o.variables['SnowFrac'][:]     = no_data_1D   # Snow covered fraction
        o.variables['SAlbedo'][:]      = no_data_1D   # Snow albedo
        o.variables['CAlbedo'][:]      = no_data_1D   # Vegetation canopy albedo
        o.variables['UAlbedo'][:]      = no_data_1D   # Urban canopy albedo
        o.variables['LAI'][:]          = no_data_1D   # Leaf area index
        o.variables['RoofSurfT'][:]    = no_data_1D   # Roof surface temperature (skin)
        o.variables['WallSurfT'][:]    = no_data_1D   # Wall surface temperature (skin)
        o.variables['RoadSurfT'][:]    = no_data_1D   # Road surface temperature (skin)
        o.variables['TairSurf'][:]     = no_data_1D   # Near surface air temperature (2m)
        o.variables['TairCanyon'][:]   = no_data_1D   # Air temperature in street canyon (bulk)
        o.variables['TairBuilding'][:] = no_data_1D   # Air temperature in buildings (bulk)
        # Sub-surface state variables **** TWO DIMENSIONAL ****
        o.variables['SoilMoist'][:,0]  = data['SoilMoist1'].values   # Average layer soil moisture
        o.variables['SoilMoist'][:,1]  = data['SoilMoist2'].values   # Average layer soil moisture
        o.variables['SoilMoist'][:,2]  = data['SoilMoist3'].values   # Average layer soil moisture
        o.variables['SoilMoist'][:,3]  = data['SoilMoist4'].values   # Average layer soil moisture
        o.variables['SoilMoist'][:,4]  = data['SoilMoist5'].values   # Average layer soil moisture
        o.variables['SoilMoist'][:,5]  = data['SoilMoist6'].values   # Average layer soil moisture
        o.variables['SoilTemp'][:,0]   = data['SoilTemp1'].values   # Average layer soil temperature
        o.variables['SoilTemp'][:,1]   = data['SoilTemp2'].values   # Average layer soil temperature
        o.variables['SoilTemp'][:,2]   = data['SoilTemp3'].values   # Average layer soil temperature
        o.variables['SoilTemp'][:,3]   = data['SoilTemp4'].values   # Average layer soil temperature
        o.variables['SoilTemp'][:,4]   = data['SoilTemp5'].values   # Average layer soil temperature
        o.variables['SoilTemp'][:,5]   = data['SoilTemp6'].values   # Average layer soil temperature

        # Evaporation components
        o.variables['TVeg'][:]         = data['TVeg'].values   # Vegetation transpiration
        o.variables['ESoil'][:]        = data['ESoil'].values   # Bare soil evaporation
        o.variables['RootMoist'][:]    = data['RootMoist'].values   # Root zone soil moisture
        o.variables['SoilWet'][:]      = no_data_1D   # Total soil wetness
        o.variables['ACond'][:]        = no_data_1D   # Aerodynamic conductance
        # Forcing data (at forcing height)
        o.variables['SWdown'][:]       = data['SWdown'].values   # Downward shortwave radiation
        o.variables['LWdown'][:]       = data['LWdown'].values   # Downward longwave radiation
        o.variables['Tair'][:]         = data['Tair'].values   # Air temperature
        o.variables['Qair'][:]         = data['Qair'].values   # Specific humidity
        o.variables['PSurf'][:]        = data['PSurf'].values   # Air pressure
        o.variables['Wind'][:]         = data['Wind'].values   # Wind speed

    return

###############################################################################
##### OTHER FUNCTIONS: These standard functions shouldn't need to be altered
###############################################################################

def set_more_info(info):
    '''sets timing and other info from forcing file and places into info dictionary.

    Inputs
    ------
    info (dictionary): script information

    Outputs
    -------
    info (dictionary): additional script information
    '''

    fpath = '%s/%s' %(info['sitepath'],info['fname_forcing'])

    with nc.Dataset(filename=fpath, mode='r', format='NETCDF4') as f:

        info['time_coverage_start']       = f.time_coverage_start
        info['time_coverage_end']         = f.time_coverage_end
        info['time_analysis_start']       = f.time_analysis_start
        info['local_utc_offset_hours']    = f.local_utc_offset_hours
        info['timestep_interval_seconds'] = f.timestep_interval_seconds
        info['timestep_number_spinup']    = f.timestep_number_spinup
        info['timestep_number_analysis']  = f.timestep_number_analysis

    # loading site parameters
    fpath = '%s/%s_sitedata_v1.csv' %(info['sitepath'], info['sitename'] )
    sitedata_full = pd.read_csv(fpath, index_col=1, delimiter=',')
    sitedata      = pd.to_numeric(sitedata_full['value'])

    info['latitude']      = sitedata['latitude']
    info['longitude']     = sitedata['longitude']

    return sitedata,info

###############################################################################

def create_empty_netcdf(info):
    '''creates empty netcdf dataset complying with Urban-PLUMBER protocol v1.0

    Inputs
    ------
    info (dictionary): script information

    '''

    timesteps = info['timestep_number_spinup'] + info['timestep_number_analysis']

    #fpath = '%s/%s' %(info['outpath'],"%s_post_processed.nc" % (info['sitename']))

    # open netcdf files (r = read only, w = write new)
    with nc.Dataset(filename=info['post_fname'], mode='w', format='NETCDF4') as o:

        # setting coordinate values
        times         = [t*int(info['timestep_interval_seconds']) for t in range(0,int(timesteps))]
        soil_layers   = [i for i in range(1,info['num_soil_layers']+1)]

        ############ create dimensions ############
        o.createDimension(dimname='time', size=timesteps)
        o.createDimension(dimname='soil_layer', size=info['num_soil_layers'])
        o.createDimension(dimname='x', size=1)
        o.createDimension(dimname='y', size=1)

        ############ create coordinates ############
        var = 'time'
        o.createVariable(var, datatype='i4', dimensions=('time'), fill_value = missing_float)
        o.variables[var].long_name     = 'Time'
        o.variables[var].standard_name = 'time'
        o.variables[var].units         = 'seconds since %s' %info['time_coverage_start']
        o.variables[var].calendar      = 'standard'
        o.variables[var][:]            = times

        var = 'soil_layer'
        o.createVariable(var, datatype='i4', dimensions=('soil_layer'), fill_value = missing_float)
        o.variables[var].long_name     = 'Soil layer number'
        o.variables[var][:]            = soil_layers

        var = 'x'
        o.createVariable(var, datatype='i4', dimensions=('x'), fill_value = missing_float)
        o.variables[var].long_name     = 'x dimension'
        o.variables[var][:]            = 1

        var = 'y'
        o.createVariable(var, datatype='i4', dimensions=('y'), fill_value = missing_float)
        o.variables[var].long_name     = 'y dimension'
        o.variables[var][:]            = 1

        ################### latidude and longitude ###################

        var = 'longitude'
        o.createVariable(var, datatype='f8', dimensions=('y', 'x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Longitude'
        o.variables[var].standard_name = 'longitude'
        o.variables[var].units         = 'degrees_east'
        o.variables[var][:]            = info['longitude']

        var = 'latitude'
        o.createVariable(var, datatype='f8', dimensions=('y', 'x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Latitude'
        o.variables[var].standard_name = 'latitude'
        o.variables[var].units         = 'degrees_north'
        o.variables[var][:]            = info['latitude']

        ##########################################################################
        ################### critical energy balance components ###################

        var = 'SWnet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net shortwave radiation (downward)'
        o.variables[var].standard_name = 'surface_net_downward_shortwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'LWnet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net longwave radiation (downward)'
        o.variables[var].standard_name = 'surface_net_downward_longwave_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qle'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Latent heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_latent_heat_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qh'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Sensible heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_sensible_heat_flux'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qstor'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Net storage heat flux in all materials (increase)'
        o.variables[var].standard_name = 'surface_thermal_storage_heat_flux'
        o.variables[var].units         = 'W/m2'

        ################### additional energy balance compoenents #################

        var = 'Qg'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Ground heat flux (downward)'
        o.variables[var].standard_name = 'downward_heat_flux_at_ground_level_in_soil'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth_Qh'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic sensible heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_sensible_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qanth_Qle'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic latent heat flux (upward)'
        o.variables[var].standard_name = 'surface_upward_latent_heat_flux_due_to_anthropogenic_energy_consumption'
        o.variables[var].units         = 'W/m2'

        var = 'Qtau'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Momentum flux (downward)'
        o.variables[var].standard_name = 'magnitude_of_surface_downward_stress'
        o.variables[var].units         = 'N/m2'

        ##################### general water balance components #####################

        var = 'Snowf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snowfall rate (downward)'
        o.variables[var].standard_name = 'snowfall_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Rainf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Rainfall rate (downward)'
        o.variables[var].standard_name = 'rainfall_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Evap'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Total evapotranspiration (upward)'
        o.variables[var].standard_name = 'surface_evapotranspiration'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qs'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface runoff (out of gridcell)'
        o.variables[var].standard_name = 'surface_runoff_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qsb'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Subsurface runoff (out of gridcell)'
        o.variables[var].standard_name = 'subsurface_runoff_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qsm'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snowmelt (solid to liquid)'
        o.variables[var].standard_name = 'surface_snow_and_ice_melt_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'Qfz'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Re-freezing of water in the snow (liquid to solid)'
        o.variables[var].standard_name = 'surface_snow_and_ice_refreezing_flux'
        o.variables[var].units         = 'kg/m2/s'

        var = 'DelSoilMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in soil moisture (increase)'
        o.variables[var].standard_name = 'change_over_time_in_mass_content_of_water_in_soil'
        o.variables[var].units         = 'kg/m2'

        var = 'DelSWE'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in snow water equivalent (increase)'
        o.variables[var].standard_name = 'change_over_time_in_surface_snow_and_ice_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'DelIntercept'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Change in interception storage (increase)'
        o.variables[var].standard_name = 'change_over_time_in_canopy_water_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'Qirrig'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Anthropogenic water flux from irrigation (increase)'
        o.variables[var].standard_name = 'surface_downward_mass_flux_of_water_due_to_irrigation'
        o.variables[var].units         = 'kg/m2/s'

        ########################## surface state variables ########################

        var = 'SnowT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow surface temperature'
        o.variables[var].standard_name = 'surface_snow_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'snow'

        var = 'VegT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation canopy temperature'
        o.variables[var].standard_name = 'surface_canopy_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'vegetation'

        var = 'BaresoilT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Temperature of bare soil'
        o.variables[var].standard_name = 'surface_ground_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'baresoil'

        var = 'AvgSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average surface temperature (skin)'
        o.variables[var].standard_name = 'surface_temperature'
        o.variables[var].units         = 'K'

        var = 'RadT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface radiative temperature'
        o.variables[var].standard_name = 'surface_radiative_temperature'
        o.variables[var].units         = 'K'

        var = 'Albedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface albedo'
        o.variables[var].standard_name = 'surface_albedo'
        o.variables[var].units         = '1'

        var = 'SWE'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow water equivalent'
        o.variables[var].standard_name = 'surface_snow_amount'
        o.variables[var].units         = 'kg/m2'

        var = 'SurfStor'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Surface water storage'
        o.variables[var].standard_name = 'surface_water_amount_assuming_no_snow'
        o.variables[var].units         = 'kg/m2'

        var = 'SnowFrac'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow covered fraction'
        o.variables[var].standard_name = 'surface_snow_area_fraction'
        o.variables[var].units         = '1'

        var = 'SAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Snow albedo'
        o.variables[var].standard_name = 'snow_and_ice_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'snow'

        var = 'CAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation canopy albedo'
        o.variables[var].standard_name = 'canopy_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'vegetation'

        var = 'UAlbedo'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Urban canopy albedo'
        o.variables[var].standard_name = 'urban_albedo'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'urban'

        var = 'LAI'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Leaf area index'
        o.variables[var].standard_name = 'leaf_area_index'
        o.variables[var].units         = '1'
        o.variables[var].subgrid       = 'vegetation'

        var = 'RoofSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Roof surface temperature (skin)'
        o.variables[var].standard_name = 'surface_roof_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'roof'

        var = 'WallSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Wall surface temperature (skin)'
        o.variables[var].standard_name = 'surface_wall_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'wall'

        var = 'RoadSurfT'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Road surface temperature (skin)'
        o.variables[var].standard_name = 'surface_road_skin_temperature'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'road'

        var = 'TairSurf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Near surface air temperature (2m)'
        o.variables[var].standard_name = 'air_temperature_near_surface'
        o.variables[var].units         = 'K'

        var = 'TairCanyon'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature in street canyon (bulk)'
        o.variables[var].standard_name = 'air_temperature_in_street_canyon'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'canyon'

        var = 'TairBuilding'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature in buildings (bulk)'
        o.variables[var].standard_name = 'air_temperature_in_buildings'
        o.variables[var].units         = 'K'
        o.variables[var].subgrid       = 'building'

        ######################## Sub-surface state variables ######################

        var = 'SoilMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','soil_layer','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average layer soil moisture'
        o.variables[var].standard_name = 'moisture_content_of_soil_layer'
        o.variables[var].units         = 'kg/m2'

        var = 'SoilTemp'
        o.createVariable(var, datatype='f8', dimensions=('time','soil_layer','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Average layer soil temperature'
        o.variables[var].standard_name = 'soil_temperature'
        o.variables[var].units         = 'K'

        ########################## Evaporation components #########################

        var = 'TVeg'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Vegetation transpiration'
        o.variables[var].standard_name = 'transpiration_flux'
        o.variables[var].units         = 'kg/m2/s'
        o.variables[var].subgrid       = 'vegetation'

        var = 'ESoil'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Bare soil evaporation'
        o.variables[var].standard_name = 'liquid_water_evaporation_flux_from_soil'
        o.variables[var].units         = 'kg/m2/s'
        o.variables[var].subgrid       = 'baresoil'

        var = 'RootMoist'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Root zone soil moisture'
        o.variables[var].standard_name = 'mass_content_of_water_in_soil_defined_by_root_depth'
        o.variables[var].units         = 'kg/m2'

        var = 'SoilWet'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Total soil wetness'
        o.variables[var].standard_name = 'relative_soil_moisture_content_above_wilting_point'
        o.variables[var].units         = '1'

        var = 'ACond'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Aerodynamic conductance'
        o.variables[var].standard_name = 'inverse_aerodynamic_resistance'
        o.variables[var].units         = 'm/s'

        ########################## forcing data variables #########################

        var = 'SWdown'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Downward shortwave radiation at measurement height'
        o.variables[var].standard_name = 'surface_downwelling_shortwave_flux_in_air'
        o.variables[var].units         = 'W/m2'

        var = 'LWdown'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Downward longwave radiation at measurement height'
        o.variables[var].standard_name = 'surface_downwelling_longwave_flux_in_air'
        o.variables[var].units         = 'W/m2'

        var = 'Tair'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air temperature at measurement height'
        o.variables[var].standard_name = 'air_temperature'
        o.variables[var].units         = 'K'

        var = 'Qair'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Specific humidity at measurement height'
        o.variables[var].standard_name = 'surface_specific_humidity'
        o.variables[var].units         = '1'

        var = 'PSurf'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Air pressure at measurement height'
        o.variables[var].standard_name = 'surface_air_pressure'
        o.variables[var].units         = 'Pa'

        var = 'Wind'
        o.createVariable(var, datatype='f8', dimensions=('time','y','x'), fill_value = missing_float)
        o.variables[var].long_name     = 'Wind speed at measurement height'
        o.variables[var].standard_name = 'wind_speed'
        o.variables[var].units         = 'm/s'

    return

#############################################################################

if __name__ == "__main__":

    for sitename in sitelist:

        main(sitename)
