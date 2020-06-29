#!/usr/bin/env python

"""
Plot visual benchmark (average seasonal cycle) of old vs new model runs.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.10.2017)"
__email__ = "mdekauwe@gmail.com"

import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import datetime
import os
import glob
from optparse import OptionParser
import string

def main(site, fname, plot_fname=None):

    df1 = read_cable_file(fname, type="CABLE")
    df1 = resample_timestep(df1, type="CABLE")

    pd.plotting.register_matplotlib_converters()

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    labels_gen = label_generator('lower', start="(", end=")")
    colours = plt.cm.Set2(np.linspace(0, 1, 7))

    ax1 = fig.add_subplot(1,1,1)

    axes = [ax1]
    vars = ["Qle"]

    props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                 ec="white")
    for v in vars:
        print(v)
        ax1.plot(df1[v].index.to_pydatetime(), df1[v].rolling(window=5).mean(),
               lw=1.5, ls="-", label=v)


    ax1.set_ylabel("LE (W m$^{-2}$)", fontsize=12)

    from matplotlib.ticker import MaxNLocator
    ax1.yaxis.set_major_locator(MaxNLocator(5))
    ax1.legend(numpoints=1, loc="best")

    if plot_fname is None:
        plt.show()
    else:
        #fig.autofmt_xdate()
        fig.savefig(plot_fname, bbox_inches='tight', pad_inches=0.1)


def read_cable_file(fname, type=None):

    f = nc.Dataset(fname)
    time = nc.num2date(f.variables['time'][:],
                        f.variables['time'].units)

    df = pd.DataFrame(f.variables['Qle'][:,0,0], columns=['Qle'])

    df['dates'] = time
    df = df.set_index('dates')

    return df



def resample_timestep(df, type=None):

    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0
    SEC_2_HLFHOUR = 1800.
    SEC_2_HOUR = 3600.

    if type == "CABLE":

        method = {'Qle':'mean'}
    elif type == "FLUX":
        # umol/m2/s -> g/C/60min
        #df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHOUR
        df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HOUR

        method = {'GPP':'sum', "Qle":"mean"}

    elif type == "MET":
        # kg/m2/s -> mm/60min
        df['Rainf'] *= SEC_2_HOUR

        method = {'Rainf':'sum'}

    df = df.resample("D").agg(method)

    return df

def label_generator(case='lower', start='', end=''):
    choose_type = {'lower': string.ascii_lowercase,
                   'upper': string.ascii_uppercase}
    generator = ('%s%s%s' %(start, letter, end) for letter in choose_type[case])

    return generator

if __name__ == "__main__":

    site = "AU-Preston"
    output_dir = "outputs"
    fname = "outputs/AU-Preston_out.nc"
    odir = "plots"
    main(site, fname)
