#!/usr/bin/env python3

"""
Script to represent time graphs 

time_graphs.py

Script purpose: the script represents time graphs for different rows
(variables) in an Excel matrix
"""

__version__ = '1.0'
_verdata = 'Apr 2020'

## Import libraries
_DEFAULT_ROUTE = 'data'
_DEFAULT_RESULTS = 'results'
_DEFAULT_FREQLIM_SUFFIXDIR = ''
_version_pandas = '0.16.0'
_version_matplotlib = '1.5.0'
_version_numpy = '1.16'

try:
    import argparse
    import re
    import sys
    import os
    import distutils.version
    from os import mkdir, listdir, rename
except ImportError:
    raise ImportError("This module requires a Linux distribution as OS")

def compare_versions(a, b):
    if a:
        a = distutils.version.LooseVersion(a)
        b = distutils.version.LooseVersion(b)
        return a >= b
    else:
        return False

def ChkDir(path):
    error = False
    try:
        os.makedirs(path)
    except(OSError):
        if not os.access(path, (os.W_OK | os.X_OK)):
            print('\nWARNING! Unable to write results in directory "' +
                path + '"...')
            error = True
    return(error)

try:
    import pandas as pd
    from pandas.plotting import register_matplotlib_converters
    register_matplotlib_converters()
except ImportError:
    raise ImportError("This module requires pandas")
else:
    if not compare_versions(pd.__version__, _version_pandas):
        raise ImportError('pandas %s or later is required; you have %s'
            % (_version_pandas, matplotlib.__version__))

try:
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    from mpl_toolkits import mplot3d
    import matplotlib.pylab as pylab
except ImportError:
    raise ImportError("This module requires matplotlib")
else:
    if not compare_versions(mpl.__version__, _version_matplotlib):
        raise ImportError(
            'matplotlib %s or later is required; you have %s' % (
                _version_matplotlib, mpl.__version__))

try:
    import numpy as np
except ImportError:
    raise ImportError("This module requires numpy")
else:
    if not compare_versions(np.__version__, _version_numpy):
        raise ImportError(
            'numpy %s or later is required; you have %s' % (
                _version_numpy, np.__version__))

## Create functions
# Read files from a given folder (return file names)
def read_file(route="data"):
    files = []
    for filename in listdir(route):
        if filename.split(".")[-1] in ['xlsx', 'xlsm', 'xltx', 'xltm']:
            files.append(filename)
    if len(files) == 0:
        print("\nFolder '" + route + "' is empty\n")
        sys.exit()
    return(files)

# Plot time graphs (dates as columns, species as rows)
def timegraph(filename, results="results", dates=True, conversion=False):
    print("\nTime graphs for", filename.split("/")[-1], "\n")
    sheets = pd.read_excel(filename, sheet_name=None, index_col=0)
    sheet_name, df = next(iter(sheets.items()))
    colnames = df.columns.values
    rownames = df.index.values
    if conversion:
        for row in rownames:
            prev = 0
            for i, col in enumerate(colnames):
                if df.loc[row, col] == 0:
                    if prev == 0:
                        df.loc[row, col] = np.NaN
                    else:
                        cont = 1
                        nex = 0
                        while nex == 0:
                            if (i + cont) == len(colnames):
                                break
                            nex = df.loc[row, colnames[i + cont]]
                            cont += 1
                        if nex != 0:
                            df.loc[row, col] = prev + (nex - prev) / cont
                            prev = df.loc[row, col]
                        else:
                            df.loc[row, col] = np.NaN
                else:
                    prev = df.loc[row, col]
    params = {'legend.fontsize': 'x-large',
            'figure.figsize': (15, 10),
            'axes.labelsize': 'x-large',
            'axes.titlesize': 30,
            'xtick.labelsize': 'x-large',
            'ytick.labelsize': 'x-large'}
    pylab.rcParams.update(params)
    for specie in rownames:
        new_df = df.loc[specie, :]
        fig, ax = plt.subplots()
        ax.plot(new_df, color="black")
        ax.set(title="Time graph " + specie,
            ylabel="Value", xlabel="Time")
        if dates:
            from matplotlib.dates import DateFormatter
            date_form = DateFormatter("%m-%y") 
            ax.xaxis.set_major_formatter(date_form)
        fig.tight_layout()
        plt.savefig(results + "/" +
            specie.replace("/", "|").replace(" ", "_") +
            "_timegraph.png")
        plt.close()
    fig, ax = plt.subplots()
    ax.set(title="Time graph " +
            filename.split("/")[-1].split(".")[0],
            ylabel="Value", xlabel="Time")
    if dates:
        from matplotlib.dates import DateFormatter
        date_form = DateFormatter("%m-%y") 
        ax.xaxis.set_major_formatter(date_form)
    for specie in rownames:
        new_df = df.loc[specie, :]
        ax.plot(new_df)
        lg = ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig(results + "/" + filename.split("/")[-1].split(".")[0] +
        "_timegraph.png", bbox_extra_artists=(lg,), bbox_inches='tight')
    plt.close()

## Main script
def main():
    ## Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Script to plot time graphs (samples as columns, sp\
ecies or variables as rows)',
        epilog=('Script to plot time graphs - by Torres A - ' +_verdata)
    )
    parser.add_argument(
        '-d', '--data',
        action='store',
        default=_DEFAULT_ROUTE,
        help=('directory path with data files (if not present, \'' +
        _DEFAULT_ROUTE + '\' will be tried), in excel format')
    )
    parser.add_argument(
        '-r', '--results',
        action='store',
        default=_DEFAULT_RESULTS,
        help=('path of the directory for saving results, relative to th\
e data directory (if not present, \'' + _DEFAULT_RESULTS + '\' will be \
tried)')
    )
    parser.add_argument(
        '-m', '--method',
        action='store',
        default="years",
        choices=["dates", "years"],
        help=('Use dates or years as axis (default: years)')
    )
    parser.add_argument(
        '-z', '--zeros',
        action='store',
        default="no",
        choices=["yes", "no"],
        help=('Convert zeros in data in values. It converts zeros in Na\
N values if they are initial or final zeros, but it converts zeros betw\
een two values in medium values (default: no)')
    )

    ## Parse arguments
    flgs = parser.parse_args()
    datapath = flgs.data
    method = flgs.method
    zeros = flgs.zeros
    _DEFAULT_FREQLIM_SUFFIXDIR = ''
    if (flgs.data == _DEFAULT_ROUTE):
        route = flgs.data + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        route = _DEFAULT_ROUTE
    if (flgs.results == _DEFAULT_RESULTS):
        results = flgs.results + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        results = _DEFAULT_RESULTS
    if (ChkDir(results)):
        results = ''
        print('Warning: the data root directory will be tried instead, \
for creating results folder')
    if (ChkDir(route)):
        route = ''
        print('Warning: the data root directory will be tried instead, \
for creating data folder')

    ## Execute functions
    if method == "dates":
        dates = True
    else:
        dates = False
    if zeros == "yes":
        conversion = True
    else:
        conversion = False
    files = read_file(route)
    for filename in files:
        filename = route + "/" + filename
        timegraph(filename, results, dates, conversion)

if __name__ == "__main__":
    main()
