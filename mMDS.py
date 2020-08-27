#!/usr/bin/env python3

"""
Script to represent metric multidimensional scaling (mMDS) analisys

mMDS.py

Script purpose: the script represents mMDS analisys graphically
"""

__version__ = '1.0'
_verdata = 'Jan 2020'

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
    from math import sqrt, log10
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

try:
    from sklearn.datasets.samples_generator import make_s_curve
    from sklearn.datasets import load_digits
    from sklearn.manifold import MDS
    from sklearn import manifold
    from sklearn.metrics import euclidean_distances
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
except ImportError:
    raise ImportError("This module requires sklearn")

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

# Standarize 
def standVals(df):
    maxelm = len(df.index.values)
    means = []
    stds = []
    for i in range(len(df.columns.values)):
        mean = df.loc[:,df.columns.values[i]].mean()
        std = df.loc[:,df.columns.values[i]].std()
        means.append(mean)
        stds.append(std)
    for i in range(len(df.index.values)):
        for j in range(len(df.columns.values)):
            df.loc[df.index.values[i], df.columns.values[j]] = (
                df.loc[df.index.values[i], df.columns.values[j]] -
                means[j]) / stds[j]
    return(df)

# Calculate unitary values
def unitVals(df):
    maxelm = len(df.index.values)
    for i in range(len(df.index.values)):
        for j in range(len(df.columns.values)):
            df.loc[df.index.values[i], df.columns.values[j]] = (
                df.loc[df.index.values[i], df.columns.values[j]] /
                sqrt(maxelm - 1))
    return(df)

# Calculate scalar product
def scalar_product(df, var1, var2):
    scalar = 0
    for i in range(len(df.index.values)):
        scalar1 = df.loc[df.index.values[i], df.columns.values[var1]]
        scalar2 = df.loc[df.index.values[i], df.columns.values[var2]]
        scalar += (scalar1 * scalar2)
    return(scalar)

# Calculate correlations matrix
def corr_matrix(df):
    maxelm = len(df.columns.values)
    df['SUM'] = df.sum(axis='columns')
    df = df.sort_values(by='SUM', ascending=False)
    df.drop('SUM', axis='columns', inplace=True)
    corrdf = np.zeros((maxelm, maxelm))
    corrdf = pd.DataFrame(corrdf, index=df.columns.values,
        columns=df.columns.values)
    for i in range(len(corrdf.index.values)):
        for j in range(len(corrdf.columns.values)):
            corrdf.loc[corrdf.index.values[i], corrdf.columns.values[j]
                ] = scalar_product(df, i, j)
    df = corrdf.iloc[:maxelm, :maxelm]
    return(df)

# Calculate distances between correlations
def dist(x):
    x = (round(x * 100000)) / 100000
    return(sqrt(1 - x) / sqrt(2))

# Calculate distances matrix
def dist_matrix(df):
    for i in range(len(df.index.values)):
        for j in range(len(df.columns.values)):
            df.loc[df.index.values[i], df.columns.values[j]] = dist(
                df.loc[df.index.values[i], df.columns.values[j]])
    return(df)

# Compute classical MDS (Francis Song)
def cmdscale(df):
    n = len(df)
    H = np.eye(n) - np.ones((n, n))/n
    B = -H.dot(df**2).dot(H)/2
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
    w, = np.where(evals > 0)
    L = np.diag(np.sqrt(evals[w]))
    V = evecs[:,w]
    Y = V.dot(L)
    return(Y, evals)

# Plot classical MDS (it represents 2D and 3D mMDS from an Excel matrix)
def plotMDS_M(filename, results="results", dim1=1, dim2=2, dim3=3,
    transpose=False, smacof=False, groups=False, sep=10, size=True,
    angle1=60, angle2=60):
    print("\nMetric MDS for", filename, "\n")
    sheets = pd.read_excel(filename, sheet_name=None, index_col=0)
    sheet_name, df = next(iter(sheets.items()))
    if groups:
        groups_col = df[df.columns.values[-1]]
        df = df.drop([df.columns.values[-1]], axis=1)
    if transpose:
        df = df.T
    print(df)
    if size == True:
        sizes = []
        for value in df.columns.values:
            sizes.append(df.loc[:,value].mean())
        minimum = min(sizes)
        if minimum < 0:
            minimum = minimum * (-1)
        length = 0
        while minimum < 1:
            length += 1
            minimum = minimum * 10
        for i in range(len(sizes)):
            sizes[i] = log10(sizes[i] * 10 ** length) * 10
    else:
        sizes = [30] * len(df.columns.values)
        length = 2
    df = standVals(df)
    df = unitVals(df)
    df = corr_matrix(df)
    df = dist_matrix(df)
    df.to_excel(results + "/" + filename.split("/")[-1].split(".")[0] +
        "_dists." + filename.split(".")[-1], sheet_name, header=True,
        index=False)
    if not smacof:
        comp, var = cmdscale(df)
        comps = []
        for i in range(len(comp)):
            comps.append(list(comp[i, [dim1, dim2, dim3]]))
        print("\nTotal variances:\n", var)
        sumv = sum(var)
        print("\nVariance explained by component", dim1, "-->",
            var[dim1-1] / sumv * 100, "%")
        print("Variance explained by component", dim2, "-->",
            var[dim2-1] / sumv * 100, "%")
        print("Variance explained by component", dim3, "-->",
            var[dim3-1] / sumv * 100, "%\n")
        principalDf = pd.DataFrame(data = comps, columns = [
            'component ' + str(dim1), 'component '
            + str(dim2), 'component ' + str(dim3)])
    else:
        EPSILON = np.finfo(np.float32).eps
        n_samples = len(df.index.values)
        embedding = MDS(n_components=3, metric=True, eps=EPSILON,
            dissimilarity="precomputed", max_iter=3000)
        df_transformed = embedding.fit_transform(df[:n_samples])
        principalDf = pd.DataFrame(data = df_transformed,
            columns = ['component ' + str(dim1), 'component '
            + str(dim2), 'component ' + str(dim3)])
    targets = features = df.index.values
    finalDf = principalDf.assign(target = features)
    X = np.array(principalDf)
    size = 30
    params = {'legend.fontsize': size,
        'legend.loc': 'center left',
        'figure.figsize': (20, 20),
        'axes.labelsize': size,
        'axes.titlesize': size,
        'xtick.labelsize': size*0.75,
        'ytick.labelsize': size*0.75,
        'axes.labelpad': size,
        'axes.titlepad': 25}
    plt.rcParams.update(params)
    colors = ["black", "green", "blue", "red", "yellow", "magenta",
        "cyan", "silver", "brown", "orange", "gold", "yellowgreen",
        "lime", "turquoise", "teal", "skyblue", "slategrey", "navy",
        "blueviolet", "indigo", "darkviolet", "purple", "pink",
        "lightgrey", "rosybrown", "mistyrose", "coral", "khaki",
        "darkseagreen", "springgreen", "lightblue", "thistle"]
    if groups:
        color = []
        grouped = []
        pos_grouped = []
        numcol = 0
        for group in range(len(groups_col)):
            if groups_col[group] not in grouped:
                grouped.append(groups_col[group])
                pos_grouped.append(group)
                color.append(colors[numcol])
                numcol += 1
            else:
                index = grouped.index(groups_col[group])
                color.append(color[pos_grouped[index]])
    else:
        color = []
        for i in range(len(targets)):
            pos = i % len(colors)
            color.append(colors[pos])
    color = np.array(color)
    sep = (sep / 10) * 100 ** (-(length))
    ## 2D Matplotlib plot
    fig, ax = plt.subplots()
    ax.set(title="2D mMDS " + filename.split("/")[-1].split(".")[0],
        ylabel='Dimension ' + str(dim2),
        xlabel='Dimension ' + str(dim1))
    ax.scatter(X[:, 0], X[:, 1], c=color, s=30)
    points = []
    names = []
    for i in range(len(targets)):
        ax.text(X[i][0] + sep, X[i][1] + sep, i+1, size=size*0.75)
        points.append(mpl.lines.Line2D(X[i], X[0], linestyle="none",
            c=color[i], marker = 'o'))
        names.append(str(i+1) + " - " + str(targets[i]))
    lg = ax.legend(points, names, bbox_to_anchor=(0, 1), loc='upper right')
    plt.savefig(results + "/" + filename.split("/")[-1].split(".")[-2] +
        "_mMDS_2D.jpg", bbox_extra_artists=(lg,), bbox_inches='tight')
    plt.close()
    ## 3D Matplotlib plot
    sizes2 = [0] * len(sizes)
    for i in range(len(sizes)):
        sizes2[i] = sizes[i] * 10
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.scatter3D(X[:, 0], X[:, 1], X[:, 2], c=color, cmap=plt.cm.jet,
        s=sizes2)
    ax.set_xlabel('Dimension ' + str(dim1))
    ax.set_ylabel('Dimension ' + str(dim2))
    ax.set_zlabel('Dimension ' + str(dim3))
    ax.set_title('3D mMDS')
    points = []
    names = []
    for i in range(len(targets)):
        ax.text(X[i][0] + sep, X[i][1] + sep, X[i][2] + sep, i+1,
            size=size*0.75)
        points.append(mpl.lines.Line2D(X[i], X[0], linestyle="none",
            c=color[i], marker = 'o'))
        names.append(str(i+1) + " - " + str(targets[i]))
    lg = ax.legend(points, names, numpoints=1, bbox_to_anchor=(0, 1),
        loc="upper right")
    ax.view_init(angle1, angle2)
    plt.savefig(results + "/" + filename.split("/")[-1].split(".")[-2] +
        "_mMDS_3D.jpg", bbox_extra_artists=(lg,), bbox_inches='tight')
    plt.close()

## Main script
def main():
    ## Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Script to plot classical metric multidimensional s\
caling (mMDS) from correlation matrix',
        epilog=('Script to plot mMDS from correlation matrix - by Torre\
s A - ' + _verdata)
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
        default="mdscale",
        choices=["smacof", "mdscale"],
        help=('Use SMACOF or mdscale method for mds (default: mds)')
    )
    parser.add_argument(
        '-g', '--groups',
        action='store',
        default="no",
        choices=["yes", "no"],
        help=('"yes" for use last column of excel files for grouping va\
riables, "no" for not grouping (default: "no"). Maximum number of group\
s that can be separated is 32.')
    )
    parser.add_argument(
        '-s', '--size',
        action='store',
        default="no",
        choices=["yes", "no"],
        help=('"yes" for representing points sizes depending on variabl\
es values (default: "no").')
    )
    parser.add_argument(
        '-t', '--transpose',
        action='store',
        default="yes",
        choices=["yes", "no"],
        help=('"yes" for transpose excel matrix, "no" for not transposi\
ng (default: "yes", transpose excel with samples as rows and variables \
as columns). If present, it also transposes groups column')
    )
    parser.add_argument(
        '-p', '--separation',
        action='store',
        default=10,
        help=('Separation between points and point names in graphs (def\
ault: 10)')
    )
    parser.add_argument(
        '-1', '--angle1',
        action='store',
        default=None,
        help=('Angle 1 for multivariate analysis representation (defaul\
t: None)')
    )
    parser.add_argument(
        '-2', '--angle2',
        action='store',
        default=None,
        help=('Angle 2 for multivariate analysis representation (defaul\
t: None)')
    )

    ## Parse arguments
    flgs = parser.parse_args()
    datapath = flgs.data
    method = flgs.method
    groups = flgs.groups
    transpose = flgs.transpose
    size = flgs.size
    sep = flgs.separation
    a1 = flgs.angle1
    a2 = flgs.angle2
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
    if method == "smacof":
        smacof = True
    else:
        smacof = False
    if groups != "no":
        groups = True
    else:
        groups = False
    if size != "no":
        size = True
    else:
        size = False
    if transpose != "no":
        transpose = True
    else:
        transpose = False
    try:
        sep = int(sep)
        if a1 != None:
            a1 = int(a1)
        if a2 != None:
            a2 = int(a2)
    except ValueError:
        print("\nSeparation and angles arguments must be integers.\n")
        sys.exit()
    files = read_file(route)
    for filename in files:
        filename = route + "/" + filename
        plotMDS_M(filename, results, dim1=1, dim2=2, dim3=3,
            transpose=transpose, smacof=smacof, groups=groups,
            size=size, sep=sep, angle1=a1, angle2=a2)

if __name__ == "__main__":
    main()
