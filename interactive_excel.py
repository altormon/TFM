#!/usr/bin/env python3

"""
Script to pretreat a matrix in an Excel file, for its use in
cmplxcruncher (cc)

interactive_excel.py

Script purpose: The script transposes the data introduced in an Excel
file, separates Excel sheets in an Excel file, converts NA in 0s,
converts absolute frequencies in relative frequencies, merge two Excel
dataframes with same dates (column names), separates Excels by dates,
separates Excel by species
"""

__version__ = '1.0'
_verdata = 'Sep 2019'

## Import libraries

_DEFAULT_ROUTE = './data'
_DEFAULT_OLD_ROUTE = './not_treated_data'
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
    from math import ceil, isnan
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
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator
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
# Move file (it moves no trasposed file(s) to another folder)
def move_file(filename, new_folder="./not_treated_data",
    new_filename=""):
    if "not_treated_data" not in listdir("."):
        mkdir(new_folder)
    if new_filename == "":
        new_filename = filename
    new_filename = new_folder + "/" + new_filename.split("/")[-1]
    rename(filename, new_filename)

# Read file (it returns a list with Excel workbooks or text files in
# folder 'data')
def read_file(route="./data", new_folder="./not_treated_data",
    new_filename=""):
    if "data" not in listdir("./"):
        print("You need to put your Excel or txt/tsv file(s) to separat\
e in a folder called 'data'")
        sys.exit()
    files = []
    for filename in listdir(route):
        if filename.split(".")[-1] in ['xlsx', 'xlsm', 'xltx', 'xltm']:
            files.append(filename)
            move_file(route + "/" + filename, new_folder = new_folder,
                new_filename=new_filename)
    return(files)

# Transpose File (it transposes Excel matrix in Excel workbooks in a
# list (input) to new Excel files in 'data' folder))
def transpose_file(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        df_trasposed = df.T
        df_trasposed.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_transposed." +
            filename.split(".")[-1], sheet_name, header=False,
            index=True)

# Separate Excel (it separates Excel sheets in Excel workbooks in a list
# (input) to new Excel files in 'data' folder)
def separate_excel(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        xls = pd.ExcelFile(old_route + "/" + filename)
        sheets = xls.sheet_names
        for sheet in sheets:
            df = pd.read_excel(old_route + "/" + filename,
                sheet_name=sheet)
            df.to_excel(new_route + "/" +
                filename.split("/")[-1].split(".")[0] + "_" + sheet +
                "." + filename.split(".")[-1], sheet, header=True,
                index=False)

# Fill file (it replaces Null (NA) values in a dataframe)
def fill_file(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        df = df.fillna("0")
        df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_filled." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Relative frequencies (it obtains relative frequencies of species for
# all dates in dataframe)
def relative_freqs(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        rownames = df.index.values
        colnames = df.columns.values
        nrow = df.shape[0]
        ncol = df.shape[1]
        sums = []
        for col in range(ncol-1):
            col += 1
            sums.append(df[colnames[col]].sum(axis = 0))
        for row in range(nrow):
            for col in range(ncol-1):
                if sums[col] != 0:
                    df.loc[rownames[row], colnames[col+1]] = (
                        float(df.loc[rownames[row],
                        colnames[col+1]]) / sums[col])
                else:
                    df.loc[rownames[row], colnames[col+1]] = 0
        df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_relative." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Merge 2 files (it merges two dataframes by column name from two files
# in Excel format)
def merge_2files(file1, file2, new_route="./data",
    old_route="./not_treated_data"):
    sheets1 = pd.read_excel(old_route + "/" + file1, sheet_name=None)
    sheets2 = pd.read_excel(old_route + "/" + file2, sheet_name=None)
    sheet_name1, df1 = next(iter(sheets1.items()))
    sheet_name2, df2 = next(iter(sheets2.items()))
    frames = [df1, df2]
    colnames1 = df1.columns.values
    colnames2 = df2.columns.values
    names = (list(df1.loc[:, colnames1[0]]) +
        list(df2.loc[:, colnames2[0]]))
    df_concat_col = pd.concat(frames, sort=False, join='inner')
    df_concat_col.insert(loc=0, column="", value=names)
    if len(file1.split("/")[-1].split(".")[0] + "_" +
        file2.split("/")[-1].split(".")[0] + "." +
        file1.split(".")[-1]) <= 31:
        if (len(sheet_name1) + len(sheet_name2) + 1) <= 31:
            df_concat_col.to_excel(new_route + "/" +
                file1.split("/")[-1].split(".")[0] + "_" +
                file2.split("/")[-1].split(".")[0] + "." +
                file1.split(".")[-1], sheet_name=(sheet_name1 + "_" +
                sheet_name2), header=True, index=False)
        else:
            df_concat_col.to_excel(new_route + "/" +
                file1.split("/")[-1].split(".")[0] + "_" +
                file2.split("/")[-1].split(".")[0] + "." +
                file1.split(".")[-1], sheet_name="sheets_merged",
                header=True, index=False)
    else:
        if (len(sheet_name1) + len(sheet_name2) + 1) <= 31:
            df_concat_col.to_excel(new_route + "/" + "files_merged." +
                file1.split(".")[-1], sheet_name=(sheet_name1 + "_" +
                sheet_name2), header=True, index=False)
        else:
            df_concat_col.to_excel(new_route + "/" + "files_merged." +
                file1.split(".")[-1], sheet_name="sheets_merged",
                header=True, index=False)

# Group by dates (it groups data by indicated dates: months, years...)
def select_by_dates(files, new_route="./data",
    old_route="./not_treated_data", days=[], months=[], years=[]):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        colnames = df.columns.values[1:]
        if len(days) == 0:
            days = list(range(32))
        else:
            for i in range(len(days)):
                try:
                    days[i] = int(days[i])
                except ValueError:
                    raise ValueError("Days must be numbers between 0 an\
d 31 inclusive.")
                if days[i] not in list(range(32)):
                    print("Days must be numbers between 0 and 31 inclus\
ive.")
                    sys.exit()
        if len(months) == 0:
            months = list(range(13))
        else:
            for i in range(len(months)):
                try:
                    months[i] = int(months[i])
                except ValueError:
                    raise ValueError("Months must be numbers between 0 \
and 12 inclusive.")
                if months[i] not in list(range(13)):
                    print("Months must be numbers between 0 and 12 incl\
usive.")
                    sys.exit()
        if len(years) == 0:
            years = list(range(10000))
        else:
            for i in range(len(years)):
                try:
                    years[i] = int(years[i])
                except ValueError:
                    raise ValueError("Years must be numbers between 0 a\
nd 31 inclusive.")
                if years[i] not in list(range(10000)):
                    print("Years must be numbers between 0 and 9999 inc\
lusive.")
                    sys.exit()
        dates = []
        for date in range(len(colnames)):
            mask = ((colnames[date].day in days) and
                (colnames[date].month in months) and
                (colnames[date].year in years))
            if mask == True:
                dates.append(colnames[date])
        new_df = df.loc[:, dates]
        names = list(df.loc[:, df.columns.values[0]])
        new_df.insert(loc=0, column="", value=names)
        new_df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_dates." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Group by species (it groups data by indicated rows of species)
def select_by_species(files, species=[], new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        if len(species) == 0:
            species = []
            for i in range(len(df.index)):
                species.append(i+1)
        else:
            for i in range(len(species)):
                try:
                    species[i] = int(species[i])-1
                except ValueError:
                    raise ValueError("Rows of species must be numbers b\
etween 0 and highest row number inclusive.")
                if species[i] not in list(range(len(df.index)+1)):
                    print("Rows of species must be numbers between 0 an\
d highest row number inclusive.")
                    sys.exit()
        new_df = df.loc[species, :]
        new_df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_species." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Group by ZRF (it groups data by indicated percent of zeros)
def select_by_ZRF(files, core=100, tail=0, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        try:
            core = int(core)
            tail = int(tail)
        except ValueError:
            raise ValueError("Percents of core / tail must be number be\
tween 0 and 100 inclusive")
        if 0 >= core >= 100 or 0 >= tail >= 100:
            print("Percents of core / tail must be number between 0 and\
 100 inclusive")
            sys.exit()
        zeros = (df == 0).sum(axis=1)
        rows = []
        for i in df.index:
            zero = (zeros[i]/len(df.columns))*100
            if zero < core and zero >= tail:
                rows.append(i)
        new_df = df.loc[rows, :]
        new_df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "ZRF." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Calculate ZRF (it calculates ZRF values)
def calculate_ZRF(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        zeros = (df == 0).sum(axis=1)
        values = []
        for i in df.index:
            zero = (zeros[i]/len(df.columns))*100
            values.append(zero)
        new_df = df.iloc[:, [0]]
        new_df["ZRF"] = values
        new_df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_ZRF_values." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

# Round up (round a given number up, with decimals given)
def roundup(x, dec=0):
    y = 1
    if dec > 0:
        for i in range(dec):
            y = y * 10
        x = ceil(x * y) / y
    else:
        x = ceil(x)
    return(x)

# Round down (round a given number down, with decimals given)
def rounddown(x, dec=0):
    y = 1
    if dec > 0:
        for i in range(dec):
            y = y * 10
        x = int(x * y) / y
    else:
        x = int(x)
    return(x)

# Matplot values (it plots values X and Y given using matplotlib)
def matplot_values(x, ylist, yerrlist,  xname="Variable1",
    ynames=["Variable2"], new_route="./data", fsize=15, rot=30):
    xnames = x
    x = np.array(range(len(x)))
    if len(ylist) > 1 and len(yerrlist) > 1:
        for i in range(len(ylist)):
            ylist[i] = np.array(ylist[i])
        for i in range(len(yerrlist)):
            yerrlist[i] = np.array(yerrlist[i])
    else:
        ylist = np.array(ylist[0])
        yerrlist = np.array(yerrlist[0])
    plt.figure(figsize=(15,10))
    plt.subplot(111)
    if len(ynames) > 1:
        for i in range(len(ylist)):
            if i == 0:
                plt.errorbar(x, ylist[i], yerr=yerrlist[i], fmt="ob",
                    label=ynames[i], capsize=5, elinewidth=2,
                    capthick=2)
            elif i == 1:
                plt.errorbar(x, ylist[i], yerr=yerrlist[i], fmt="vr",
                    label=ynames[i], capsize=5, elinewidth=2,
                    capthick=2)
            elif i == 2:
                plt.errorbar(x, ylist[i], yerr=yerrlist[i], fmt="sg",
                    label=ynames[i], capsize=5, elinewidth=2,
                    capthick=2)
            else:
                plt.errorbar(x, ylist[i], yerr=yerrlist[i], fmt="xk",
                    label=ynames[i], capsize=5, elinewidth=2,
                    capthick=2)
    else:
        plt.errorbar(x, ylist, yerr=yerrlist, fmt="ok", label=ynames[0],
            capsize=5, elinewidth=2, capthick=2)
    if len(ynames) > 1:
        miny = 1
        maxy = 0
        maxyerr = 0
        for y in ylist:
            for i in y:
                if i > maxy and isnan(i) == False:
                    maxy = i
                elif i < miny and isnan(i) == False:
                    miny = i
        for yerr in yerrlist:
            for j in yerr:
                if j > maxyerr and isnan(j) == False:
                    maxyerr = j
        maxy += maxyerr
        if (miny - maxyerr) > 0:
            miny = miny - maxyerr
    else:
        miny = 1
        maxy = 0
        maxyerr = 0
        for y in ylist:
            if y > maxy and isnan(y) == False:
                maxy = y
            elif y < miny and isnan(y) == False:
                miny = y
        for yerr in yerrlist:
            if yerr > maxyerr and isnan(yerr) == False:
                maxyerr = yerr
        if (miny - maxyerr) <= 0:
            miny = 0
        else:
            miny = miny - maxyerr
    maxy = roundup(maxy, dec=1)
    miny = rounddown(miny, dec=1)
    plt.yticks((miny, (miny + maxy)/2, maxy), size=fsize)
    yname = ""
    if len(ynames) > 1:
        for i in ynames[:-1]:
            yname += i + "_"
        yname += ynames[-1]
    else:
        yname = ynames[0]
    if len(xnames) <= 20:
        plt.xticks(range(len(xnames)), xnames, rotation=rot, size=fsize)
    else:
        try:
            for i in range(len(xnames)):
                xnames[i] = int(xnames[i])
        except ValueError:
            for i in range(len(xnames)):
                xnames[i] = i+1
        bins = []
        for i in range(0, len(xnames), ceil(len(xnames)/12)):
            bins.append(xnames[i])
        bins_sep = []
        bins_index = 0
        for i in range(len(xnames)):
            if bins[bins_index] == xnames[i]:
                bins_sep.append(bins[bins_index])
                bins_index += 1
            else:
                bins_sep.append("")
            if bins_index == len(bins):
                while i <= len(xnames):
                    bins_sep.append("")
                    i += 1
                break
        plt.locator_params(axis="x", tight=True, nbins=12)
        plt.xticks(range(len(bins_sep)), bins_sep, rotation=rot,
            size=fsize)
    plt.xlabel(xname, fontsize=fsize)
    if len(ynames) > 1:
        plt.suptitle('Plotting ' + yname + ' with time', fontsize=fsize)
        plt.legend(loc='upper right', fontsize=fsize)
        plt.grid(True)
    plt.ylabel(yname, fontsize=fsize)
    plt.savefig(new_route + "/" + yname + "_vs_" +  xname +
        "_matplotlib.png")

# Plot values (it plots values X and Y given using pandas)
def plot_values(x, ylist, yerrlist,  xname="Variable1",
    ynames=["Variable2"], new_route="./data"):
    lists = []
    lists.append(x)
    if len(ylist) > 1 and type(ylist) == list:
        for i in ylist:
            lists.append(i)
    else:
        lists.append(ylist[0])
    if len(yerrlist) > 1 and type(yerrlist) == list:
        for i in yerrlist:
            lists.append(i)
    else:
        lists.append(yerrlist[0])
    names = []
    names.append(xname)
    if len(ynames) > 1 and type(ynames) == list:
        for i in ynames:
            names.append(i)
        for i in ynames:
            names.append(i + "_err")
    else:
        names.append(ynames[0])
        names.append(ynames[0] + "_err")
    dfobj = pd.DataFrame(lists, index = names, columns = range(len(x)))
    new_dfobj = dfobj.T
    if len(ynames) > 1 and type(ynames) == list:
        for i in range(len(ynames)):
            fig = new_dfobj.plot(style=".", x=xname, y=ynames[i],
                yerr=(ynames[i] + "_err")).get_figure()
            fig.savefig(new_route + "/" + ynames[i] + "_vs_" + xname +
                "_pandas.png")
    else:
        fig = new_dfobj.plot(style=".", x=xname, y=ynames[0],
            yerr=(ynames[0] + "_err")).get_figure()
        fig.savefig(new_route + "/" + ynames[0] + "_vs_" + xname +
            "_pandas.png")

# Plot Variability and Beta (it plots V values with error and Beta
# values with errors in separated plots, from Excel results obtained by
# cmplxcruncher)
def plot_V_B(files, min_times=0, min_elements=0, new_route="./data",
    old_route="./not_treated_data", V=[], V_err=[], B=[], B_err=[],
    times=[]):
    separate_excel(files, new_route=new_route, old_route=old_route)
    files = read_file(new_route, old_route)
    _DF_COLS = ['Elements', 'Times', 'AbsFreq_mean', 'AbsFreq_std',
        'AbsFreq_sum', 'V', 'V_err', 'beta', 'beta_err', 'R^2',
        'pcorr^2', 'model', 'xW_V', 'xW_V_err', 'xW_beta',
        'xW_beta_err', 'xW_R^2', 'xW_pcorr^2', 'xW_model', 'xW_V_stan',
        'xW_V_err_stan', 'xW_beta_stan', 'xW_beta_err_stan']
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        os.remove(old_route + "/" + filename)
        if list(df.columns)[1:] != _DF_COLS:
            print("Excel input must be cmplxcruncher Excel results")
            sys.exit()
        try:
            min_times = int(min_times)
            min_elements = int(min_elements)
        except ValueError:
            raise ValueError("Number of minimum times / elements must b\
e number >= 0")
        times.append(filename)
        if len(df.index) > 0:
            elements = False
            numtimes = False
            for i in range(len(df.columns)):
                if df.columns[i] == "Elements":
                    if float(df.iloc[0,i]) > min_elements:
                        elements = True
                if df.columns[i] == "Times":
                    if float(df.iloc[0,i]) > min_times:
                        numtimes = True
                if elements == True and numtimes == True:
                    if df.columns[i] == "xW_V":
                        V.append(float(df.iloc[0,i]))
                    if df.columns[i] == "xW_V_err":
                        V_err.append(float(df.iloc[0,i]))
                    if df.columns[i] == "xW_beta":
                        B.append(float(df.iloc[0,i]))
                    if df.columns[i] == "xW_beta_err":
                        B_err.append(float(df.iloc[0,i]))
            if elements == False or numtimes == False:
                V.append(np.NaN)
                V_err.append(np.NaN)
                B.append(np.NaN)
                B_err.append(np.NaN)
        else:
            V.append(np.NaN)
            V_err.append(np.NaN)
            B.append(np.NaN)
            B_err.append(np.NaN)
    for i in range(len(times)):
        time = times[i]
        time = re.sub("[M|m]onth[_1|1]0", "October", time)
        time = re.sub("[M|m]onth[_1|1]1", "November", time)
        time = re.sub("[M|m]onth[_1|1]2", "December", time)    
        time = re.sub("[M|m]onth[_1|1]", "January", time)
        time = re.sub("[M|m]onth[_2|2]", "February", time)
        time = re.sub("[M|m]onth[_3|3]", "March", time)
        time = re.sub("[M|m]onth[_4|4]", "April", time)
        time = re.sub("[M|m]onth[_5|5]", "May", time)
        time = re.sub("[M|m]onth[_6|6]", "June", time)
        time = re.sub("[M|m]onth[_7|7]", "July", time)
        time = re.sub("[M|m]onth[_8|8]", "August", time)
        time = re.sub("[M|m]onth[_9|9]", "September", time)
        time = re.sub("[Y|y]ear", "", time)
        time = re.sub("[D|d]ay", "", time)        
        times[i] = time.split("_")[-1].split(".")[0]
    if sorted(times) == ["April", "August", "December", "February",
        "January", "July", "June", "March", "May", "November",
        "October", "September"]:
        months = ["January", "February", "March", "April", "May", 
            "June", "July", "August", "September", "October",
            "November", "December"]
        for i in range(len(times)):
            for j in range(len(months)):
                if times[i] == months[j]:
                    times[i] = j+1
        V = [x for y, x in sorted(zip(times,V))]
        V_err = [x for y, x in sorted(zip(times,V_err))]
        B = [x for y, x in sorted(zip(times,B))]
        B_err = [x for y, x in sorted(zip(times,B_err))]
        times = months
    else:  ## if we do not have to order by months, we order by names
        V = [x for y, x in sorted(zip(times,V))]
        V_err = [x for y, x in sorted(zip(times,V_err))]
        B = [x for y, x in sorted(zip(times,B))]
        B_err = [x for y, x in sorted(zip(times,B_err))]
        times = sorted(times)
    matplot_values(times, [V, B], [V_err, B_err], xname="Time",
        ynames=["Variability", "Beta"], new_route=new_route)
    matplot_values(times, [V], [V_err], xname="Time",
        ynames=["Variability"], new_route=new_route)
    matplot_values(times, [B], [B_err], xname="Time", ynames=["Beta"],
        new_route=new_route)
    plot_values(times, [V, B], [V_err, B_err], xname="Tiempo",
        ynames=["Variability", "Beta"], new_route=new_route)
    plot_values(times, [V], [V_err], xname="Time",
        ynames=["Variability"], new_route=new_route)
    plot_values(times, [B], [B_err], xname="Time", ynames=["Beta"],
        new_route=new_route)

# Round up 2 (round a given number up, with positions instead of
# decimals given)
def roundup2(x, pos=0):
    y = 1
    if pos > 0:
        for i in range(pos):
            y = y * 10
        x = ceil(x / y) * y
    return(x)

# Scientific (convert a given number to scientific format, with decimal
# positions given)
def scientific(x, pos):
    return('%.0E' % x)

# Plot Frequencies (it plots absolute frequencies values with error)
def plot_freqs(files, new_route="./data",
    old_route="./not_treated_data", freqs=[], times=[], errors=[]):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        df = df.drop(df.columns[0], axis=1)
        print(df)
        df = df *30
        print(df)
        freq = df.to_numpy().sum() / len(df.columns)
        error = df.to_numpy().std()
        errors.append(error)
        freqs.append(freq)
        times.append(filename.split("_")[-1].split(".")[0])
    if sorted(times) == ["April", "August", "December", "February",
        "January", "July", "June", "March", "May", "November",
        "October", "September"]:
        months = ["January", "February", "March", "April", "May", 
            "June", "July", "August", "September", "October",
            "November", "December"]
        for i in range(len(times)):
            for j in range(len(months)):
                if times[i] == months[j]:
                    times[i] = j+1
        freqs = [x for y, x in sorted(zip(times,freqs))]
        errors = [x for y, x in sorted(zip(times,errors))]
        times = months
    else:  ## if we do not have to order by months, we order by names
        freqs = [x for y, x in sorted(zip(times,freqs))]
        errors = [x for y, x in sorted(zip(times,errors))]
        times = sorted(times)
    df = pd.DataFrame([times, freqs, errors],
        index=["Times", "Frequencies", "Deviations"],
        columns=range(len(times)))
    new_df = df.T
    maxy = roundup2(max(freqs) + max(errors),
        pos=len(str(ceil(max(freqs))))-1)
    miny = 0
    ticks = [miny, (miny + maxy)/3, 2*(miny + maxy)/3, maxy]
    plt = new_df.plot(kind="bar", x="Times", y="Frequencies",
        yerr="Deviations", yticks=ticks, legend=False, figsize=(10,10),
        capsize=5, fontsize=10)
    plt.tick_params(axis="x", labelsize=10, rotation=30)
    plt.set_xlabel("Months")
    plt.set_ylabel("Absolute Frequency")
    plt.set_title("Histogram of monthly absolute frequencies")
    scientific_formatter = FuncFormatter(scientific)
    plt.yaxis.set_major_formatter(scientific_formatter)
    fig = plt.get_figure()
    fig.savefig(new_route + "/Freqs_vs_Time.png")

# Species Frequencies (it obtains absolute frequencies by species)
def species_freqs(files, new_route="./data",
    old_route="./not_treated_data"):
    for filename in files:
        sheets = pd.read_excel(old_route + "/" + filename,
            sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        names = list(df.loc[:, df.columns.values[0]])
        df = df.drop(df.columns[0], axis=1)
        new_df = df.sum(axis=1) / len(df.columns)
        new_df.to_excel(new_route + "/" +
            filename.split("/")[-1].split(".")[0] + "_freqs_species." +
            filename.split(".")[-1], sheet_name, header=True,
            index=False)

## Main script
def main():
    ## Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Interactive preprocessing for Complex systems anal\
yzer Python3 code (new Excel generated in data folder, old Excels in da\
ta folder to not_treated_data folder)',
        epilog='interactive Excel preprocessing for python tool cmplxcr\
uncher, developed by J. M. Marti (DLS team) - by Torres A - ' + _verdata
    )
    parser.add_argument(
        '-d', '--data',
        action='store',
        default=_DEFAULT_ROUTE,
        help=('directory path with data files (if not present, \'' +
            _DEFAULT_ROUTE + '\' will be tried)')
    )
    parser.add_argument(
        '-m', '--method',
        action='store',
        default="preprocessing",
        choices=["preprocessing", "plotting"],
        help=('Use this script to preprocess Excel matrix before execut\
ing cmplxcruncher, or to plot results (default: preprocessing)')
    )
    parser.add_argument(
        '-r', '--rawdata',
        action='store',
        default=_DEFAULT_OLD_ROUTE,
        help=('path of the directory to move raw data relative to the d\
ata directory (if not present, \'' + _DEFAULT_OLD_ROUTE + '\' will be t\
ried)')
    )

    ## Parse arguments
    flgs = parser.parse_args()
    datapath = flgs.data
    method = flgs.method
    _DEFAULT_FREQLIM_SUFFIXDIR = ''
    if (flgs.data == _DEFAULT_ROUTE):
        route = flgs.data + _DEFAULT_FREQLIM_SUFFIXDIR + '/'
    else:
        route = _DEFAULT_ROUTE + "/"
    if (flgs.rawdata == _DEFAULT_OLD_ROUTE):
        old_route = flgs.rawdata + _DEFAULT_FREQLIM_SUFFIXDIR + '/'
    else:
        old_route = _DEFAULT_OLD_ROUTE + "/"
    if (ChkDir(old_route)):
        old_route = ''
        print('Warning: the data root directory will be tried instead, \
for writing not treated data')

    ## Execute functions
    if method == "preprocessing":
        print("\nPreprocessing:\n")
        separate = input("Do you want to separate your Excel sheets? (o\
nly one sheet for each Excel file) (y/n): ")
        while separate not in ["y", "Y", "n", "N"]:
            separate = input("Do you want to separate your Excel sheets\
? (only one sheet for each Excel file) (y/n): ")
        if separate in ["Y", "y"]:
            files = read_file(route, old_route)
            separate_excel(files, new_route=route, old_route=old_route)
        
        transpose = input("Do you want to transpose your data? (species\
 as rows, dates as colums) (y/n): ")
        while transpose not in ["y", "Y", "n", "N"]:
            transpose = input("Do you want to transpose your data? (spe\
cies as rows, dates as colums) (y/n): ")
        if transpose in ["Y", "y"]:
            files = read_file(route, old_route)
            transpose_file(files, new_route=route, old_route=old_route)
        
        fill = input("Do you want to fill NA in your data? (replace NA \
values with 0s) (y/n): ")
        while fill not in ["y", "Y", "n", "N"]:
            fill = input("Do you want to fill NA in your data? (replace\
 NA values with 0s) (y/n): ")
        if fill in ["Y", "y"]:
            files = read_file(route, old_route)
            fill_file(files, new_route=route, old_route=old_route)
        
        freqs = input("Do you want to calculate relative frequencies? (\
y/n): ")
        while freqs not in ["y", "Y", "n", "N"]:
            freqs = input("Do you want to calculate relative frequencie\
s? (y/n): ")
        if freqs in ["Y", "y"]:
            files = read_file(route, old_route)
            relative_freqs(files, new_route=route, old_route=old_route)
        
        merge = input("Do you want to merge the two dataframes? (y/n): \
")
        while merge not in ["y", "Y", "n", "N"]:
            merge = input("Do you want to merge the two dataframes? (y/\
n): ")
        if merge in ["Y", "y"]:
            files = read_file(route, old_route)
            if len(files) == 2:
                merge_2files(files[0], files[1], new_route=route,
                    old_route=old_route)
            else:
                print("You need only two files, with the same extension\
 (Excel or txt/tsv), to merge.")
                sys.exit()
        
        dates = input("Do you want to group data by dates? (y/n): ")
        while dates not in ["y", "Y", "n", "N"]:
            dates = input("Do you want to group data by dates? (y/n): ")
        if dates in ["Y", "y"]:
            files = read_file(route, old_route)
            filenames = listdir(old_route)
            days = input("Select the days you want to select, separated\
 by commas, without spaces (0-31, <ENTER> for no separation, 'a' for al\
l days): ")
            if days == "":
                days = []
            elif days != "" and days != "a":
                days = days.split(",")
            months = input("Select the months you want to select, separ\
ated by commas, without spaces (0-12, <ENTER> for no separation, 'a' fo\
r all months): ")
            if months == "":
                months = []
            elif months != "" and months != "a":
                months = months.split(",")
            years = input("Select the years you want to select, separat\
ed by commas, without spaces (<ENTER> for no separation, 'a' for all ye\
ars): ")
            if years == "":
                years = []
            elif years == "a":
                year1 = input("Introduce the lowest year of your data: \
")
                year2 = input("Introduce the higgest year of your data:\
 ")
                try:
                    year1 = int(year1)
                    year2 = int(year2)
                except ValueError:
                    raise ValueError("Years must be numbers between 0 a\
nd 9999 inclusive.")
            else:
                years = years.split(",")
            if days != "a" and months != "a" and years != "a":
                select_by_dates(files, new_route=route, days=days,
                    months=months, years=years, old_route=old_route)
            elif days == "a" and months != "a" and years != "a":
                for i in range(31):
                    select_by_dates(files, new_route=route, days=[i+1],
                        months=months, years=years, old_route=old_route)
                    new_files = read_file(route, old_route,
                        new_filename=(
                        filenames[0].split("/")[-1].split(".")[0]) +
                        "_day" + str(i+1) + "." +
                        filenames[0].split("/")[-1].split(".")[-1])
            elif days == "a" and months == "a" and years != "a":
                for i in range(31):
                    for j in range(12):
                        select_by_dates(files, new_route=route,
                            days=[i+1], months=[j+1], years=years,
                            old_route=old_route)
                        new_files = read_file(route, old_route,
                            new_filename=(
                            filenames[0].split("/")[-1].split(".")[0]) +
                            "_day" + str(i+1) + "_month" + str(j+1) +
                            "." +
                            filenames[0].split("/")[-1].split(".")[-1])
            elif days == "a" and months == "a" and years == "a":
                for i in range(31):
                    for j in range(12):
                        for k in range(year1, year2+1):
                            select_by_dates(files, new_route=route,
                                days=[i+1], months=[j+1], years=[k],
                                old_route=old_route)
                            new_files = read_file(route, old_route,
                                new_filename=(
                                filenames[0].split("/")[-1].split(".\
")[0]) + "_day" + str(i+1) + "_month" + str(j+1) + "_year" + str(k) +
                                "." + filenames[0].split("/\
")[-1].split(".")[-1])
            elif days != "a" and months == "a" and years != "a":
                for j in range(12):
                    select_by_dates(files, new_route=route, days=days,
                        months=[j+1], years=years, old_route=old_route)
                    new_files = read_file(route, old_route,
                        new_filename=(
                        filenames[0].split("/")[-1].split(".")[0]) +
                        "_month" + str(j+1) + "." +
                        filenames[0].split("/")[-1].split(".")[-1])
            elif days != "a" and months == "a" and years == "a":
                for j in range(12):
                    for k in range(year1, year2+1):
                        select_by_dates(files, new_route=route,
                            days=days, months=[j+1], years=[k],
                            old_route=old_route)
                        new_files = read_file(route, old_route,
                            new_filename=(
                            filenames[0].split("/")[-1].split(".")[0]) +
                            "_month" + str(j+1) + "_year" + str(k) +
                            "." +
                            filenames[0].split("/")[-1].split(".")[-1])
            elif days != "a" and months != "a" and years == "a":
                for k in range(year1, year2+1):
                    select_by_dates(files, new_route=route, days=days,
                        months=months, years=[k], old_route=old_route)
                    new_files = read_file(route, old_route,
                        new_filename=(
                        filenames[0].split("/")[-1].split(".")[0]) +
                        "_year" + str(k) + "." +
                        filenames[0].split("/")[-1].split(".")[-1])

        species = input("Do you want to group data by species? (y/n): ")
        while species not in ["y", "Y", "n", "N"]:
            species = input("Do you want to group data by species? (y/n\
): ")
        if species in ["Y", "y"]:
            files = read_file(route, old_route)
            species = input("Select the rows of species you want to sel\
ect (column names is row 0), separated by commas, without spaces ('g' f\
or a group of species): ")
            if species == "g":
                species1 = input("Introduce the lowest row of your data\
: ")
                species2 = input("Introduce the higgest row of your dat\
a: ")
                try:
                    species1 = int(species1)
                    species2 = int(species2)
                except ValueError:
                    raise ValueError("Rows of species must be numbers b\
etween 0 and 9999 inclusive.")
                species = []
                for i in range(species1, species2):
                    species.append(i)
            else:
                species = species.split(",")
            select_by_species(files, species=species, new_route=route,
                old_route=old_route)

        zeros = input("Do you want to group data by ZRF? (y/n): ")
        while zeros not in ["y", "Y", "n", "N"]:
            zeros = input("Do you want to group data by ZRF? (y/n): ")
        if zeros in ["Y", "y"]:
            files = read_file(route, old_route)
            core = input("Percent of core for the data (for example, 5 \
for data with <5 percent of zeros) (<ENTER> for nothing) (cmplxcruncher\
 uses 20 by default): ")
            if core != "":
                select_by_ZRF(files, core=core, new_route=route,
                    old_route=old_route)
            else:
                tail = input("Percent of tail for the data (for example\
, 95 for data with >=95 percent of zeros) (<ENTER> for nothing) (cmplxc\
runcher uses 60 by default for >4 time series, 50 for 4 time series, an\
d 40 for <4 time series): ")
                if tail != "":
                    select_by_ZRF(files, tail=tail, new_route=route,
                        old_route=old_route)

    elif method == "plotting":
        print("\nAnalyze results:\n")
        plotVB = input("Do you want to plot V and Beta results of your \
data ? (V and Beta vs Time or Species) (cmplxcruncher Excel results mus\
t be given) (y/n): ")
        while plotVB not in ["y", "Y", "n", "N"]:
            plotVB = input("Do you want to plot V and Beta results of y\
our data? (V and Beta vs Time or Species) (cmplxcruncher Excel results \
must be given) (y/n): ")
        if plotVB in ["Y", "y"]:
            files = read_file(route, old_route)
            plot_V_B(files, new_route=route, old_route=old_route)
            input("\nPress <ENTER> to finish.\n")
            sys.exit()
        
        plotfreqs = input("Do you want to plot absolute frequencies of \
your data? (quantity vs Time or Species) (cmplxcruncher input Excel mus\
t be given) (y/n): ")
        while plotfreqs not in ["y", "Y", "n", "N"]:
            plotfreqs = input("Do you want to plot absolute frequencies\
 of your data? (Freqeuncy vs Time or Species) (cmplxcruncher input Exce\
l must be given) (y/n): ")
        if plotfreqs in ["Y", "y"]:
            files = read_file(route, old_route)
            plot_freqs(files, new_route=route, old_route=old_route)
            input("\nPress <ENTER> to finish.\n")
            sys.exit()

        zrf = input("Do you want to calculate ZRF values? (cmplxcrunche\
r input Excel must be given) (y/n): ")
        while zrf not in ["y", "Y", "n", "N"]:
            zrf = input("Do you want to calculate ZRF values? (cmplxcru\
ncher input Excel must be given) (y/n): ")
        if zrf in ["Y", "y"]:
            files = read_file(route, old_route)
            calculate_ZRF(files, new_route=route, old_route=old_route)

        spg = input("Do you want to obtain absolute freqs of species? (\
cmplxcruncher input Excel must be given) (y/n): ")
        while spg not in ["y", "Y", "n", "N"]:
            spg = input("Do you want to obtain absolute freqs of specie\
s? (cmplxcruncher input Excel must be given) (y/n): ")
        if spg in ["Y", "y"]:
            files = read_file(route, old_route)
            species_freqs(files, new_route=route, old_route=old_route)

if __name__ == "__main__":
    main()

input("\nPress <ENTER> to finish.\n")