#!/usr/bin/env python3

"""
Script to obtain Krona results from Excel counts data.

counts_to_krona.py

Script purpose: the script obtains Krona results from Excel counts data,
with taxon names or tax ids as rows and samples as columns
"""

__version__ = '1.0'
_verdata = 'Jun 2020'

## Import libraries
_DEFAULT_ROUTE = 'data'
_DEFAULT_RESULTS = 'results'
_DEFAULT_TAXONOMY = 'taxonomy'
_DEFAULT_FREQLIM_SUFFIXDIR = ''
_version_pandas = '0.16.0'
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
except ImportError:
    raise ImportError("This module requires pandas")
else:
    if not compare_versions(pd.__version__, _version_pandas):
        raise ImportError('pandas %s or later is required; you have %s'
            % (_version_pandas, matplotlib.__version__))

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
        files.append(filename)
    if len(files) == 0:
        print("\nFolder '" + route + "' is empty\n")
        sys.exit()
    return(files)

# Represent Excel counts matrix with Krona (it obtains as many Krona
# results as samples present in Excel matrix, a Krona summary, and a
# text file with unclassified taxa)
def excel_to_krona(route="data", results="results", taxonomy="taxonomy",
    update=False, rows="taxname", maxcounts = 10000):
    files = read_file(route)
    if update:
        os.system("mkdir -p ~/krona/taxonomy")
        os.system("git clone https://github.com/marbl/Krona.git")
        os.system("Krona/KronaTools/updateTaxonomy.sh ~/krona/taxonomy")
        os.system("gzip -d taxonomy.tab.gz")
        os.system("mv taxonomy.tab ~/krona/taxonomy")
    for filename in files:
        sheets = pd.read_excel(route + "/" + filename, sheet_name=None)
        sheet_name, df = next(iter(sheets.items()))
        rownames = df.index.values
        colnames = df.columns.values
        columns = []
        nrow = df.shape[0]
        ncol = df.shape[1]
        for row in range(nrow):
            for col in range(ncol):
                if row == 0:
                    columns.append([])
                columns[col].append(df.loc[rownames[row],colnames[col]])
        taxa = columns[0]
        counts = columns[1:]
        if rows == "taxid":
            maxcounts = maxcounts - len(taxa)
            for i in range(len(counts)):
                sumcounts = 0
                factor = 1
                for j in range(len(counts[i])):
                    if j in index_selected:
                        sumcounts += counts[i][j]
                if sumcounts > maxcounts:
                    factor = sumcounts / maxcounts
                new_filename = (results + "/" + filename.split(".")[0] +
                    "_" + str(colnames[i+1]).split(" ")[0] + ".krona")
                new_file = open(new_filename, "w")
                for j in range(len(counts[i])):
                    newcount = int(round(float(counts[i][j]) / factor))
                    newtaxid = str(taxa[j])
                    if newcount > 0:
                        for k in range(newcount):
                            new_file.write("read_" + k + "_" + taxa[j] +
                                "\t" + newtaxid + "\n")
                    else:
                        not_selected_species.append(taxa[k])
                new_file.close()
                new_filename2 = (results + "/" + filename.split(".")[0] +
                    str(colnames[i+1]).split(" ")[0] +
                    "_species_discarded.txt")
                new_file2 = open(new_filename2, "w")
                if factor > 1:
                    print("Too many reads (> 10.000 reads) in " +
                    new_filename + ", taxa with less than " +
                    str(int(round(factor))) + " reads won't be plot " +
                    "(not selected species in " + new_filename2 + ")")
                    for specie in not_selected_species:
                        new_file2.write(specie + "\n")
                new_file2.close()
                os.system("ktImportTaxonomy " + new_filename +
                    " -o " + new_filename + ".html")
                #os.system("google-chrome " + new_filename + ".html")
        elif rows == "taxname":
            tax_files = read_file(taxonomy)
            if taxonomy[-1] == "/":
                taxonomy = taxonomy[:-1]
            if "names.dmp" in tax_files:
                names = open(taxonomy + "/names.dmp", "r").readlines()
            else:
                print("names.dmp file not in taxonomy specified folder")
                sys.exit()
            taxids = []
            taxnames = []
            taxsins = []
            for line in names:
                taxids.append(line.split("\t|\t")[0].strip())
                taxnames.append(line.split("\t|\t")[1].strip().lower())
                taxsins.append(line.split("\t|\t")[2].strip().lower())
            taxids_selected = []
            index_selected = []
            no_taxnames = []
            replaced_taxnames = []
            for i in range(len(taxa)):
                if taxa[i].lower() in taxnames:
                    j = taxnames.index(taxa[i].lower())
                    taxids_selected.append(taxids[j])
                    index_selected.append(i)
                elif taxa[i].lower() in taxsins:
                    j = taxsins.index(taxa[i].lower())
                    taxids_selected.append(taxids[j])
                    index_selected.append(i)
                else:
                    if taxa[i] not in no_taxnames:
                        no_taxnames.append(taxa[i])
                    if len(taxa[i].split(" ")) == 2:
                        old_taxa = taxa[i]
                        taxa[i] = ("unclassified " +
                            taxa[i].split(" ")[0])
                        if taxa[i].lower() in taxnames:
                            j = taxnames.index(taxa[i].lower())
                            taxids_selected.append(taxids[j])
                            index_selected.append(i)
                            replaced_taxnames.append(old_taxa)
                        elif taxa[i].lower() in taxsins:
                            j = taxsins.index(taxa[i].lower())
                            taxids_selected.append(taxids[j])
                            index_selected.append(i)
                            replaced_taxnames.append(old_taxa)
            if len(no_taxnames) > 0:
                print("No present taxa are show in output file " +
                    results + "/not_present_taxa.txt")
                print("If taxa is specie (genera_name + specie_name), '\
unclassified genera_name' will be tried")
                no_taxnames_file = open(results +
                    "/not_present_taxa.txt", "w")
                for notax in no_taxnames:
                    if notax in replaced_taxnames:
                        no_taxnames_file.write(notax +
                            " (replaced by 'unclassified " +
                            notax.split(" ")[0] + "')\n")
                    else:
                        no_taxnames_file.write(notax + "\n")
                no_taxnames_file.close()
            maxcounts = maxcounts - len(taxids_selected)
            for i in range(len(counts)):
                sumcounts = 0
                factor = 1
                for j in range(len(counts[i])):
                    if j in index_selected:
                        sumcounts += counts[i][j]
                if sumcounts > maxcounts:
                    factor = sumcounts / maxcounts
                not_selected_species = []
                new_filename = (results + "/" + filename.split(".")[0] +
                    "_" + str(colnames[i+1]).split(" ")[0] + ".krona")
                new_file = open(new_filename, "w")
                for j in range(len(counts[i])):
                    if j in index_selected:
                        k = index_selected.index(j)
                        newcount = int(round(
                            float(counts[i][j]) / factor))
                        newtaxid = str(taxids_selected[k])
                        if newcount > 0:
                            for l in range(newcount):
                                new_file.write("read_" + str(l+1) +
                                    "_" + taxa[j] + "\t" + newtaxid +
                                    "\n")
                        else:
                            not_selected_species.append(taxa[k])
                new_file.close()
                new_filename2 = (results + "/" + filename.split(".")[0] +
                    str(colnames[i+1]).split(" ")[0] +
                    "_species_discarded.txt")
                new_file2 = open(new_filename2, "w")
                if factor > 1:
                    print("Too many reads (> 10.000 reads) in " +
                    new_filename + ", taxa with less than " +
                    str(int(round(factor))) + " reads won't be plot " +
                    "(not selected species in " + new_filename2 + ")")
                    for specie in not_selected_species:
                        new_file2.write(specie + "\n")
                new_file2.close()
                os.system("ktImportTaxonomy " + new_filename + " -o " +
                    new_filename + ".html")
            #os.system("google-chrome " new_filename + ".html")
        krona_results = read_file(results)
        all_results = ""
        for filename2 in krona_results:
            if filename2.split(".")[-1] == "krona":
                all_results += results + "/" + filename2 + " "
        os.system("ktImportTaxonomy " + all_results + "-o " + results +
            "/" + filename.split(".")[0] + "_krona_results.html")

## Main script
def main():
    ## Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Script to obtain Krona results from Excel counts d\
ata, with taxon names or taxa ids as rows and samples as columns',
        epilog=('Script to obtain Krona results from Excel counts data \
, with taxon names or taxa ids as rows and samples as columns - by Torr\
es A - ' + _verdata)
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
        '-t', '--taxonomy',
        action='store',
        default=_DEFAULT_TAXONOMY,
        help=('path of the directory with taxonomy "names.dmp" file , r\
elative to the data directory (if not present, \'' + _DEFAULT_TAXONOMY +
'\' will be tried)')
    )
    parser.add_argument(
        '-D', '--download',
        action='store',
        choices=["yes", "no"],
        default="no",
        help=('download taxonomy database, instead of using -t paramete\
r (default: no)')
    )
    parser.add_argument(
        '-R', '--rowstype',
        action='store',
        choices=["taxname", "taxid"],
        default="taxname",
        help=('Type of data used as row names, taxa names complete name\
s or taxa ids (default: taxname)')
    )
    parser.add_argument(
        '-u', '--updatetaxonomy',
        action='store',
        choices=["yes", "no"],
        default="no",
        help=('Update or not Taxonomy database for Krona result (defaul\
t: no)')
    )

    ## Parse arguments
    flgs = parser.parse_args()
    download = flgs.download
    rows = flgs.rowstype
    update = flgs.updatetaxonomy
    _DEFAULT_FREQLIM_SUFFIXDIR = ''
    if (flgs.data == _DEFAULT_ROUTE):
        route = flgs.data + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        route = flgs.route
    if (flgs.results == _DEFAULT_RESULTS):
        results = flgs.results + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        results = flgs.results
    if (ChkDir(results)):
        results = ''
        print('Warning: the data root directory will be tried instead, \
for creating results folder')
    if (ChkDir(route)):
        route = ''
        print('Warning: the data root directory will be tried instead, \
for creating data folder')
    if download == "yes":
        os.system("mkdir taxonomy")
        os.system("cd taxonomy")
        os.system("wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump\
.tar.gz")
        os.system("tar xvzf taxdump.tar.gz")
        os.system("cd ..")
        taxonomy = _DEFAULT_TAXONOMY
    else:
        if (flgs.taxonomy == _DEFAULT_TAXONOMY):
            taxonomy = flgs.taxonomy + _DEFAULT_FREQLIM_SUFFIXDIR
        else:
            taxonomy = flgs.taxonomy
    if update == "yes":
        update = True
    else:
        update = False

    ## Execute functions
    excel_to_krona(route=route, results=results, taxonomy=taxonomy,
        update=update, rows=rows)

if __name__ == "__main__":
    main()
