#!/usr/bin/env python3

"""
Script to obtain an Excel matrix from kreport, including or not only
microalgae results, with selected taxa level, and normalize counts by
genome size

centrifuge_normalize.py

Script purpose: the script filtrates Centrifuge kreports for
cmplxcruncher analysis. It obtains filtered and/or genome size
normalized counts matrix in Excel 
"""

__version__ = '1.0'
_verdata = 'Jun 2020'

## Import libraries
_DEFAULT_KREPORT = 'kreports'
_DEFAULT_TSV = 'tsv'
_DEFAULT_OUT = 'tab'
_DEFAULT_RESULTS = 'results'
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
# Read files from a given folder (it returns filenames)
def read_file(route="data"):
    files = []
    for filename in listdir(route):
        files.append(filename)
    if len(files) == 0:
        print("There are no files in folder '" + route + "'\n")
        sys.exit()
    return(files)

# Select taxa ids (it returns taxa ids, taxa names and number of reads
# from a kreport file and a given taxa)
def select_taxids(kreport="kreport.txt", taxa="species",
    taxon=["unclassified"]):
    kreport_file = open(kreport, "r").readlines()
    print("Unclassified reads in file '" + kreport + "': " +
        str(kreport_file[0].split("\t")[1].strip()) + " (" +
        str(kreport_file[0].split("\t")[0].strip()) + " %)")
    tax_letters = ["S", "G", "F", "O", "C", "P", "K", "D", "U"]
    save = False
    taxa_letter = "U"
    species = []
    species_taxids = []
    reads = []
    taxids = []
    for i in range(len(kreport_file)):
        read = int(kreport_file[i].split("\t")[1].strip())
        letter = kreport_file[i].split("\t")[3].strip()
        taxid = kreport_file[i].split("\t")[4].strip()
        specie = kreport_file[i].split("\t")[5].strip().lower()
        if letter == taxa_letter:
            save = False
        if specie in taxon and letter in tax_letters:
            taxa_letter = letter
            save = True
        if save == True:
            reads.append(read)
            taxids.append(taxid)
            if taxa == letter:
                species.append(specie)
                species_taxids.append(taxid)
    return(taxids, species_taxids, species, reads)

# Filter reads (it returns number of reads filtered for selected taxa
# ids, from a given Centrifuge output file)
def filter_reads(out="out.tab", length=1, mhl=1, quality=1, taxids=[],
    results="results"):
    out_file = open(out, "r").readlines()
    out_results = open(results + "/" +
        out.split("/")[-1].split(".")[0] + "_filtered." +
        out.split(".")[-1], "w")
    out_results.write(out_file[0])
    counts = [0] * len(taxids)
    for line in out_file[1:]:
        cols = line.split("\t")
        file_taxid = cols[2]
        file_quality = int(cols[3])
        file_mhl = int(cols[5])
        file_length = int(cols[6])
        file_matches = int(cols[7])
        if (file_matches == 1 and file_length >= length and file_quality
            >= quality and file_mhl >= mhl and file_taxid in taxids):
            out_results.write(line)
            ind = taxids.index(file_taxid)
            counts[ind] += 1
    out_results.close()
    return(counts)

# Normalize reads (it returns number of reads normalized by genome size
# for selected taxa ids, from a given Centrifuge tsv file)
def normalize_reads(tsv="tsv.tsv", taxids=[], reads=[], rounded=False):
    tsv_file = open(tsv, "r").readlines()
    counts = [0] * len(taxids)
    for line in tsv_file[1:]:
        cols = line.split("\t")
        file_taxid = cols[1]
        file_genome = int(cols[3])
        if file_taxid in taxids and file_genome > 0:
            ind = taxids.index(file_taxid)
            if len(reads) == 0:
                file_reads = int(cols[5])
            else:
                file_reads = reads[ind]
            if rounded:
                count = int(round((file_reads / file_genome) * 1000))
            else:
                count = (file_reads / file_genome) * 1000
            counts[ind] = count
    return(counts)

# Convert Cenrtifuge results to an Excel matrix for cmplxcruncher (it
# returns an Excel matrix with raw reads for selected taxa ids, an Excel
# matrix with filtered reads, and an Excel matrix with normalized reads)
def centrifuge_to_cmplxcruncher(kreport="kreports", tsv="tsv",
    out="tab", results="results", length=1, mhl=1, quality=1,
    taxa="species", taxon="unclassified"):
    taxa_list1 = ["species", "genus", "family", "order", "class",
        "phylum", "kingdom", "domain", "unclassified"]
    taxa_list2 = ["S", "G", "F", "O", "C", "P", "K", "D", "U"]
    for i in range(len(taxa_list1)):
        if taxa_list1[i] == taxa:
            taxa = taxa_list2[i]
    taxon = taxon.lower().strip().split(",")
    kreports = read_file(kreport)
    tsvs = read_file(tsv)
    outs = read_file(out)
    filenames = []
    extensions = []
    for i in kreports:
        for j in tsvs:
            for k in outs:
                if (i.split(".")[0] == k.split(".")[0] and
                    i.split(".")[0] == j.split(".")[0]):
                    filenames.append(i.split(".")[0])
                    if len(extensions) < 3:
                        extensions.append(i.split(".")[-1])
                        extensions.append(j.split(".")[-1])
                        extensions.append(k.split(".")[-1])
    numfile = 0
    for filename in filenames:
        numfile += 1
        kreport_name = kreport + "/" + filename + "." + extensions[0]
        tsv_name = tsv + "/" + filename + "." + extensions[1]
        out_name = out + "/" + filename + "." + extensions[2]
        print("Selecting tax ids for " + filename + "...")
        taxids, species_taxids, species, reads = select_taxids(
            kreport=kreport_name, taxa=taxa, taxon=taxon)
        print("Filtering reads for " + filename + "...")
        filtered_reads = filter_reads(out=out_name, length=length,
            mhl=mhl, quality=quality, taxids=taxids, results=results)
        print("Normalizing reads for " + filename + "...")
        normalized_reads = normalize_reads(tsv=tsv_name, taxids=taxids,
            reads=filtered_reads)
        print("Saving results for " + filename + "...")
        reads_species = [0]*len(species)
        filtered_reads_species = [0]*len(species)
        normalized_reads_species = [0]*len(species)
        for i in range(len(species_taxids)):
            if species_taxids[i] in taxids:
                ind = taxids.index(species_taxids[i])
                reads_species[i] = reads[ind]
                filtered_reads_species[i] = filtered_reads[ind]
                normalized_reads_species[i] = normalized_reads[ind]
        data = {filename: reads_species}
        data_filt = {filename: filtered_reads_species}
        data_filt_norm = {filename: normalized_reads_species}
        df = pd.DataFrame (data, index=species, columns=[filename])
        df_filt = pd.DataFrame (data_filt, index=species, columns=[
            filename])
        df_filt_norm = pd.DataFrame (data_filt_norm, index=species,
            columns=[filename])
        if numfile == 1:
            df1 = df
            df2 = df_filt
            df3 = df_filt_norm
        else:
            df1 = pd.concat([df1, df], axis=1, sort=True)
            df2 = pd.concat([df2, df_filt], axis=1, sort=True)
            df3 = pd.concat([df3, df_filt_norm], axis=1, sort=True)
    df1 = df1.fillna(0)
    df2 = df2.fillna(0)
    df3 = df3.fillna(0)
    df1.to_excel(results + "/results_raw_reads.xlsx",
        sheet_name="results", header=True, index=True)
    df2.to_excel(results + "/results_filtered_reads.xlsx",
        sheet_name="results", header=True, index=True)
    df3.to_excel(results + "/results_normalized_reads.xlsx",
        sheet_name="results", header=True, index=True)

## Main script
def main():
    ## Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Script to obtain an Excel matrix from Centrifuge r\
esults, filtering and normalizing by genome size. Input filenames must \
be the same (except extension)',
        epilog=('Script to obtain an Excel matrix from Centrifuge resul\
ts, filtering and normalizing by genome size - by Torres A - ' +
            _verdata)
    )
    parser.add_argument(
        '-k', '--kreport',
        action='store',
        default=_DEFAULT_KREPORT,
        help=('Path of kreport files directory (if not present, \'' +
            _DEFAULT_KREPORT + '\' will be tried)')
    )
    parser.add_argument(
        '-t', '--tsv',
        action='store',
        default=_DEFAULT_TSV,
        help=('Path of tab separated values (tsv) files directory, obta\
ined from Centrifuge (if not present, \'' + _DEFAULT_TSV + '\' will be \
tried).')
    )
    parser.add_argument(
        '-o', '--out',
        action='store',
        default=_DEFAULT_OUT,
        help=('Path of Centrifuge output files (with tab separated value\
s) directory, obtained from Centrifuge (if not present, \'' +
            _DEFAULT_OUT + '\' will be tried)')
    )
    parser.add_argument(
        '-r', '--results',
        action='store',
        default=_DEFAULT_RESULTS,
        help=('Path of the directory for saving results, relative to th\
e data directory (if not present, \'' + _DEFAULT_RESULTS + '\' will be \
tried)')
    )
    parser.add_argument(
        '-l', '--length',
        action='store',
        default=1,
        help=('Minimum reads length for filtering (default: 1). As an e\
xample for 100 bp reads, 50-60 bp filter is used in Centrifuge article,\
 and 60 bp in Recentrifuge')
    )
    parser.add_argument(
        '-m', '--mhl',
        action='store',
        default=1,
        help=('Minimum hits lentgh (MHL) for filtering (default: 1). As\
 an example, 75 MHL value is used in Recentrifuge.')
    )
    parser.add_argument(
        '-q', '--quality',
        action='store',
        default=1,
        help=('Minimum quality reads (score) for filtering (default: 1)\
. As an example, 300 minimum score is used in Centrifuge article.')
    )
    parser.add_argument(
        '-s', '--selectedtaxa',
        action='store',
        default="species",
        choices=["species", "genus", "family", "order", "class",
        "phylum", "kingdom", "domain", "unclassified"],
        help=('Taxonomy level for Excel results (default: species)')
    )
    parser.add_argument(
        '-a', '--algaetaxa',
        action='store',
        choices=["yes", "no"],
        default="yes",
        help=('Use only phylums, domains and classes of microalgae as t\
axon groups, under which the indicated taxon level must be selected (de\
fault: yes). The taxon groups are the following: cyanobacteria, Bacilla\
riophyta, Bolidophyceae, Charophyceae, Chlorarachniophyceae, Chlorokybo\
phyceae, Chlorophyta, Chrysophyceae, Coleochaetophyceae, Cryptophyceae,\
 Dictyochophyceae, Dinophyceae, Euglenida, Eustigmatophyceae, Glaucocys\
tophyceae, Haptophyta, Klebsormidiophyceae, Mesostigmatophyceae, Aurear\
enophyceae, Xanthophyceae, Chrysoparadoxa, Phaeothamniophyceae, Raphido\
phyceae, Bangiophyceae, Compsopogonophyceae, Rhodellophyceae, Stylonema\
tophyceae, Synchromophyceae, Synurophyceae, Xanthophyceae, Zygnemophyce\
ae')
    )

    ## Parse arguments
    flgs = parser.parse_args()
    kreport = flgs.kreport
    tsv = flgs.tsv
    out = flgs.out
    length = flgs.length
    mhl = flgs.mhl
    quality = flgs.quality
    taxa = flgs.selectedtaxa.lower()
    algae = flgs.algaetaxa.lower()
    if (flgs.kreport == _DEFAULT_KREPORT):
        kreport = flgs.kreport + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        kreport = flgs.kreport
    if (flgs.tsv == _DEFAULT_TSV):
        tsv = flgs.tsv + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        tsv = flgs.tsv
    if (flgs.out == _DEFAULT_OUT):
        out = flgs.out + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        out = flgs.out
    if (flgs.results == _DEFAULT_RESULTS):
        results = flgs.results + _DEFAULT_FREQLIM_SUFFIXDIR
    else:
        results = flgs.results
    if (ChkDir(results)):
        results = ''
        print('Warning: the data root directory will be tried instead, \
for creating results folder')
    if (ChkDir(kreport)):
        kreport = ''
        print('Warning: the data root directory will be tried instead, \
for creating kreports folder')
    if (ChkDir(tsv)):
        tsv = ''
        print('Warning: the data root directory will be tried instead, \
for creating tsv folder')
    if (ChkDir(out)):
        out = ''
        print('Warning: the data root directory will be tried instead, \
for creating outputs (tab) folder')
    try:
        length = int(length)
        quality = int(quality)
        mhl = int(mhl)
    except ValueError:
        print("\nLength, mhl and quality must be integers\n")
        sys.exit()
    if algae == "yes":
        taxon = "cyanobacteria,Bacillariophyta,Bolidophyceae,Charophyce\
ae,Chlorarachniophyceae,Chlorokybophyceae,Chlorophyta,Chrysophyceae,Col\
eochaetophyceae,Cryptophyceae,Dictyochophyceae,Dinophyceae,Euglenida,Eu\
stigmatophyceae,Glaucocystophyceae,Haptophyta,Klebsormidiophyceae,Mesos\
tigmatophyceae,Aurearenophyceae,Xanthophyceae,Chrysoparadoxa,Phaeothamn\
iophyceae,Raphidophyceae,Bangiophyceae,Compsopogonophyceae,Rhodellophyc\
eae,Stylonematophyceae,Synchromophyceae,Synurophyceae,Xanthophyceae,Zyg\
nemophyceae"
    else:
        taxon = "unclassified"

    ## Execute functions
    print("")
    print("\nIt is needed a folder with Centrifuge output files (tab), \
another folder with Centrifuge tsv files, and another folder with krepo\
rt files. Files from same sample must have same file names (except exte\
nsion). For example:\n\
\n\
tab\n\
|\n\
|--> Sample1.tab\n\
|--> Sample2.tab\n\
\n\
tsv\n\
|\n\
|--> Sample1.tsv\n\
|--> Sample2.tsv\n\
\n\
kreports\n\
|\n\
|--> Sample1.txt\n\
|--> Sample2.txt\n\
\n\
- Example for output (tab) files (first three lines):\n\
\n\
readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumM\
atches\n\
SRR2962414.3\tunclassified\t0\t0\t0\t0\t202\t1\n\
SRR2962414.8\tunclassified\t0\t0\t0\t0\t202\t1\n\
\n\
- Example for tsv files (first three lines):\n\
\n\
name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n\
root\t1\tno rank\t0\t57935\t57935\t0.0\n\
Bacteria\t2\tsuperkingdom\t0\t20867\t20867\t0.0\n\
\n\
- Example for kreport files (first three lines):\n\
\n\
 81.14\t9111494\t9111494\tU\t0\tunclassified\n\
 18.86\t2117903\t57951\t-\t1\troot\n\
 17.95\t2015867\t0\t-\t131567\t  cellular organisms\n\
\n")
    
    centrifuge_to_cmplxcruncher(kreport=kreport, tsv=tsv, out=out,
        results=results, length=length, mhl=mhl, quality=quality,
        taxa=taxa, taxon=taxon)

if __name__ == "__main__":
    main()
