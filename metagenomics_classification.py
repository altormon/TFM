#!/usr/bin/env python3

"""
Script to classificate reads with Centrifuge

metagenomics_classification.py

Script purpose: the script classificates reads with Centrifuge. It also
repairs reads that became disordered or had some mates eliminated, based
on BBMap GitHub project
"""

__version__ = '1.0'
_verdata = 'Feb 2020'


_ROUTE = 'filtered_fastq'
_DEFAULT_FREQLIM_SUFFIXDIR = ''

## Import moduls
from argparse import ArgumentParser
import os
import subprocess

## Create functions
def ChkDir(path):
    error = False
    try:
        os.makedirs(path)
    except(OSError):
        if not os.access(path, (os.W_OK | os.X_OK)):
            print('\nWARNING! Unable to write results in directory "'
                + path + '"...')
            error = True
    return(error)

## Main Program
def main():
        ## Argument Parser Configuration
        parser = ArgumentParser()
        parser.add_argument(
        '-f', '--fastq',
        action='store',
        default=_ROUTE,
        help=('directory path with fastq files (if not present, \''
                + _ROUTE + '\' will be tried)')
        )
        parser.add_argument(
        '-r', '--repair',
        action='store',
        choices=["yes", "no"],
        default="yes",
        help=('Repair paired-end samples (default: yes)')
        )
        parser.add_argument(
        '-i', '--index',
        action='store',
        default="nt",
        help=('Name of Centrifuge index (default: nt)')
        )
        parser.add_argument(
        '-c', '--create',
        action='store',
        choices=["yes", "no"],
        default="no",
        help=('Create Centrifuge index (default: no)')
        )

        ## Parse arguments
        flgs = parser.parse_args()
        repair = flgs.repair
        index = flgs.index
        create = flgs.create
        if (flgs.fastq.replace("./", "").replace("/", "")
                != _ROUTE):
                route = (flgs.fastq.replace("./", "").replace("/", "")
                        + _DEFAULT_FREQLIM_SUFFIXDIR + '/')
        else:
                route = _ROUTE + "/"

        ## Create Centrifuge index
        # if a filtered NCBI-nt is needed, see
        # https://github.com/GW-HIVE/HIVE-lab/tree/master/Filtered_nt
        # it is necessary to use this script in a supercomputer
        if create == "yes":
                c1 = "mkdir nt"
                c2 = "cd nt"
                c3 = "wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz"
                c4 = "gunzip nt.gz"
                c5 = "mkdir taxonomy"
                c6 = "cd taxonomy"
                c7 = "wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxd\
ump.tar.gz"
                c8 = "tar xvzf taxdump.tar.gz"
                c9 = "cd .."
                c10 = "wget ”ftp://ftp.ncbi.nih.gov/pub/taxonomy/access\
ion2taxid/nucl_*.accession2taxid.gz”"
                c11 = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accessi\
on2taxid/pdb.accession2taxid.gz"
                c12 = "gunzip -c *.accession2taxid.gz | awk -v OFS='\t'\
 '{print $2, $3}' >> acc2tax.map"
                c13 = "mv nt.fa nt_unmasked.fa"
                c14 = "dustmasker -infmt fasta -in nt_unmasked.fa -leve\
l 20 -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g' > nt.fa"
                c15 = "centrifuge-build --ftabchars=14 -p 32 --bmax 134\
2177280 --conversion-table acc2tax.map --taxonomy-tree taxonomy/nodes.d\
mp --name-table taxonomy/names.dmp nt.fa nt"
                os.system(c1)
                os.system(c2)
                os.system(c3)
                os.system(c4)
                os.system(c5)
                os.system(c6)
                os.system(c7)
                os.system(c8)
                os.system(c9)
                os.system(c10)
                os.system(c11)
                os.system(c12)
                os.system(c13)
                os.system(c14)
                os.system(c15)

        ## Repair (paired ends)
        if repair == "yes":
                c16 = ("git clone https://github.com/BioInfoTools/BBMap.\
git")
                c17 = ("mv BBMap/sh/rep")
                c18 = ("mv BBMap/sh/repair.sh  ./")
                c19 = ("rm -r BBMap/")
                os.system(c16)
                os.system(c17)
                os.system(c18)
                #os.system(c19)
                for fastq in os.listdir(route):
                        if "_1_" in fastq:
                                c20 = ("repair.sh in=" + route + fastq +
                                        " in2=" + route +
                                        fastq.replace("_1_", "_2_") +
                                        " out=" + route +
                                        fastq.split("_1_")[0] +
                                        "_1_fixed.fastq" + " out2=" +
                                        route + fastq.split("_1_")[0] +
                                        "_2_fixed.fastq" +
                                        " outsingle=" + route +
                                        "singletons_" +
                                        fastq.split("_1_")[0] +
                                        ".fastq")
                                os.system(c20)
        
        ## Classificate (paired ends)
        for fastq in os.listdir(route):
                if "_1_fixed" in fastq:
                        c21 = ("centrifuge -q -p 32 -x " + index +
                                " -1 " + route + fastq + " -2 " +
                                route + fastq.replace("_1_", "_2_") +
                                " -S " + fastq.split("_1_")[0] + ".tab"+
                                " --out-fmt tab -k 1 --report-file " +
                                fastq.split("_1_")[0] + ".tsv")
                        c22 = ("centrifuge-kreport -x " + index + " " +
                                fastq.split("_1_")[0] + ".tab > " +
                                "kreport_" + fastq.split("_1_")[0] +
                                ".txt")
                        os.system(c21)
                        os.system(c22)

if __name__ == '__main__':
        main()
