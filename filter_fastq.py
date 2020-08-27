#!/usr/bin/env python3

"""
Script to filter fastq files

filter_fastq.py

Script purpose: the script filters fastq files downloaded from ENA
repository, and obtains quality files (FastQC and MultiQC)
"""

__version__ = '1.0'
_verdata = 'Feb 2020'


_DEFAULT_ROUTE = 'fastq'
_NEW_ROUTE = 'filtered_fastq'
_DEFAULT_FREQLIM_SUFFIXDIR = ''

## Import moduls
from argparse import ArgumentParser
import os
import subprocess

## Create functions
def rev_com(seq=None):
    if seq is not None:
        rev_seq = seq[::-1]
        rev_seq = rev_seq.translate(str.maketrans("ATGCUN", "TACGAN"))
        return(rev_seq)

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
        default=_DEFAULT_ROUTE,
        help=('directory path with fastq files (if not present, \''
                + _DEFAULT_ROUTE + '\' will be tried)')
        )
        parser.add_argument(
        '-r', '--results',
        action='store',
        default=_NEW_ROUTE,
        help=('path of the results directory relative to the data\
 directory (if not present, \'' + _NEW_ROUTE + '\' will be tried)')
        )
        parser.add_argument(
        '-a', '--adapter',
        action='store',
        default="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
        help=('Sequence of the adapter for trimming (rev comp is\
 obtained and analyzed automatically)')
        )
        parser.add_argument(
        '-c', '--cores',
        action='store',
        default=1,
        help=('cores used for analysis (default: 1)')
        )
        parser.add_argument(
        '-l', '--length',
        action='store',
        default=101,
        help=('Maximum length of sequences for filtering ' +
                '(minimum length = maximum length - adapter length)')
        )
        parser.add_argument(
        '-q', '--quality',
        action='store',
        default=10,
        help=('Quality for sequence filtering')
        )
        parser.add_argument(
        '-p', '--program',
        action='store',
        choices=["cutadapt", "fastx", "both"],
        default="both",
        help=('Program to use for preprocesing, cutadapt and/or FastX t\
oolkit (default: both)')
        )

        ## Parse arguments
        flgs = parser.parse_args()
        seq = flgs.adapter.upper()
        rev = rev_com(seq)
        cores = flgs.cores
        length = flgs.length
        quality = flgs.quality
        program = flgs.program
        if (flgs.fastq.replace("./", "").replace("/", "")
                != _DEFAULT_ROUTE):
                route = (flgs.fastq.replace("./", "").replace("/", "")
                        + _DEFAULT_FREQLIM_SUFFIXDIR + '/')
        else:
                route = _DEFAULT_ROUTE + "/"
        if (flgs.results.replace("./", "").replace("/", "")
                != _NEW_ROUTE):
                res = (flgs.results.replace("./", "").replace("/", "")
                        + _DEFAULT_FREQLIM_SUFFIXDIR + '/')
        else:
                res = _NEW_ROUTE + "/"
        if (ChkDir(res)):
                res = ''
                print('Warning: the filtered_fastq root directory will \
be tried instead, for writing results')

        ## Filtering fastqs using cutadapt
        paireds = []
        if program == "cutadapt" or program == "both":
                for fastq in os.listdir(route):
                        ## Single end
                        if "_1" not in fastq and "_2" not in fastq:
                                command1 = ("cutadapt -j " + cores +
                                        " -g " + seq + " -o " +
                                        res + fastq.split(".")[0] +
                                        "_trim1." +
                                        fastq.split(".")[-1] +
                                        " " + route + fastq)
                                command2 = ("cutadapt -j " + cores +
                                        " -a " + rev + " -o " +
                                        res + fastq.split(".")[0] +
                                        "_trim2." +
                                        fastq.split(".")[-1] + " " +
                                        res + fastq.split(".")[0] +
                                        "_trim1." +
                                        fastq.split(".")[-1])
                                command3 = ("cutadapt -j " + cores +
                                        " -m " + length + " -M " +
                                        length + " -q " +
                                        quality + " -o " + res +
                                        fastq.split(".")[0] +
                                        "_filtered_cutadapt." +
                                        fastq.split(".")[-1] + " " +
                                        res + fastq.split(".")[0] +
                                        "_trim2." +
                                        fastq.split(".")[-1])
                                os.system(command1)
                                os.system(command2)
                                os.system(command3)
                        ## Paired end
                        else:
                                if fastq.split("_")[0] not in paireds:
                                        paireds.append(
                                                fastq.split("_")[0])
                                        command4 = ("cutadapt -j " +
                                                cores + " -g " + seq +
                                                " -G " + seq +
                                                " -o " + res + 
                                                fastq.split("_")[0] +
                                                "_1_trim1." +
                                                fastq.split(".")[-1] +
                                                " -p " + res +
                                                fastq.split("_")[0] +
                                                "_2_trim1." +
                                                fastq.split(".")[-1] +
                                                " " + route +
                                                fastq.split("_")[0] +
                                                "_1." +
                                                fastq.split(".")[-1] +
                                                " " + route +
                                                fastq.split("_")[0] + 
                                                "_2." +
                                                fastq.split(".")[-1])
                                        command5 = ("cutadapt -j " +
                                                cores + " -a " + rev +
                                                " -A " + rev +
                                                " -o " + res + 
                                                fastq.split("_")[0] +
                                                "_1_trim2." +
                                                fastq.split(".")[-1] +
                                                " -p " + res +
                                                fastq.split("_")[0] +
                                                "_2_trim2." +
                                                fastq.split(".")[-1] +
                                                " " + res +
                                                fastq.split("_")[0] +
                                                "_1_trim1." +
                                                fastq.split(".")[-1] +
                                                " " + res +
                                                fastq.split("_")[0] + 
                                                "_2_trim1." +
                                                fastq.split(".")[-1])
                                        command6 = ("cutadapt -j " +
                                                cores + " -m " +
                                                length + " -M " +
                                                length + " -q " +
                                                quality + " -o " + res +
                                                fastq.split("_")[0] +
                                                "_1_filtered_cutadapt."+
                                                fastq.split(".")[-1] +
                                                " " + res +
                                                fastq.split("_")[0] +
                                                "_1_trim2." +
                                                fastq.split(".")[-1])
                                        command7 = ("cutadapt -j " +
                                                cores + " -m " +
                                                length + " -M " +
                                                length + " -q " +
                                                quality + " -o " + res +
                                                fastq.split("_")[0] +
                                                "_2_filtered_cutadapt."+
                                                fastq.split(".")[-1] +
                                                " " + res +
                                                fastq.split("_")[0] +
                                                "_2_trim2." +
                                                fastq.split(".")[-1])
                                        os.system(command4)
                                        os.system(command5)
                                        os.system(command6)
                                        os.system(command7)

        ## Filtering fastqs using FastX toolkit
        if program == "fastx" or program == "both":
                for fastq in os.listdir(route):
                        command8 = ("fastx_trimmer -i " + route +
                                fastq + " -o " + res +
                                fastq.split(".")[0] + "_trim." +
                                fastq.split(".")[-1])
                        command9 = ("fastx_clipper -a " + seq + " -i " +
                                res + fastq.split(".")[0] + "_trim." +
                                fastq.split(".")[-1] + " -o " + res +
                                fastq.split(".")[0] + "_clip1." +
                                fastq.split(".")[-1])
                        command10 = ("fastx_clipper -a " + rev +
                                " -i " + res +
                                fastq.split(".")[0] + "_clip1." +
                                fastq.split(".")[-1] + " -o " + res +
                                fastq.split(".")[0] + "_clip2." +
                                fastq.split(".")[-1])
                        command11 = ("fastq_quality_trimmer -t " +
                                quality + " -l " + length + " -i " +
                                res + fastq.split(".")[0] + "_clip2." +
                                fastq.split(".")[-1] + " -o " + res +
                                fastq.split(".")[0] + "_qtrim." +
                                fastq.split(".")[-1])
                        command12 = ("fastq_quality_filter -q " +
                                quality + " -i " + res +
                                fastq.split(".")[0] + "_qtrim." +
                                fastq.split(".")[-1] + " -o " + res +
                                fastq.split(".")[0] +
                                "_filtered_fastx." +
                                fastq.split(".")[-1])
                        os.system(command8)
                        os.system(command9)
                        os.system(command10)
                        os.system(command11)
                        os.system(command12)

        ## Save FastQC and MultiQC results of filtered fastqs
        if "filtered_fastqc_fastx" not in os.listdir("./"):
                os.mkdir("filtered_fastqc_fastx")
        if "filtered_fastqc_cutadapt" not in os.listdir("./"):
                os.mkdir("filtered_fastqc_cutadapt")
        for fastq in os.listdir(res):
                if fastq.split("_")[-1] == "fastx.fastq":
                        command13 = ("fastqc -o filtered_fastqc_fastx/ "
                                + res + fastq)
                        os.system(command13)
                elif fastq.split("_")[-1] == "cutadapt.fastq":
                        command14 = ("fastqc -o filtered_fastqc_cutadap\
t/ " + res + fastq)
                        os.system(command14)
        command15 = "multiqc filtered_fastqc_cutadapt/"
        command16 = "mv multiqc_report.html filtered_fastqc_cutadapt/mu\
ltiqc_report.html"
        command17 = "multiqc filtered_fastqc_fastx/"
        command18 = "mv multiqc_report_1.html filtered_fastqc_fastx/mul\
tiqc_report.html"
        os.system(command15)
        os.system(command16)
        os.system(command17)
        os.system(command18)

if __name__ == '__main__':
        main()
