#!/usr/bin/env python3

"""
Script to download fastq files from ENA repository

download_fastq_ENA.py

Script purpose: the script downloads fastq files from ENA repository,
given a project accesion, and obtains quality files (FastQC and MultiQC)
"""

__version__ = '1.0'
_verdata = 'Feb 2020'

## Import libraries
from argparse import ArgumentParser
import os
import subprocess

## Main script
def main():
        ## Setting Arguments
        parser = ArgumentParser()

        ## Parameter taxonomy
        parser.add_argument(
                '--project', '-p',
                action ='store',
                required =True,
                help='Project Primary Study Accession.'
        )

        ## Process arguments
        args = parser.parse_args()
        study = args.project

        ## Download fastqs (paired end)
        url = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accessi\
on=' + study + '&result=read_run&fields=study_accession,sample_accessio\
n,secondary_sample_accession,experiment_accession,run_accession,tax_id,\
scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,\
submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_i\
ndex_galaxy&format=tsv&download=true')
        command0 = ("""awk 'FS="\\t", OFS="\\t" { print }' ./""" + study
+  """.txt | cut -f1,9 | awk -F '\\t' 'NR > 1, OFS="\\n" {print $2}'""")
        command1 = 'wget -O ' + study + '.txt "' + url + '"'
        command2 = 'mkdir ./fastq/'
        os.system(command1)
        os.system(command2)
        end = subprocess.check_output(command0, shell=True)
        end = str(end).strip()
        if "PAIRED" in end:
                command3 = ("""awk 'FS="\\t", OFS="\\t" { print }' ./"""
+ study + """.txt | cut -f5,15 | awk -F '\\t' 'OFS="\\n" {print $2}' | \
awk NF | cut -f4,5,6 --delimiter='/' | awk -F '/' 'OFS="\\t" {print $1,\
 $2, $3}' | awk 'NR > 1, OFS="\\n" {print "wget ftp://ftp.sra.ebi.ac.uk\
/vol1/fastq/" "" $1 "/" $2 "/" $3 "/" $3 "_1.fastq.gz;\\n" "wget ftp://\
ftp.sra.ebi.ac.uk/vol1/fastq/" "" $1 "/" $2 "/" $3 "/" $3 "_2.fastq.gz;\
\\n" "mv " $3 "_1.fastq.gz fastq/" $3 "_1.fastq.gz;\\n" "mv " $3 "_2.fa\
stq.gz fastq/" $3 "_2.fastq.gz;\\n" "gzip -d fastq/" $3 "_1.fastq.gz;\\\
n" "gzip -d fastq/" $3 "_2.fastq.gz;\\n"'} > download_""" + study + """\
.sh""")
        else:
                command3 = ("""awk 'FS="\\t", OFS="\\t" { print }' ./\
""" + study + """.txt | cut -f5,15 | awk -F "\\t" 'OFS="\\n" {print $2}\
' | awk NF | cut -f4,5,6 --delimiter="/" | awk -F "/" 'OFS="\\t" {print\
 $1, $2, $3}' | awk 'NR > 1, OFS="\\n" {print "wget ftp://ftp.sra.ebi.a\
c.uk/vol1/fastq/" "" $1 "/" $2 "/" $3 "/" $3 ".fastq.gz;\\n" "mv " $3 "\
.fastq.gz fastq/" $3 ".fastq.gz;\\n" "gzip -d fastq/" $3 ".fastq.gz;\\n\
"'} > download_""" + study + """.sh""")
        os.system(command3)

        ## Generate FastQCs
        command4 = 'mkdir ./fastqc_res/'
        if end == "PAIRED":
                command5 = ("""awk 'FS="\\t", OFS="\\t" { print }' ./\
""" + study + """.txt | cut -f5,15 | awk -F "\\t" 'OFS="\\n" {print $2}\
' | awk NF | cut -f4,5,6 --delimiter="/" | awk -F "/" 'OFS="\\t" {print\
 $1, $2, $3}' | awk 'NR > 1, OFS="\\n" {print "fastqc -o fastqc_res/ fa\
stq/" $3 "_1.fastq;\\n" "fastqc -o fastqc_res/ fastq/" $3 "_2.fastq;\\n\
"}' >> download_""" + study + """.sh""")
        else:
                command5 = ("""awk 'FS="\\t", OFS="\\t" { print }' ./"""
+ study + """.txt | cut -f5,15 | awk -F "\\t" 'OFS="\\n" {print $2}' | \
awk NF | cut -f4,5,6 --delimiter="/" | awk -F "/" 'OFS="\\t" {print $1,\
 $2, $3}' | awk 'NR > 1, OFS="\\n" {print "fastqc -o fastqc_res/ fastq/\
" $3 ".fastq;\\n"}' >> download_""" + study + """.sh""")
        command6 = 'sh ./download_' + study + '.sh'
        command7 = 'rm ./download_' + study + '.sh'
        command8 = 'rm ./' + study + ".txt"
        command9 = 'multiqc fastqc_res/'
        command10 = 'mv multiqc_report.html fastqc_res/multiqc_report.h\
tml'
        command11 = 'rm -r multiqc_data/'
        os.system(command4)
        os.system(command5)
        os.system(command6)
        os.system(command7)
        os.system(command8)
        os.system(command9)
        os.system(command10)
        os.system(command11)

if __name__ == '__main__':
        main()
