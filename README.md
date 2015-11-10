# fastq2biom
**PE Illumina FASTQ to biom superscript alpha v0.4**

A draft Z shell script automating steps processing Illumina paired-end FASTQ to a BIOM format.

Developed in Z-shell (zsh)

## To use:

* Place all and only the paired-end Illumina files of the this project in this local folder.
* You may need to modify the filenames or the script (in function otu) which parses and joins the PE files based on the naming scheme of your Illumina FASTQ files.

## Dependencies:

The following commands need to be in the $PATH variable, and named as such:"
* mothur - from www.mothur.org"
* uc2otutab.py - from www.drive5.com"
* biom - from biom-format.org"
* usearch8 - from www.drive5.com"

This script requires Python 2 due to the use of uc2otutab.py from drive5


