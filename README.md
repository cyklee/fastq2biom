```
    ___                       ______  _     _             
   / __)             _       (_____ \| |   (_)            
 _| |__ _____  ___ _| |_ ____  ____) ) |__  _  ___  ____  
(_   __|____ |/___|_   _) _  |/ ____/|  _ \| |/ _ \|    \ 
  | |  / ___ |___ | | || |_| | (_____| |_) ) | |_| | | | |
  |_|  \_____(___/   \__)__  |_______)____/|_|\___/|_|_|_|
                           |_|                            
```
# fastq2biom
**PE Illumina FASTQ to biom superscript Beta v1.1**

A draft Z shell script automating steps processing Illumina paired-end FASTQ to a BIOM format.
Developed in Z-shell (zsh).

## What's new?

Now with optional taxonomic assignment, downstream qiime analysis, and a cool title banner!


## To use:

* Place all and only the paired-end Illumina files of the this project in this local folder.
* You may need to modify the filenames or the script (in function `otu()`) which parses and joins the PE files based on the naming scheme of your Illumina FASTQ files.
* I used Anaconda `conda create -n python2 python=2.7 anaconda` to create a Python 2 envrionment which is activated via `source activate python2` in the script. If you're not natively in Python 3, you can probably comment out the lines (starting with `source`) that activate and deactivates the environment.

## Dependencies:

The following commands need to be in the $PATH variable, and named as such:
* mothur - from http://www.mothur.org
* uc2otutab.py - from http://www.drive5.com
* biom - from http://biom-format.org
* usearch8 - from http://www.drive5.com
* Optional QIIME scripts

This script requires Python 2 due to the use of uc2otutab.py from drive5 and QIIME v1.9.1

## TODO:
* Generate FASTQC report of the input FASTQ files.
* Report sequence size distribution (perhaps via mothur summary.seqs?) <-- R script written
* Perhaps use the size distribution to calculate trim.seq parameters (e.g. adjust to capture a portion of the distribution, within reason).
* Report pre-processing read counts vs. post-processing read counts <- simple FASTQ counting script availble perhaps vs. merged read count *2 ?

## Issues
* Log files need to also log files being worked on (e.g. during merging) so they can be linked to the diagnostic output



