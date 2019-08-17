```
    ___                       ______  _     _             
   / __)             _       (_____ \| |   (_)            
 _| |__ _____  ___ _| |_ ____  ____) ) |__  _  ___  ____  
(_   __|____ |/___|_   _) _  |/ ____/|  _ \| |/ _ \|    \ 
  | |  / ___ |___ | | || |_| | (_____| |_) ) | |_| | | | |
  |_|  \_____(___/   \__)__  |_______)____/|_|\___/|_|_|_|
                           |_|                            
```

* NOTE fastq2biom is no longer being maintained as I have migrated towards ASV-based approach with data2.

# fastq2biom
**PE Illumina FASTQ to biom superscript**

A draft Z shell script automating steps processing Illumina paired-end FASTQ to a BIOM format.
Developed for Z-shell (zsh).

## What's new?

1. ITS (usearch9) support added.
2. 16S (usearch10) support added.
3. 16S relaxed merge mode added.
4. In usearch10, unoise3 support and updated codes added.


## To use:

* Place all and only the paired-end Illumina files of the this project in this local folder.
* The ID of the sample in the file name should be deliminted by underscore; e.g. "SampleID_L001_R1.fastq".
* I used Anaconda `conda create -n python2 python=2.7 anaconda` to create a Python 2 envrionment which is activated via `source activate python2` in the script. If you're not natively in Python 3, you can probably comment out the lines (starting with `source`) that activate and deactivates the environment.

## Dependencies:

The following commands need to be in the $PATH variable, and named as such:
* mothur - from http://www.mothur.org
* biom - from http://biom-format.org
* usearch[n] - from http://www.drive5.com (Note the version-dependent scripts)
* Optional QIIME scripts

## TODO:
* Generate FASTQC report of the input FASTQ files.
* [fastx_learn](http://www.drive5.com/usearch/manual/cmd_fastx_learn.html) error rate determination.
* Report sequence size distribution (gist: cyklee/fastx_sequence_length.sh).
* Report pre-processing read counts (simple grep for unfiltered fastq).
* Testing 16S_usearch10 script for unoise3 and singleton inclusion (alpha=Y).

## Possibly?
* Snakemake?
* Perhaps use the size distribution to calculate trim.seq parameters (e.g. adjust to capture a portion of the distribution).

## Issues
* Previous log file issue solved due to logging being supported by newer usearch.



