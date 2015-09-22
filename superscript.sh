#!/usr/bin/env zsh
echo "PE Illumina FASTQ to biom superscript alpha v0.2.1"
echo "--------------------------------------------------------"
echo "Dependencies:"
echo "Currently using environment path in **zshrc**"
echo "mothur" from www.mothur.org"
echo "uc2otutab.py from www.drive5.com"
echo "biom-format from biom-format.org"
echo "usearch8 from www.drive5.com"
echo "--------------------------------------------------------" 
echo "This script Python 2 due to the use of uc2otutab.py from drive5"
echo "Place all and only the paired-end Illumina files of the this project in this local folder."
echo "--------------------------------------------------------" 
echo "Enter the project name that will be given to the output files as an identifier:"
read project
echo "Project name: $project"

echo "Merging reads..."
for file in *R1_001.fastq
do
    if [ -f "$file" ]; then
        sampleID=$(echo "$file" | cut -d "_" -f1) # "_" as delimiter
        usearch8 -fastq_mergepairs $file -reverse ${file%R1_001.fastq}R2_001.fastq -fastqout "${sampleID}_merged.fastq"
    fi
done &> merge_log.txt

echo "Convert FASTQ to FASTA with 1.0 max expected error..."
for file in *merged.fastq; do
    sampleID=$(echo "$file" | cut -d "_" -f1)
    usearch8 -fastq_filter "${sampleID}_merged.fastq" -fastaout "${sampleID}_filtered.fasta" -fastq_maxee 1.0
done &> filter_log.txt

echo "FASTA header renaming..."
for filename in $(ls *_filtered.fasta); do 
samplename=$(echo $filename | cut -d '_' -f1)
sed "-es/^>\(.*\)/>\1;barcodelabel=$samplename;/" < "$filename" > ${samplename}_processed.fasta
done

for filename in *processed.fasta; do
cat $filename >> "${project}_sequence.fasta"
done

# Hardcoded mothur trimming parameters

echo "Sequence trimming, dereplication, singleton removal, and clustering"
mothur "#trim.seqs(fasta="${project}_sequence.fasta", minlength=200, maxlength=500, maxhomop=6)"
usearch8 -derep_fulllength "${project}_sequence.trim.fasta" -fastaout "${project}_seq_derep.fasta" -sizeout
usearch8 -sortbysize "${project}_seq_derep.fasta" -fastaout "${project}_seq_derep2.fasta" -minsize 2
usearch8 -cluster_otus "${project}_seq_derep2.fasta"  -otus "${project}_rep_set.fasta"

echo "Header stripping..."
fasta_number.py "${project}_rep_set.fasta" "OTU_" > "${project}_rep_numbered.fasta"

#------- Abundance calculation by using unfiltered but trimmed reads ------

echo "Generating unfiltered reads for abundance information..."

for file in *merged.fastq; do
    sampleID=$(echo "$file" | cut -d "_" -f1)
    usearch8 -fastq_filter "${sampleID}_merged.fastq" -fastaout "${sampleID}_unfiltered.fasta" #-fastq_maxee 1.0 <- note the exclusion
done &> unfiltered_log.txt

for filename in $(ls *_unfiltered.fasta); do 
samplename=$(echo $filename | cut -d '_' -f1)
sed "-es/^>\(.*\)/>\1;barcodelabel=$samplename;/" < "$filename" > ${samplename}_unprocessed.fasta # <- continue with the un-theme
done
for filename in $(ls *unprocessed.fasta); do
cat $filename >> "${project}_abundance.fasta" # <- new name!
done
mothur "#trim.seqs(fasta="${project}_abundance.fasta", minlength=200, maxlength=500)"

echo "Mapping the abundance data..."
usearch8 -usearch_global "${project}_abundance.trim.fasta" -db "${project}_rep_numbered.fasta" -id 0.97 -strand plus -uc "${project}_readmap.uc"

#-------- End of abundance calculation section ----------------------------

echo "Activating python2 environment for .uc --> txt script" 
uc2otutab.py "${project}_readmap.uc" > "${project}_readmap.txt"

echo "Output HDF5 biom OTU table"
biom convert -i "${project}_readmap.txt" -o "${project}.biom" --table-type="OTU table" --to-hdf5

