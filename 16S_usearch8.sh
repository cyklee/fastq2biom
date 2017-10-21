#!/usr/bin/env zsh
# Prints welcome message and takes project name
welcome() {
echo "PE Illumina FASTQ to biom superscript for USEARCH v8.x"
echo "Please note due to command changes, this script will NOT work with <v8, v9, or v10."
echo "Developed for Z-shell (zsh)"
echo "--------------------------------------------------------"
echo "Dependencies:"
echo "The following commands need to be in the \$PATH variable, and named as such:"
echo "mothur - from www.mothur.org"
echo "uc2otutab.py - from www.drive5.com"
echo "biom - from biom-format.org"
echo "usearch8 - from www.drive5.com"
echo "Optional qiime integration (in Python 2)."
echo "--------------------------------------------------------" 
echo "This script requires Python 2 due to the use of uc2otutab.py from drive5"
echo "Place all and only the paired-end Illumina files of the this project in this local folder."
echo "--------------------------------------------------------" 
echo "Enter the project name that will be given to the output files as an identifier:"
read project
echo "Project name: $project"
read -q "response?Do you want to conduct downstream QIIME analysis? [y/N]"
# instead of interactive Q/A perhaps convert these into switches (-q for qiime and -a for alpha diversity mode (include singleton))
printf "\n"
}
# Picks representative OTU sequences using usearch8 and mothur, they need to be in $PATH
otu() {
echo "Merging reads..."
for file in *R1_001.fastq
do
    if [ -f "$file" ]; then
        sampleID=$(echo "$file" | cut -d "_" -f1) # "_" as delimiter
        usearch8 -fastq_mergepairs $file -reverse ${file%R1_001.fastq}R2_001.fastq -fastqout "${sampleID}_merged.fastq"
    fi
done &> merge_log.txt # The shell report redirect may be replaced with -report switch (introduced in usearch v8.1.1859)

echo "Convert FASTQ to FASTA with 1.0 max expected error..."
for file in *merged.fastq; do
    sampleID=$(echo "$file" | cut -d "_" -f1)
    usearch8 -fastq_filter "${sampleID}_merged.fastq" -fastaout "${sampleID}_filtered.fasta" -fastq_maxee 1.0
done &> filter_log.txt
# No good inbuilt report function in usearch yet, but should list filtered filename in the collated report for traceability.

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
}

#Abundance calculation by using unfiltered but trimmed reads
abundance() {
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
}

cleanup() {
echo "Cleaning up..."
mkdir log
mv unfiltered_log.txt log/
mv merge_log.txt log/
mv filter_log.txt log/

mkdir original_sequences
mv *R1*fastq original_sequences/
mv *R2*fastq original_sequences/

mkdir intermediary_sequences
mv *.fasta intermediary_sequences/
mv *.fastq intermediary_sequences/
cp "intermediary_sequences/${project}_rep_numbered.fasta" ./
}

biom_conversion() {
# Make biom file from classic OTU table
# This is dependent on my anaconda python2 environment set up, rename to your env name
echo "Activating python2 environment for .uc --> txt script" 
source activate python2

echo "Converting uc output to classic OTU table"
uc2otutab.py "${project}_readmap.uc" > "${project}_readmap.txt"

# convert "classic" otu table to biom file
echo "Output HDF5 biom OTU table"
biom convert -i "${project}_readmap.txt" -o "${project}.biom" --table-type="OTU table" --to-hdf5

echo "Deactivate python2 environment, may be activated again for QIIME downstream."
source deactivate
}


# Optional QIIME analysis, these won't work with non-16S data (e.g., ITS, 18S, rbcL)
qiime() {
echo "Activating python2 envrionment for QIIME..."
source activate python2

echo "Assigning taxonomy, this might take a while."
assign_taxonomy.py -i "${project}_rep_numbered.fasta" -m rdp -o rdp_tax

echo "Adding taxonomy to biom file, use this biom file for downstream analysis from here on."
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp "rdp_tax/${project}_rep_numbered_tax_assignments.txt"  -i "${project}.biom" -o "${project}_tax.biom"

echo "Making taxonomy stack bar graphs. Using the catchall 'summarize_taxa_through_plots.py' script"
summarize_taxa_through_plots.py -i "${project}_tax.biom" -o "${project}_tax_summary"

echo "Conducting three steps to produce a phylogenetic tree, for phylogeny-based analysis downstream (Unifrac, PD_whole_tree)."
echo "Aligning representative OTU sequence."
align_seqs.py -i "${project}_rep_numbered.fasta"
echo "Filtering the alignment."
filter_alignment.py -i "pynast_aligned/${project}_rep_numbered_aligned.fasta" -o pynast_aligned/
echo "Making the phylogenetic tree."
make_phylogeny.py -i "pynast_aligned/${project}_rep_numbered_aligned_pfiltered.fasta"

echo "All done! Now exit python2 environment."
source deactivate
}

welcome
otu
abundance
cleanup
biom_conversion

if [[ "$response" =~ ^[Yy]$ ]]
then
qiime
fi


