#!/usr/bin/env zsh
# Prints welcome message and takes project name
welcome() {
echo "PE Illumina FASTQ to biom superscript alpha USEARCH v9 branch"
echo "Developed for Z-shell (zsh)"
echo "--------------------------------------------------------"
echo "Dependencies:"
echo "The following commands need to be in the \$PATH variable, and named as such:"
echo "mothur - from www.mothur.org"
echo "uc2otutab.py - from www.drive5.com"
echo "biom - from biom-format.org"
echo "usearch9 - from www.drive5.com"
echo "Optional qiime integration (in Python 2)."
echo "--------------------------------------------------------" 
echo "This script requires Python 2 due to the use of uc2otutab.py from drive5"
echo "Place all and only the paired-end Illumina files of the this project in this local folder."
echo "--------------------------------------------------------" 
echo "Enter the project name that will be given to the output files as an identifier:"
read project
echo $(date +%Y-%m-%d\ %H:%M) "Project name: $project"
printf "\n"
read -q "relaxed?Do you want to relax the FASTQ merge pair criteria? [y/N]"
printf "\n"
read -q "qiime?Do you want to conduct downstream QIIME analysis? [y/N]"
printf "\n"
read -q "alpha?Do you want to keep singletons for alpha diversity analysis? [y/N]"
printf "\n"
read -q "denoise?Do you want to conduct UNOISE high resolution analysis? [y/N]"
printf "\n"
}
# Picks representative OTU sequences using usearch9 and mothur, they need to be in $PATH
otu() {
echo $(date +%Y-%m-%d\ %H:%M) "Merging reads with 5 maximum mismatches or 5% maximum mismatch in the overlapped region."

if [[ "$relaxed" =~^[Yy]$ ]]; then
	echo "Running relaxed paired-end FASTQ merging  (-fastq_maxdiffs 10 -fastq_maxdiffpct 10 -fastq_trunctail 5)"
	for file in *R1_001.fastq; do
		if [ -f "$file" ]; then
			sampleID=$(echo "$file" | cut -d "_" -f1) # "_" as delimiter
			usearch9 -fastq_mergepairs $file -fastqout "${sampleID}_merged.fastq" -log "${sampleID}_merge.log"
		fi
	done
else
	echo "Running standard paired-end FASTQ merging"
	for file in *R1_001.fastq; do
		if [ -f "$file" ]; then
			sampleID=$(echo "$file" | cut -d "_" -f1) # "_" as delimiter
			usearch9 -fastq_mergepairs $file -fastqout "${sampleID}_merged.fastq" -log "${sampleID}_merge.log"
		fi
	done
fi

echo $(date +%Y-%m-%d\ %H:%M) "Converting FASTQ to FASTA with 1.0 max expected error."

# Attempting to use the >SampleID.readnumber syntax here, SampleID will have to be alphanumeric
# Here @ = $sampleID as done by usearch (using _ and . delimiters)
# An alternative to -relabel is to use `-sample string` and append sample=string to the read label
for file in *merged.fastq; do
	sampleID=$(echo "$file" | cut -d "_" -f1)
	usearch9 -fastq_filter "${sampleID}_merged.fastq" -fastq_maxee 1.0 -fastqout "${sampleID}_filtered.fastq" -sample $sampleID -log ${sampleID}_filter.log
	usearch9 -fastq_filter "${sampleID}_merged.fastq" -fastaout "${sampleID}_unfiltered.fasta" -sample $sampleID -log ${sampleID}_unfilter.log
done

# Concatenate the files in to a single file
for filename in *filtered.fastq; do
	cat $filename >> "${project}_sequence.fastq"
done

for filename in *unfiltered.fasta; do
	cat $filename >> "${project}_unfiltered_seq.fasta"
done
}

unoise() {
echo $(date +%Y-%m-%d\ %H:%M) "EXPERIMENTAL: Conducting UNOISE high resolution analysis (http://www.drive5.com/usearch/manual/unoise_pipeline.html)."
time usearch9 -fastx_uniques "${project}_sequence.fastq" -sizeout -fastqout "${project}_uniques.fastq" -log uniques.log
usearch9 -unoise "${project}_uniques.fastq" -tabbedout unoise.txt -fastaout "${project}_denoised.fasta"
time usearch9  -usearch_global "${project}_unfiltered_seq.fasta" -db "${project}_denoised.fasta" -id 1.00 -strand plus -uc "${project}_denoised.uc" -biomout "${project}_denoised.biom"

# Make copies of the required file for QIIME anlaysis in the anticipated file names
if [["$qiime" =~ ^[Yy]$ ]]; then
cp "${project}_denoised.fasta" "${project}_rep_numbered.fasta"
cp "${project}_denoised.biom" "${project}.biom"
fi
}

normal(){
# Note singleton is not removed here
echo $(date +%Y-%m-%d\ %H:%M) "Dereplicating and sorting the unique reads by abundance."
usearch9 -fastx_uniques "${project}_sequence.fastq" -sizeout -fastaout "${project}_uniques.fasta" -log uniques.log

# Hardcoded mothur filtering parameters. This is mostly an insurance step to remove highly unusual sequences that are unlikely to be biological.
# usearch only does minlength discard, perhaps I will migrate the mothur command to vsearch (with min/maxseqlength)
# https://github.com/jimmyodonnell/banzai has an implementation with grep and awk
	
echo $(date +%Y-%m-%d\ %H:%M) "Sequence trimming, dereplication, singleton removal, and clustering."
mothur "#trim.seqs(fasta="${project}_uniques.fasta", minlength=200, maxlength=500, maxhomop=6)"

if [[ "$alpha" =~ ^[Yy]$ ]]; then
	usearch9 -cluster_otus "${project}_uniques.trim.fasta"  -otus "${project}_rep_numbered.fasta" -relabel OTU
else
# No need to invoke an external python script to relabel the reads? I like this
	usearch9 -cluster_otus "${project}_uniques.trim.fasta" -minsize 2 -otus "${project}_rep_numbered.fasta" -relabel OTU
fi
time usearch9  -usearch_global "${project}_unfiltered_seq.fasta" -db "${project}_rep_numbered.fasta" -id 0.97 -strand plus -uc "${project}_readmap.uc" -biomout "${project}.biom"
}

cleanup() {
echo $(date +%Y-%m-%d\ %H:%M) "Cleaning up..."
mkdir log
mv *.log log/

mkdir original_sequences
mv *R1*fastq original_sequences/
mv *R2*fastq original_sequences/

mkdir intermediary_sequences
mv *.fasta intermediary_sequences/
mv *.fastq intermediary_sequences/

cp "intermediary_sequences/${project}_rep_numbered.fasta" ./
}

# biom_conversion() is unnecessary for usearch9. Keeping it for now for debug purposes.
biom_conversion() {
# Make biom file from classic OTU table
# This is dependent on my anaconda python2 environment set up, rename to your env name
echo "Activating python2 environment for .uc --> txt script" 
source activate python2

echo "Converting uc output to classic OTU table"
uc2otutab.py "${project}_readmap.uc" > "${project}_readmap.txt"
#convert "classic" otu table to biom file
echo "Output HDF5 biom OTU table"
biom convert -i "${project}_readmap.txt" -o "${project}.biom" --table-type="OTU table" --to-hdf5

echo "Deactivate python2 environment, may be activated again for QIIME downstream."
source deactivate
}

qiime_analysis() {
printf "\n"
echo $(date +%Y-%m-%d\ %H:%M) "- Downstream analysis with QIIME."; printf "\n"
echo "Activating python2 envrionment for QIIME."; printf "\n"
source activate python2

echo $(date +%Y-%m-%d\ %H:%M)  "Assigning taxonomy, this might take a while."
time assign_taxonomy.py -i "${project}_rep_numbered.fasta" -m rdp -o rdp_tax

echo $(date +%Y-%m-%d\ %H:%M)  "Adding taxonomy to biom file, use this biom file for downstream analysis from here on."
time biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp "rdp_tax/${project}_rep_numbered_tax_assignments.txt"  -i "${project}.biom" -o "${project}_tax.biom"

echo $(date +%Y-%m-%d\ %H:%M) "Making taxonomy stack bar graphs. Using the catchall 'summarize_taxa_through_plots.py' script"
time summarize_taxa_through_plots.py -i "${project}_tax.biom" -o "${project}_tax_summary"

echo $(date +%Y-%m-%d\ %H:%M) "Conducting three steps to produce a phylogenetic tree, for phylogeny-based analysis downstream (Unifrac, PD_whole_tree)."

echo "Aligning representative OTU sequence."
align_seqs.py -i "${project}_rep_numbered.fasta"
echo "Filtering the alignment."
filter_alignment.py -i "pynast_aligned/${project}_rep_numbered_aligned.fasta" -o pynast_aligned/
echo "Making the phylogenetic tree."
make_phylogeny.py -i "pynast_aligned/${project}_rep_numbered_aligned_pfiltered.fasta"
}

welcome
otu

if [[ "$denoise" =~ ^[Yy]$ ]]; then
	echo "Running UNOISE ZOTU mode"
	unoise
else
	echo "Running 97% OTU clustering mode"
	normal
fi

cleanup

if [[ "$qiime" =~ ^[Yy]$ ]]; then
	qiime_analysis
fi

