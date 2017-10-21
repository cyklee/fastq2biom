#!/usr/bin/env zsh

# Prototyping ITS processing script using:
# cutadapt
# USEARCH v9
# QIIME

echo "Enter project codename that will be given to the output files as an idenfitier:"
read project
echo $(date +%Y-%m-%d\ %H:%M) "Project name: $project"

echo $(date +%Y-%m-%d\ %H:%M) "Merging reads with 5 maximum mismatches or 5% maximum mismatch in the overlapped region."
for file in *R1_001.fastq; do
	 if [ -f "$file" ]; then
		sampleID=$(echo "$file" | cut -d "_" -f1) # "_" as delimiter
		usearch9 -fastq_mergepairs $file -fastqout "${sampleID}_merged.fastq" -log "${sampleID}_merge.log"
	 fi
done

# 1. Filter FASTQ for representative sequence picking and
# 2. convert FASTQ to unfiltered FASTA to map back onto the rep seq 
for i in *_merged.fastq; do
	sampleID=$(echo "$i" | cut -d "." -f1)
	usearch9 -fastq_filter $i -fastq_maxee 1.0 -fastqout "${sampleID}_filtered.fastq" -sample $sampleID -log ${sampleID}_filtered.log
	usearch9 -fastq_filter $i -fastaout "${sampleID}_unfiltered.fasta" -sample $sampleID -log ${sampleID}_unfiltered.log
done

# Concatenate the filtered reads into a single FASTQ file
for file in *_filtered.fastq; do
	cat $file >> "${project}_sequence.fastq"
done

# Concatenate the unfiltered reads into a single FASTA file
for file in *_unfiltered.fasta; do
	cat $file >> "${project}_unfiltered_seq.fasta"
done

# Dereplicate
usearch9 -fastx_uniques "${project}_sequence.fastq" -sizeout -fastaout "${project}_uniques.fasta" -log uniuqes.log

# Some filtering
mothur "#trim.seqs(fasta="${project}_uniques.fasta", minlength=200, maxlength=500, maxhomop=20)"

# OTU clustering
usearch9 -cluster_otus "${project}_uniques.trim.fasta" -minsize 2 -otus "${project}_rep_numbered.fasta" -relabel OTU

# Generating readmap and biom file
usearch9 -usearch_global "${project}_unfiltered_seq.fasta" -db "${project}_rep_numbered.fasta" -id 0.97 -strand plus -uc "${project}_readmap.uc" -biomout "${project}.biom"

# USEARCH SINTAX taxonomy classification
# The output can be integrated into analysis using RDPUtils R package
usearch9 -sintax "${project}_rep_numbered.fasta" -db ~/Bioinformatics/RDP_ITS/rdp_its_v2.fa -strand both -tabbedout "${project}_rdp_sintax.txt" -log "${project}_rdp_sintax.log"

