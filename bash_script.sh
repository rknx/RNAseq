#!/bin/sh

########################################################################

Parallelized RNA-seq transcript abundance pipeline
Author: Anuj Sharma
Email: rknx@outlook.com
Github: https://github.com/rknx/rna-seq-deseq2
Modified: May 13th, 2021

########################################################################

## Prepare environment --------------------------------------------------

# Timestamp function
ts() { echo -e [`date +%m/%d` `date +%H:%M:%S`]: $1; }
export -f ts

# Load required modules to check if they are present.
# If working on standalone PC, install programs and set alias
ts "Loading modules."
module load samtools parallel hisat2/2.2.0 htseq fastqc multiqc flexbar

# Workspace
ts "Creating wokspace sub-folders."
mkdir refs fastq align htseq trim untrimmed tmp

# Set tmp directory
export TMPDIR=`pwd`/tmp/

# Move fastq files to right place
mv *.fastq.gz fastq/

## Prepare reference files -------------------------------------------------

ts "Preparing reference files."

# Links for reference files
# Reference genome fasta (compressed gz)
export refgenome="$1"

# Annotation gtf for reference genome (compressed gz)
export refannotation="$2"

# Name for reference genome
export refname="$3"

# Reference genome
[[ ! -f refs/genome.fna ]] && \
    ts "Downloading and extracting reference genome." && ( \
        [[ `wget --spider $refgenome 2> /dev/null` -eq 0 ]] && \
            wget -O - $refgenome | gzip -d > refs/genome.fna || \
            ts "Reference genome link is not valid." \
    ) || ts "Reference genome already present.\nRun 'rm -r refs/*' before this script to rebsuild reference files."

# Annotation file
[[ ! -f refs/genes.gtf ]] && \
    ts "Downloading reference annotation/" && ( \
        [[ `wget --spider $refannotation 2> /dev/null` -eq 0 ]] && \
            wget -O - $refannotation | gzip -d | grep -v "unknown_transcript" > refs/genes.gtf || \
            ts "Reference annotation link is not valid." \
    ) || ts "Reference annotation already present."

# Preapre annotation table
[[ -f refs/genes.gtf && ! -f refs/annotation.table ]] && \
    grep 'gbkey "Gene"' refs/genes.gtf | cut -f2,4 -d '"' | sed 's/"GeneID:/\t/g' | sort > refs/annotation.table


## Indexing ---------------------------------------------------------------

# Load required modules
module load hisat/2.2.0

# Graph-building may need more memory.

ts "Indexing the reference genome."

# Build list of splice sites. This function is built-in hisat2
[[ ! -f refs/splicesites.tsv ]] && \
    extract_splice_sites.py refs/genes.gtf > refs/splicesites.tsv

# Build list of exons. This function is built-in hisat2
[[ ! -f refs/exons.tsv ]] && \
    extract_exons.py refs/genes.gtf > refs/exons.tsv

# Index genome
[[ ! -s refs/"$refname".2.ht2 ]] && ( \
    [[ ! -f refs/genome.fna || ! -f refs/splicesites.tsv || ! -f refs/exons.tsv ]] && \
        ts "Some reference files are not present." && exit 1
    ) && \
    hisat2-build -p 8 \
        --ss refs/splicesites.tsv \
        --exon refs/exons.tsv \
        refs/genome.fna \
        refs/"$refname" \
    || \
    ts "Reference genome is already indexed."

# ## Checking for strandness (optional) -------------------------------------

# Load required modules
module load parallel rseqc hisat2/2.2.0

# mkdir strandness

# # Download some tools
# wget -o strandness/gtfToGenePred http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
# chmod a+x strandness/gtfToGenePred

# wget -o strandness/genePredToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
# chmod a+x strandness/genePredToBed

# # Make bed file from annotation
# strandness/gtfToGenePred refs/genes.gtf strandness/genes.genePred
# strandness/genePredToBed strandness/genes.genePred strandness/genes.bed12

# # Create small subset of any fastq files
# fqlist=( `find fastq/ -name *.fastq.gz | sed 's/^fastq\///' | rev | cut -d "_" -f1,2 --complement | rev | sort -u` )

# zcat fastq/"${fqlist[1]}"_R1_001.fastq.gz | head -n 800000 > strandness/test_R1.fq
# zcat fastq/"${fqlist[1]}"_R2_001.fastq.gz | head -n 800000 > strandness/test_R2.fq

# # Align the sample fastq files
# hisat2 -p 8 -q -x refs/"$refname" -1 strandness/test_R1.fq -2 strandness/test_R2.fq -S strandness/test.sam

# # Infer type of strandness (From RseQC)
# infer_experiment.py -r strandness/genes.bed12 -i strandness/test.sam
# sed -n 3p strandness/result.txt

## Quality control --------------------------------------------------------

# Load required modules
module load parallel fastqc multiqc

# Get list of all fastq files. Extract condition_sample#_lane#
# from condition_sample#_lane#_readpair#_filepart#.fastq.gz
fqlist0=( `find fastq/ -name *.fastq.gz | sed 's/^fastq\///' | rev | cut -d "_" -f1 --complement | rev` )
ts "Checking quality for ${#fqlist0[@]} fastq files." 

# Function to analyze FASTQ seqeunce quality.
qcheck() {

    ts "Checking $1."

    # Check quality of individual files
    [[ ! -f fastq/"$1"_001_fastqc.html ]] && \
        fastqc -o fastq/ fastq/"$1"_001.fastq.gz || \
        ts "QC result already present. Skipping."
}

# Make the function global
export -f qcheck

# Parallelize the alignment step with GNU parallel
parallel --jobs 8 qcheck ::: ${fqlist0[@]}

# Combine all QC into one file
[[ ! -f fastq/multiqc_report.html ]] && \
    multiqc -o fastq/ fastq/ || \
    ts "MultiQC report already exists.\nTo rerun, run 'rm fastq/multiqc_report.html' before this script."

# View multiQC output manually or use lynx browser

## Trimming adapters (if necessary) ---------------------------------------

# Load required modules
module load parallel flexbar

# Get illumina adapter list
[[ ! -f trim/adapters.fa ]] && \
    ts "Downloading Illumina adapters." && \
    wget -O - http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa > trim/adapters.fa || \
    ts "Adapter sequence is already present."

# Get list of all fastq files. Extract condition_sample#_lane#
# from condition_sample#_lane#_readpair#_filepart#.fastq.gz
# rev is used to deal with subfolder with '_'.

fqlist=( `find fastq/ -name *.fastq.gz | sed 's/^fastq\///' | rev | cut -d "_" -f1,2 --complement | rev | sort -u` )
ts "Trimming FASTQ files. ${#fqlist[@]} files found."
ts "After trimming FASTQ files, consider rerunning fastqc again manually to verify quality."
ts "To analyze untrimmed files, run 'rm trim/*; mv untrimmed/* trim/'."

# Trim reads using adapters
trim() {
    ts "Trimming $1."

    flexbar -ao 7 -n 8 -u 4 -m 50 -z GZ \
        -a trim/adapters.fa \
        -r fastq/"$1"_R1_001.fastq.gz \
        -p fastq/"$1"_R2_001.fastq.gz \
        -t trim/$1
    
    # Move original fastq to backup directory and new file to fastq directory
    [[ -s trim/"$1"_1.fastq.gz && -s trim/"$1"_2.fastq.gz ]] && \
        mv fastq/"$1"_* untrimmed/ && \
        mv trim/"$1"_1.fastq.gz fastq/"$1"_R1_001.fastq.gz && \
        mv trim/"$1"_2.fastq.gz fastq/"$1"_R2_001.fastq.gz && \
        ts "$1 successfully moved."
}

# Make the function global
export -f trim

# Parallelize the alignment step with GNU parallel
parallel --jobs 1 trim ::: ${fqlist[@]}

## Alignment -----------------------------------------------------------------

# Load required modules
module load parallel hisat/2.2.0 samtools

# Get list of all fastq files. Extract condition_sample#_lane#
# from condition_sample#_lane#_readpair#_filepart#.fastq.gz
fqlist2=( `find fastq/ -name *.fastq.gz | sed 's/^fastq\///' | rev | cut -d "_" -f1,2 --complement | rev | sort -u` )
ts "Aligning reads to reference genome. ${#fqlist2[@]} files found." 

# Function to align with hisat (not stranded) followed by binary conversion.
align() {

    ts "Aligning $1."

    # Alignemnt step
    hisat2 -p 8 --dta \
        -x refs/"$refname" \
        -1 fastq/"$1"_R1_001.fastq.gz \
        -2 fastq/"$1"_R2_001.fastq.gz \
        -S align/$1.sam
    
    # Convert to binary
    ts "Converting $1 to binary format."
    [[ -f align/$1.sam ]] && \
        samtools view -@ 8 -bSh align/$1.sam > align/$1.bam || \
        ts "$1.sam not found."

    # Remove large SAM file is BAM is present.
    ts "Removing $1.sam."
    [[ -s align/$1.bam ]] && \
        rm align/$1.sam || \
        ts "$1.bam not found."
}

# Make the function global
export -f align

# Parallelize the alignment step with GNU parallel
parallel --jobs 1 align ::: ${fqlist2[@]}

# Raise warning if some SAM files still present
[[ `ls align/*.sam 2> /dev/null | wc -l` -gt 0 ]] && \
    ts "Some SAM files were not converted into BAM files."

## Merge the bam files --------------------------------------------------

# Load required modules
module load parallel samtools

# Get list of all bam files.
# Process condition_sample#_lane#.bam to condition only
smlist=( `find align/ -name *_L*.bam | sed 's/^align\///' | rev | cut -d "_" -f1 --complement | rev | sort -u` )
ts "Merging, sorting and indexing BAM files. ${#smlist[@]} samples found." 

# Fuction to merge multilane BAM files, following by sorting and indexing.
mergebam() {

    ts "Precessing $1."

    # If multiple lanes present, merge the BAM files.
    [[ `ls align/$1*_L*.bam 2> /dev/null | wc -l` -gt 1 ]] && \
        ts "Multiple lanes found. Merging." && \
        samtools cat -@ 4 -o align/$1.bam align/$1*_L*.bam && \
        ( [[ -s align/$1.bam ]] && rm align/$1*_L*.bam )
    
    # If single lane present, rename the BAM file.
    [[ `ls align/$1*_L*.bam 2> /dev/null | wc -l` -eq 1 ]] && \
        ts "Single lane found. Renaming." && \
        mv align/$1*_L*.bam align/$1.bam

    # Sort the BAM file.
    [[ -s align/$1.bam ]] && \
        ts "Sorting $1.bam" && \
        samtools sort -@ 4 -o align/"$1"_sort.bam align/$1.bam || \
        ts "Combined BAM not found. Sorting skipped."
    
    # Remove unsorted BAM file.
    [[ -s align/"$1"_sort.bam ]] && \
        rm align/"$1".bam

    # Index the sorted BAM file.
    [[ -s align/"$1"_sort.bam ]] && \
        ts "Indexing "$1"_sort.bam" && \
        samtools index -@ 4 -b align/"$1"_sort.bam || \
        ts "Sorting operation was incomplete. Indexing step skipped."
}

# Make the function global
export -f mergebam

# Paralellize merging and sorting opeartion
parallel --jobs 2 mergebam ::: ${smlist[@]}

# Raise warning if some lane specific BAM files remain.
[[ `ls align/*_L*.bam 2> /dev/null | wc -l` -gt 0 ]] && \
    ts "Some BAM files mere not combined."

## Transcript abundance by HTseq2 ----------------------------------

# Load required modules
module load parallel htseq

# Get list of all merged bam files
smlist2=( `find align/ -name *_sort.bam | sed 's/^align\///' | rev | cut -d "_" -f1 --complement | rev | sort -u` )
ts "Calculating transcript abundance for ${#smlist2[@]} samples." 

# Function to get transcript abundance with htseq
htseq() {

    ts "Processing $1."

    # Run HTseq
    htseq-count -f bam -r pos -m intersection-strict -a 1 \
        align/"$1"_sort.bam \
        refs/genes.gtf > htseq/$1.counts.tsv

    # Separate metadat from counts
    grep '__' htseq/$1.counts.tsv > htseq/$1.meta.tsv
    sed -i '/^__/d' htseq/$1.counts.tsv
}

# Make function global
export -f htseq

# Parallelize transcript abundance
parallel --jobs 8 htseq ::: ${smlist2[@]}

ts "Job completed."
