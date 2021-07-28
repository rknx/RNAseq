#!/bin/sh

########################################################################

Parallelized RNA-seq transcript abundance pipeline
Author: Anuj Sharma
Email: rknx@outlook.com
Github: https://github.com/rknx/RNAseq
Modified: July 27th, 2021

########################################################################

#!/bin/sh

## Prepare environment ---------------------------------------------------------

# Timestamp function
__ts() { echo -e [`date "+%m/%d %H:%M:%S"`]: $1; }
export -f __ts

# Load required modules to check if they are present.
# If working on standalone PC, install programs and set alias
__ts "Loading modules."
module load parallel hisat2/2.2.0 fastqc multiqc trim_galore samtools htseq

# Workspace
__ts "Creating wokspace sub-folders."
mkdir refs fastq align htseq tmp

# Set tmp directory
export TMPDIR=`pwd`/tmp/

# Move fastq files to right place
mv *.fastq.gz fastq/

## Prepare reference files -----------------------------------------------------

__ts "Preparing reference files."

# Links for reference files
# Reference genome fasta (compressed gz)
export refgenome="$1"

# Annotation gtf for reference genome (compressed gz)
export refannotation="$2"

# Name for reference genome
export refname="$3"

# Prokaryotic (=prok) or eukarytic (=euk)
export refdomain="$4"

# Nunber of CPU for parallel processing
export ncpu="$5"
__ts "$ncpu CPUS specified for processing."

# Reference genome
[[ ! -s refs/genome.fna ]] && \
    __ts "Downloading and extracting reference genome." && ( \
        [[ `wget --spider $refgenome 2> /dev/null` -eq 0 ]] && \
            wget -O - $refgenome | gzip -d > refs/genome.fna || \
            __ts "Reference genome link is not valid." \
    ) || \
    __ts "Reference genome already present.\n\t\tRun 'rm -r refs/*' \
        before this script to rebuild reference files."

# Annotation file
[[ ! -s refs/genes.gtf ]] && \
    __ts "Downloading reference annotation/" && ( \
        [[ `wget --spider $refannotation 2> /dev/null` -eq 0 ]] && \
            wget -O - $refannotation | gzip -d | grep -v "unknown_transcript" |\
                awk 'FS=OFS="\t" {if ($4>$5) {s=$4; $4=$5; $5=s}; print }' > \
                    refs/genes.gtf || \
            __ts "Reference annotation link is not valid." \
    ) || __ts "Reference annotation already present."

# Prepare annotation table
[[ -s refs/genes.gtf && ! -s refs/annotation.table ]] && \
    grep 'gbkey "Gene"' refs/genes.gtf | \
        cut -f2,4 -d '"' | \
        sed 's/"GeneID:/\t/g' | \
        sort > refs/annotation.table


## Indexing --------------------------------------------------------------------

# Load required modules
module -q reset; module load hisat2/2.2.0

__ts "Indexing the reference genome."

# Function to index prokaryotic genome
__index_prokar() {

    [[ ! -s refs/"$refname".2.ht2 ]] && \
        __ts "Prokaryotic organism specified. No splicing expected." && \
        __ts "Starting hisat2 indexing." && \
        hisat2-build -p $ncpu \
            refs/genome.fna \
            refs/"$refname" || \
        __ts "Reference genome is already indexed. Skipping indexing."

}

# Function to index eukaryotic genome
__index_eukar() {

    # Build list of splice sites. This function is built-in hisat2
    [[ ! -s refs/splicesites.tsv ]] && \
        __ts "Building table of splice sites" && \
        extract_splice_sites.py refs/genes.gtf > refs/splicesites.tsv

    # Build list of exons. This function is built-in hisat2
    [[ ! -s refs/exons.tsv ]] && \
        __ts "Building table of exons" && \
        extract_exons.py refs/genes.gtf > refs/exons.tsv

    # Check for necessary files
    [[ ! -s refs/splicesites.tsv && ! -s refs/exons.tsv ]] && \
        __ts "Some reference files are not present." && exit 1

    # Index eukaryotic genome
    [[ ! -s refs/"$refname".2.ht2 ]] && \
        __ts "Eukaryotic organism specified." && \
        __ts "Starting hisat2 indexing." && \
        hisat2-build -p $ncpu \
            --ss refs/splicesites.tsv \
            --exon refs/exons.tsv \
            refs/genome.fna \
            refs/"$refname" || \
        __ts "Reference genome is already indexed. Skipping indexing."
}

# Run the indexing function
[[ $refdomain = "prok" ]] && __index_prokar || __index_eukar

## Quality control -------------------------------------------------------------

# Load required modules
module -q reset; module load parallel fastqc multiqc

# Get list of all fastq files. Extract condition_sample#_lane#
# from condition_sample#_lane#_readpair#_filepart#.fastq.gz
fqlist0=( $(
    find fastq/ -name *.fastq.gz | \
    sed 's/^fastq\///' | \
    rev | cut -d "_" -f1 --complement | rev 
) )
__ts "Checking quality for ${#fqlist0[@]} fastq files." 

# Function to analyze FASTQ seqeunce quality.
__check_qual() {

    __ts "Checking $1."

    # Check quality of individual files
    [[ ! -s fastq/"$1"_001_fastqc.html ]] && \
        fastqc -o fastq/ fastq/"$1"_001.fastq.gz || \
        __ts "QC result already present. Skipping."
}

# Make the function global
export -f __check_qual

# Parallelize the alignment step with GNU parallel
parallel --jobs $ncpu __check_qual ::: ${fqlist0[@]}

# Combine all QC into one file
[[ ! -s fastq/multiqc_report.html ]] && \
    multiqc -o fastq/ fastq/ || \
    __ts "MultiQC report already exists.\n\t\t \
        To rerun, run 'rm fastq/multiqc_report.html' before this script."

# View multiQC output manually or use lynx browser
# The trimming step below should be done manually if necessary.

## Trimming adapters (if necessary) with trim_galore ---------------------------

# Load required modules
module -q reset; module load parallel trim_galore

# Get list of all fastq files. Extract condition_sample#_lane#
fqlist=( $(
    find fastq/ -name *.fastq.gz | \
    sed 's/^fastq\///' | \
    rev | cut -d "_" -f1,2 --complement | rev | \
    sort -u
) )
__ts "Trimming FASTQ files. ${#fqlist[@]} files found."
__ts "After trimming FASTQ files, consider re-running fastqc manually."

# Trimming function
__trim_adapt() {

    __ts "Trimming $1."

    # Trimming step
    [[ ! -s fastq/"$1"_R1_001.fastq.gz.og && \
        ! -s fastq/"$1"_R2_001.fastq.gz.og ]] && \
        __ts "Trimming $1." && \
        trim_galore -j 8 --paired \
            -o fastq/ \
            fastq/"$1"_R1_001.fastq.gz \
            fastq/"$1"_R2_001.fastq.gz || \
        __ts "$1 seems to be trimmed already. Skipping."

    # Swap trimmed and untrimmed files
    [[ -s fastq/"$1"_R1_001_val_1.fq.gz && \
        -s fastq/"$1"_R2_001_val_2.fq.gz ]] && \
        __ts "Backing up original FASTQ files for $1." && \
        mv fastq/"$1"_R1_001.fastq.gz fastq/"$1"_R1_001.fastq.gz.og && \
        mv fastq/"$1"_R2_001.fastq.gz fastq/"$1"_R2_001.fastq.gz.og && \
        find fastq/ -name "$1"*fastqc* -exec \
            sh -c 'x="{}"; mv "$x" "${x}.og"' \; && \
        __ts "Original fastq files are backup up with .og extenion." && \
        __ts "Renaming trimmed files for $1." && \
        mv fastq/"$1"_R1_001_val_1.fq.gz fastq/"$1"_R1_001.fastq.gz && \
        mv fastq/"$1"_R2_001_val_2.fq.gz fastq/"$1"_R2_001.fastq.gz && \
        mv fastq/*_trimming_report.txt fastq/trim_report/

}

# Make the function global
export -f __trim_adapt

# Parallelize the alignment step with GNU parallel
parallel --jobs $(($ncpu/8)) __trim_adapt ::: ${fqlist[@]}

## Alignment with hisat2 -------------------------------------------------------

# Load required modules
module -q reset; module load parallel samtools hisat2/2.2.0

# Get list of all fastq files. Extract condition_sample#_lane#
# from condition_sample#_lane#_readpair#_filepart#.fastq.gz
fqlist2=( $(
    find fastq/ -name *.fastq.gz | \
    sed 's/^fastq\///' | \
    rev | cut -d "_" -f1,2 --complement | rev | \
    sort -u
) )
__ts "Aligning reads to reference genome. ${#fqlist2[@]} files found."

# Alignment function followed by binary conversion.
__align_reads() {

    __ts "Aligning $1."

    shrt=$( echo $1 | sed 's/_L00.//' ) 

    # Alignemnt step
    [[ ! -s align/"$shrt"_sort.bam ]] && \
        hisat2 -p 8 --dta \
            -x refs/"$refname" \
            -1 fastq/"$1"_R1_001.fastq.gz \
            -2 fastq/"$1"_R2_001.fastq.gz \
            -S align/$1.sam || \
        __ts "Sorted bam file already exists for $1. Skipping alignment."
    
    # Convert to binary
    __ts "Converting $1 to binary format."
    [[ -s align/$1.sam ]] && \
        samtools view -@ 8 -bSh align/$1.sam > align/$1.bam || \
        __ts "$1.sam not found."

    # Remove large SAM file is BAM is present.
    __ts "Removing $1.sam."
    [[ -s align/$1.bam ]] && \
        rm align/$1.sam || \
        __ts "$1.bam not found."

}

# Make the function global
export -f __align_reads

# Parallelize the alignment step with GNU parallel
parallel --jobs $(($ncpu/8)) __align_reads ::: ${fqlist2[@]}

# Raise warning if some SAM files still present
[[ `ls align/*.sam 2> /dev/null | wc -l` -gt 0 ]] && \
    __ts "Some SAM files were not converted into BAM files."

## Merge the bam files ---------------------------------------------------------

# Load required modules
module -q reset; module load parallel samtools

# Get list of all bam files.
# Process condition_sample#_lane#.bam to condition only
smlist=( $(
    find align/ -name *_L*.bam | \
    sed 's/^align\///' | \
    rev | cut -d "_" -f1 --complement | rev | \
    sort -u
) )
__ts "Merging, sorting and indexing BAM files. ${#smlist[@]} samples found." 

# Fuction to merge multilane BAM files, following by sorting and indexing.
__merge_sort() {

    __ts "Precessing $1."

    # If multiple lanes present, merge the BAM files.
    [[ `find align/ -size +0 -name "$1*_L*.bam" | wc -l` -gt 1 ]] && \
        __ts "Multiple lanes found. Merging." && \
        samtools cat -@ 4 -o align/$1.bam align/$1*_L*.bam && \
        ( [[ -s align/$1.bam ]] && rm align/$1*_L*.bam )
    
    # If single lane present, rename the BAM file.
    [[ `find align/ -size +0 -name "$1*_L*.bam" | wc -l` -eq 1 ]] && \
        __ts "Single lane found. Renaming." && \
        mv align/$1*_L*.bam align/$1.bam

    # Sort the BAM file.
    [[ -s align/$1.bam ]] && \
        __ts "Sorting binary file for $1." && \
        samtools sort -@ 4 -o align/"$1"_sort.bam align/$1.bam && \
        __ts "Writing alignment metrics for $1." && \
        samtools flagstat -@ 4 align/"$1"_sort.bam > align/"$1"_summary.txt && \
        __ts "Indexing sorted bam for $1." && \
        samtools index -@ 4 -b align/"$1"_sort.bam || \
        __ts "Combined BAM not found. Sorting skipped."
    
    # Remove unsorted BAM file.
    [[ -s align/"$1"_sort.bam.bai ]] && \
        rm align/"$1".BAM
}

# Make the function global
export -f __merge_sort

# Paralellize merging and sorting opeartion
parallel --jobs $(($ncpu/4)) __merge_sort ::: ${smlist[@]}

# Raise warning if some lane specific BAM files remain.
[[ `ls align/*_L*.bam 2> /dev/null | wc -l` -gt 0 ]] && \
    __ts "Some BAM files mere not combined."

## Transcript abundance by HTseq2 ----------------------------------------------

# Load required modules
module -q reset; module load parallel htseq

# Get list of all merged bam files
smlist2=( $(
    find align/ -name *_sort.bam | \
    sed 's/^align\///' | \
    rev | cut -d "_" -f1 --complement | rev | \
    sort -u
) )
__ts "Calculating transcript abundance for ${#smlist2[@]} samples."

# Function to get transcript abundance with htseq
__count_genes() {

    # Run HTseq
    [[ ! -s htseq/$1.meta.tsv ]] && \
        __ts "Counting transcripts for $1." && \
        htseq-count -f bam -s no -r pos -a 1 -t $2 -i gene_id \
            -m intersection-strict \
            align/"$1"_sort.bam \
            refs/genes.gtf > htseq/$1.counts.tsv || \
        __ts "Count abundance is already quantified for $1. Skipping."

    # Separate metadata from counts
    [[ ! -s htseq/$1.meta.tsv ]] && \
        grep '__' htseq/$1.counts.tsv > htseq/$1.meta.tsv && \
        sed -i '/^__/d' htseq/$1.counts.tsv && \
        awk 'BEGIN{OFS="\t"} {c+=$2} END {print "__aligned_count", c}' \
            htseq/$1.counts.tsv >> htseq/$1.meta.tsv && \
        cat htseq/$1.meta.tsv | grep -v "__total" | \
            awk -vOFMT='%.2f' \
                -vtotl=$(grep -v "__total" htseq/$1.meta.tsv | \
                    awk '{x+=$2} END {print x}') \
                'OFS="\t", ORS="%\n" {print $1, $2, $2/totl*100 } \
            END {print "__total_count", totl, "100"}' > \
            htseq/$1.meta.tsv

}

# Make function global
export -f __count_genes

# Parallelize transcript abundance
[[ $refdomain = "prok" ]] && httype="gene" || httype="exon"
parallel --jobs $(($ncpu/8)) __count_genes ::: ${smlist2[@]} ::: $httype

__ts "Job completed."

############################## Alternative module ##############################

# ## Indexing with bowtie2 -----------------------------------------------------

# # Load required modules
# module -q reset; module load bowtie2

# # Index reference genome
# [[ ! -s refs/"$refname".2.bt2 ]] && \
#     __ts "Start bowtie2 indexing." && \
#     bowtie2-build --threads $ncpu \
#         refs/genome.fna \
#         refs/"$refname" || \
#     __ts "Reference genome is already indexed. Skipping indexing."

############################## Alternative module ##############################

# ## Trimming adapters (if necessary) with flexbar --------------------------

# # Load required modules
# module -q reset; module load parallel flexbar

# # Get list of all fastq files. Extract condition_sample#_lane#
# fqlist=( $(
#     find fastq/ -name *.fastq.gz | \
#     sed 's/^fastq\///' | \
#     rev | cut -d "_" -f1,2 --complement | rev | \
#     sort -u
# ) )
# __ts "Trimming FASTQ files. ${#fqlist[@]} files found."
# __ts "After trimming FASTQ files, consider rerunning fastqc manually."

# # Get illumina adapter list
# [[ ! -s fastq/adapters.fa ]] && \
#     __ts "Downloading Illumina adapters." && \
#     wget -O - http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa > \
#         fastq/adapters.fa || \
#     __ts "Adapter sequence is already present."

# # Trimming function
# __trim_adapt() {
#     __ts "Trimming $1."

#     # Trimming step
#     [[ ! -s "$1"_R1_001.fastq.gz.og && ! -s "$1"_R2_001.fastq.gz.og ]] && \
#         flexbar -ao 7 -n 8 -u 4 -m 50 -z GZ \
#             -a fastq/adapters.fa \
#             -r "$1"_R1_001.fastq.gz \
#             -p "$1"_R2_001.fastq.gz \
#             -t $1 || \
#         __ts "$1 has already been trimmed. Skipping."
    
#     # Swap trimmed and untrimmed files
#     [[ -s "$1"_R1_001_val_1.fq.gz && -s "$1"_R2_001_val_2.fq.gz ]] && \
#         __ts "Backing up original FASTQ files for $1." && \
#         mv "$1"_R1_001.fastq.gz "$1"_R1_001.fastq.gz.og && \
#         mv "$1"_R2_001.fastq.gz "$1"_R2_001.fastq.gz.og && \
#         find fastq/ -name "$1"*fastqc* -exec \
#             sh -c 'x="{}"; mv "$x" "${x}.og"' \; && \
#         __ts "Original fastq files are backup up with .og extenion." && \
#         __ts "Renaming trimmed files for $1." && \
#         mv "$1"_1.fastq.gz "$1"_R1_001.fastq.gz && \
#         mv "$1"_2.fastq.gz "$1"_R2_001.fastq.gz
# }

# # Make the function global
# export -f __trim_adapt

# # Parallelize the alignment step with GNU parallel
# parallel --jobs $(($ncpu/8)) __trim_adapt ::: ${fqlist[@]}

############################## Alternative module ##############################

# ## Checking for strandness (optional) ----------------------------------------

# # Load required modules
# module load rseqc hisat2/2.2.0

# # Make output directory
# mkdir strdir

# # Download some tools
# wget -O \
#     - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred > \
#     strdir/gtfToGenePred
# chmod a+x strdir/gtfToGenePred

# wget -O \
#     - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed > \
#     strdir/genePredToBed
# chmod a+x strdir/genePredToBed

# # Make bed file from annotation
# strdir/gtfToGenePred refs/genes.gtf strdir/genes.genePred
# strdir/genePredToBed strdir/genes.genePred strdir/genes.bed12

# # Create small subset of any fastq files
# fqlist=( $(
#     find fastq/ -name *.fastq.gz | \
#     sed 's/^fastq\///' | \
#     rev | cut -d "_" -f1,2 --complement | rev | \
#     sort -u
# ) )
# zcat fastq/"${fqlist[1]}"_R1_001.fastq.gz | \
#     head -n 800000 > strdir/test_R1.fq
# zcat fastq/"${fqlist[1]}"_R2_001.fastq.gz | \
#     head -n 800000 > strdir/test_R2.fq

# # Align the sample fastq files
# hisat2 -p 8 -q \
#     -x refs/$refname \
#     -1 strdir/test_R1.fq \
#     -2 strdir/test_R2.fq \
#     -S strdir/test.sam

# # Infer type of strandness (From RseQC)
# infer_experiment.py -r strdir/genes.bed12 -i strdir/test.sam
# sed -n 3p strdir/result.txt

############################## Alternative module ##############################

# ## Alignment with TopHat2 ----------------------------------------------------

# # Load required modules
# module -q reset; module load parallel tophat

# # Get list of all fastq files. Extract condition_sample#_lane#
# # from condition_sample#_lane#_readpair#_filepart#.fastq.gz
# fqlist2=( $(
#     find fastq/ -name *.fastq.gz | \
#     sed 's/^fastq\///' | \
#     rev | cut -d "_" -f1,2 --complement | rev | \
#     sort -u
# ) )
# __ts "Aligning reads to reference genome. ${#fqlist2[@]} files found."

# # Alignemnt function
# __align_reads() {

#     __ts "Aligning $1."

#     # Alignment step
#     [[ ! -s align/"$1"_sort.bam ]] && \
#         tophat -p 8 \
#             -G refs/genes.gtf \
#             -o align/$1 \
#             refs/"$refname" \
#             fastq/"$1"_R1_001.fastq.gz \
#             fastq/"$1"_R2_001.fastq.gz || \
#         --ts "Sorted bam file already exists for $1. Skipping alignment."
    
#     # Moving the bam file
#     [[ ! -s align/$1/accepted_hits.bam ]] && \
#         mv align/$1/accepted_hits.bam align/$1.bam

# }

# # Make the function global
# export -f __align_reads

# # Parallelize the alignment step with GNU parallel
# parallel --jobs $(($ncpu/8)) __align_reads ::: ${fqlist2[@]}
