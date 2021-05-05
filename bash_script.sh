#!/bin/sh
#SBATCH --account <account>
#SBATCH --qos <qos>
#SBATCH --job-name <rnaseq>
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user <email>
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=17
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --output=rnaseq_%j.log

date

## Prepare environment --------------------------------------------------

# Load required modules
# If working on standalone PC, install programs and set alias
module load parallel fastqc multiqc hisat2 samtools htseq #stringtie

# Workspace
mkdir refs fastq align htseq #stringtie

## Prepare reference files -------------------------------------------------

# Links for C. sinensis
refgenome="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/317/415/GCF_000317415.1_Csi_valencia_1.0/GCF_000317415.1_Csi_valencia_1.0_genomic.fna.gz"
refannotation="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/317/415/GCF_000317415.1_Csi_valencia_1.0/GCF_000317415.1_Csi_valencia_1.0_genomic.gtf.gz"
refname="ref_genome" # E.g., Csi_valencia, Sly_Heinz

# Reference genome
wget -O - $refgenome | gzip -d > refs/genome.fna

# Annotation file
wget -O - $refannotation | gzip -d | grep -v "unknown_transcript" > refs/genes.gtf

# Preapre annotation table
grep 'gbkey "Gene"' refs/genes.gtf | cut -f2,4 -d '"' | sed 's/"GeneID:/\t/g' | sort > refs/annotation.table

## Indexing ---------------------------------------------------------------

# Run slurm with --=mem 16g --ncpu 9 --time 2:00:00
# Graph-building may need more memory.

# Build list of splice sites
# This is bilt-in with hisat2
extract_splice_sites.py refs/genes.gtf > refs/splicesites.tsv

# Build list of exons
# This is bilt-in with hisat2
extract_exons.py refs/genes.gtf > refs/exons.tsv

# Index genome
hisat2-build -p 8 --ss refs/splicesites.tsv --exon refs/exons.tsv refs/genome.fna refs/"$refname"

## Checking for strandness (optional) -------------------------------------

ml rseqc
mkdir strandness

# Download some tools
wget -o strandness/gtfToGenePred http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
chmod a+x strandness/gtfToGenePred

wget -o strandness/genePredToBed http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
chmod a+x strandness/genePredToBed

# Make bed file from annotation
strandness/gtfToGenePred refs/genes.gtf strandness/genes.genePred
strandness/genePredToBed strandness/genes.genePred strandness/genes.bed12

# Create small subset of any fastq files
fqpre="C1-1_S7_L001"

zcat fastq/"$fqpre"_R1_001.fastq.gz | head -n 800000 > strandness/test_R1.fq
zcat fastq/"$fqpre"_R2_001.fastq.gz | head -n 800000 > strandness/test_R2.fq

# Align the sample fastq files
hisat2 -p 8 -q -x refs/Csi_valencia -1 strandness/test_R1.fq -2 strandness/test_R2.fq -S strandness/test.sam

# Infer type of strandness
# From RseQC
infer_experiment.py -r strandness/genes.bed12 -i strandness/test.sam
sed -n 3p strandness/result.txt

## Quality control

# Check quality of individual files
fastqc -o fastqc fastq/*.fastq.gz

# Combine all QC into one file
multiqc ./fastq

# View multiQC output manually or use lynx browser

## Trimming adapters (if necessary) ---------------------------------------

ml flexbar
mkdir trim

# Get illumina adapter list
wget -o trim/adapters.fa http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa

# Trim reads using adapters
for i in `ls fastq/*.fastq.gz | sed 's/R._001.fastq.gz//g' | tr ' ' '\n' | sort -u | tr '\n' ' '`
do
    flexbar \
        --adapter-min-overlap 7 \
        --adapter-trim-end RIGHT \
        --adapters trim/adapters.fa \
        --pre-trim-left 13 \
        --max-uncalled 300 \
        --min-read-length 25 \
        --threads 8 \
        --zip-output GZ \
        --reads "$i"R1_001.fastq.gz \
        --reads2 "$i"R2_001.fastq.gz \
        --target trim
done

## Alignment -----------------------------------------------------------------

# Run slurm with --=mem 17g --ncpu 17 --time 5:00:00
# 2*8 CPUs for alignment and 1 extra CPU for the rest.
# As citrus does not have a lot of repeats, 1 GB memory /CPU is fine.

# Get list of all fastq files. Process
# condition_sample#_lane#_readpair#_filepart#.fastq.gz to
# condition_sample#_lane# only
mv *.fastq.gz fastq/
fqlist=`basename -s fastq.gz fastq/*.fastq.gz | cut -d "_" -f1-3 | uniq`

# Function to align with hisat
# Note that reads are not stranded
# Reference genome has already been indexed
# samtools converts output .sam to binary .bam
align() {
    hisat2 -p 8 \
        -x refs/"refname" \
        --dta \
        -1 fastq/"$1"_R1_001.fastq.gz \
        -2 fastq/"$1"_R2_001.fastq.gz \
        -S align/$1.sam
    samtools view -bSh align/$1.sam > align/$1.bam; \
}

# Make the function global
export -f align

# Parallelize the alignment step with GNU parallel
parallel --jobs 2 align ::: ${fqlist[@]}

# Remove sam file if bam file found.
for file in `basename -s .bam align/*.bam`; do rm align/$file.sam; done
[[ `find ./align/ -name *.sam | wc -l` -gt 0 ]] && echo "Some SAM files were not converted into BAM files."

## Merge the bam files --------------------------------------------------

# Run slurm with --=mem 2 Gb --ncpu 5 --time 0:15:00

# Get list of all bam files.
# Process condition_sample#_lane#.bam to condition only
conds=`basename -s .bam align/*.bam | cut -d "-" -f1 | uniq`

# Fuction to merge bams by conditions + sorting
mergebam() {
    samtools merge align/$1.merged.bam align/$1*.bam
    samtools sort -o align/$1.sort.bam align/$1.merged.bam
    samtools index -b align/$1.sort.bam
}

# Mke the function global
export -f mergebam

# Paralellize merge
parallel --jobs 4 mergebam ::: ${conds[@]}

# Delete previous bam files
newbam=`basename -s .merged.bam align/*.merged.bam`
for file in ${newbam[@]}; do rm -r align/$file*_L00*.bam; done

## (Deprecated) Transcript abundance by stringtie -------------------------

# --=mem 9 Gb --ncpu 9 --time 3:00:00

# bamlist=`basename -s .sort.bam align/*.sort.bam`
# expr() {
#     stringtie -p 2 -e -B\
#         -G refs/genes.gtf \
#         -o stringtie/$1.transcripts.gtf \
#         -A stringtie/$1.abundance.tsv \
#         align/$1.sort.bam
# }
# export -f expr
# parallel --jobs 4 expr ::: ${bamlist[@]}

## Transcript abundance by HTseq2 ----------------------------------

# --=mem 5 Gb --ncpu 5 --time 1:00:00

# Get list of all merged bam files
bamlist=`basename -s .sort.bam align/*.sort.bam`

# Function to get transcript abundance with htseq
htseq() {
    htseq-count \
        --format bam \
        --order pos \
        --mode intersection-strict \
        --minaqual 1 \
        --type exon \
        --idattr gene_id \
        align/$1.sort.bam \
        refs/genes.gtf > htseq/$1.counts.tsv
    grep '__' htseq/$1.counts.tsv > htseq/$1.meta.tsv
    sed -i '/^__/d' htseq/$1.counts.tsv
}

# Make function global
export -f htseq

# Parallelize transcript abundance
parallel --jobs 4 htseq ::: ${bamlist[@]}

# Remove metadata from count files
for i in `basename -s .counts.tsv htseq/*.counts.tsv`
do
    grep '__' htseq/$i.counts.tsv > htseq/$i.meta.tsv
    sed -i '/^__/d' htseq/$i.counts.tsv
done
