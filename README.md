# Parallelized RNA-seq analysis pipeline

### Installation
1. Download zip and extract to root of your working directory OR Clone repo to your working diectory.
2. `chomd +x bash_script.sh`

### Use
1. Copy/move all your fastq files to root of your working directory.
2. If running stand alone, install the following together with their dependencies: fastqc, multiqc, flexbar(v3.x), hisat2 (v2.x.x), samtools, htseq and gnu-parallel. If running in HPC, make sure the modules are available.
3. To generate raw transcript count tables, run bash_script.sh OR submit batch to SLURM.
Run standalone as `./bash_script.sh <link to reference genome> <link to anotation> <reference identifier>`
If runnning in HPC, update the arguments in batch file and `sbatch batch`. Compressed GZ file links are expected for reference genome and annotation.
4. Use R_script.rmd to do differential expression and enrichment analysis. All necssary packages are available in CRAN, Github, or Bioconductor and will be installed by the script if necessary.

### Detail

| Activity | Software |
| --- | --- |
| Alignment | HISAT2 |
| Abundance | HTSeq |
| Differential expression | DESeq2 |
| Enrichement | Clusterprofiler |

### Citation
**Sharma A**, Ference CM, Shantharaj D, Baldwin EA, Manthey JA, Jones JB. 2022. Transcriptomic analysis of changes in *Citrus × microcarpa* gene expression post *Xanthomonas citri* subsp. *citri* infection. European Journal of Plant Pathology 162:163–181. doi: http://dx.doi.org/10.1007/s10658-021-02394-6.
