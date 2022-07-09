<html>
  <head>
    <meta name="google-site-verification" content="FC8xEq44yQepn6GHyTtgwE8snz8H4lnnfgcJE16CIHY" />
  </head>
</html>

# SHARE-seq-alignmentV2
Pipeline for demultiplexing and aligning both ATAC and RNA data generated in SHARE-seqV2.\
The V2 pipeline is much faster than previous version and less demanding. \
Note: The fastq and output files are reformated and different from V1 pipeline.\
**Author: Sai Ma. masai.zju@gmail.com**

# Important note
The ATAC and RNA barcodes would in the format of R1.xxx,R2.xxx,R3.xxx,P1.yy, where xxx is in the range of 01-192, yy is in the range of 01-96.\
The P1.xx indicates the primers that were used to amplify each ATAC and RNA sub-libraries. It is expected to be different between RNA and ATAC assay. The rest part of barcode (R1.xxx,R2.xxx,R3.xxx) would be the same for ATAC and RNA library.\
Depending on the downstream analysis pipeline, sometimes the comma in the barcode would be converted to period (R1.xxx,R2.xxx,R3.xxx,P1.yy --> R1.xxx.R2.xxx.R3.xxx.P1.yy).\
A barcode translation table indicating the P1.xx used in the assay would be beneficial, when submitting data to GEO.

# Installation
This pipeline requires following packages to be properly installed and added to system path: GNU parallel, Bcl2fastq, fastp, zcat, STAR, bowtie2, python3 (all python code upadted to python3), umi_tools, samtools, picard (2.14.1, newer version may result in error), R, featureCounts, read_distribution.py from RSeQC, bedtools. 

The SHARE-seq-alignment scripts can be directly downloaded from the github website.\
[https://github.com/masai1116/SHARE-seq-alignmentV2/](https://github.com/masai1116/SHARE-seq-alignmentV2/)

After downloading all scripts, update the general configuration section in main script "Share_seqV2_example.sh":
1) myPATH # where the SHARE-seq scripts are installed. e.g. myPATH='/mnt/users/Script/share-seq-github-v2/'
2) pythohPATH # where python3 is installed e.g. pythohPATH='/usr/bin/python/' 
3) picardPATH # where picard is installed e.g. picardPATH='/mnt/bin/picard/picard.jar'

The pipeline also requres gtf files and aligner index files to be download and unziped into the right location.\
GTF files can be downloaded [here](https://drive.google.com/file/d/1j4MOLRz4033EEyj2Oo_CG0FcZQDoBYwQ/view?usp=sharing).\
Bowtie2 index files (Hg19 and mm10) can be downloaded [here](https://drive.google.com/file/d/1bXIxznwirsZ6DZhqK1gw6ZKlj-UjFRhn/view?usp=sharing).\
Assuming SHARE-seq aligment scripts are installed to "/home/SHARE-seq-alignment/", the gtf files should be placed in the "/home/SHARE-seq-alignment/gtf/" folder.\
Reference genome files can be downloaded [here](https://drive.google.com/file/d/1PW5FD9GfhEOl5kQgAFxZSF57kcXtlrJG/view?usp=sharing).\
Assuming SHARE-seq aligment scripts are installed to "/home/SHARE-seq-alignment/", the reference genome files should be placed in the "/home/SHARE-seq-alignment/refGenome/" folder.\
If user prefer to use their own genome build, the bowtie2 index files should be placed in the "/home/SHARE-seq-alignment/refGenome/bowtie2/" folder.\
Four sets of index files (hg38, hg19, mm10 and hg19-mm10 combined genome) for star aligner should be prepared according to star aligner [manual](https://github.com/alexdobin/STAR), or downloded from here: [hg19](https://drive.google.com/file/d/1IXI4DP-mjh2qc-EQe1WnWJQOCVEI4KVX/view?usp=sharing), [mm10](https://drive.google.com/file/d/1n0UwzOeUbX7TIBOrcBbXjgH3i0UH-Ka5/view?usp=sharing), [combined genome](https://drive.google.com/file/d/15Z2YMUDiavYG0s9zLFAbbwqA0VhVNu-f/view?usp=sharing).\
The unziped index files should be placed in the "/home/SHARE-seq-alignment/refGenome/star/hg38", "/home/SHARE-seq-alignment/refGenome/star/hg19", "/home/SHARE-seq-alignment/refGenome/star/mm10", and "/home/SHARE-seq-alignment/refGenome/star/both", respectively.\
The index file for hg19-mm10 combined genome can be downloaded from [10x Genomics website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).

# How to run the script?
A small set of fastq files for testing are in the test_fastq_nova/ folder
Before running, three sections in the main script "Split_seq_example.sh" need to be updated for each run, inlcuding 
A) paths B) sample configuration C) fastq configuration. After update all specific information in Share_seqV2_example.sh and config.example.yaml, run the script by ```./Share_seqV2_example.sh```
## A) paths
1) rawdir=./example_fastq/ # where the raw data is
2) dir=./test/ # where output data will be stored; ./example_output/ shows the output from the example fastqs 
3) yaml=./config_example.ymal # where the ymal configuration file is

## B) sample configuration
1) Project=(BMMC.RNA BMMC.ATAC) # use differnt name for each sample 
2) Type=(RNA ATAC)  # ATAC or RNA
3) Genomes=(hg38 hg38) # both mm10 hg19 hg38 \
RawReadsPerBarcode and ReadsPerBarcode options are designed to remove barcodes with too few reads and speed up processing. 
4) ReadsPerBarcode=(10 10) # reads cutoff to barcodes: 100 for full run; 10 for QC run
5) keepMultiMapping=(F F)  # default F; F for species mixing or cell lines, T for low yield tissues (only keep the primarily aligned reads), doesn't matter for ATAC. Allowing multi-mapping redas will increase the percent of mito reads
6) keepmito=(F F) ## default F, remove mito reads for ATAC analysis

## C) fastq configuration
1) Start=Fastq # Bcl or Fastq
2) Runtype=QC # QC or full, QC only analyze 12M reads to get a quick sense of data
3) chem=fwd # rev or fwd: speficy the chemistry used in sequencing, nova1.5 & nextseq use rev; nova1.0 uses fwd

## RNA-seq options
The pipeline also offers flexible RNA-seq specific options for advanced users. 
1) removeSingelReadUMI=F # T or F; default is F. If T, UMIs with single read will be removed.
2) keepIntron=T # T or F; default is T. If F, intronic RNA reads will be discarded.
3) cores=16
4) genename=gene_name # default gene_name; gene_name (official gene symbol) or gene_id (ensemble gene name)
5) mode=fast # fast or regular; default fast; fast: dedup with custom  script; regular: dedup with umitools
fast mode gives more UMIs because taking genome position into account when dedup. It doesn't collapse UMIs map to different position. The lib size estimation is not accurate.

## Other option
cleanup=T # this will remove intermediate fastqs and save disk space. When using the pipeline the first time, set it to F could help troubleshoot.

# Sample barcode table
SHARE-seq allows mutiplexing samples in one run. We use ymal file to store PCR barcode information (P1.xx, refers to the Ad1.xx primers used in the PCR step). See ```config_example.yaml``` as an example. This file needs to be updated for each sample and each sequencing run. When multiple sublibraries are sequenced at the same time, simply add addtional P1.xx to the yaml file. (e.g. P1.13)
```
---
Project1:
  Name: BMMC.RNA
  Primer:
  - P1.12
  - P1.13
  Type: RNA
Project2:
  Name: BMMC.ATAC
  Primer:
  - P1.04
  Type: ATAC
```        
P1.xx can be P1.01, P1.02, ..., P1.96.\
The detialed information about these barcode can be found in [SHARE-seq manuscript](https://www.sciencedirect.com/science/article/pii/S0092867420312538).

# Pipeline output
After successfully running the pipeline, data for each project will be kept in the folder with the same name as the project. 
Several QC plots will be generated.\
For ATAC, the most important files are the fragment file after removing duplicates, mito reads and shifting Tn5 cutting position (.fragments.tsv.gz) and a summarize report (.counts.csv.gz).\
For RNA, the most important files are UMIxCell matrix (.gene.bc.matrices.h5) and a summarize report (.counts.csv.gz).\
This pipeline currently keeps many intermedia files. If preferred, they can be manually removed afterwards.

# Read data example
A set of fastq files for human bone marrow cells experiment can be downloaded [here](https://drive.google.com/drive/folders/1d0gfb7qrBL76MMh0JPRX3z9-1Wqxt0vd?usp=sharing). 
It only takes a few minutes to run the pipeline on the test data.
More SHARE-seqV2 data will be released soon. 


# Cite us
For more details, please refer to [Ma et al. Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin, Cell 2020](https://www.sciencedirect.com/science/article/pii/S0092867420312538)
