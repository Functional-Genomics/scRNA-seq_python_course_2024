## Introduction to Cell Ranger for scRNA-seq data analysis using Python

This section is prepared and delivered by [Jiawei Wang](http://jiawei.nohup.cc/) at EMBL-EBI, with the assistance of `Jinzheng Ren` from the Australian National University.

### 0. `Demonstrations/cellranger_setup.md`
This file contains the instructions for setting up **Cell Ranger**.

### 1. `Data/cellranger_data/`
This directory contains the files necessary for demonstration practices. Our focus will primarily be on two folders:

#### 1.1 `references/`
   * `Homo_sapiens.GRCh38.dna.chromosome.21.fa`: A FASTA file for your genome reference.
   * `gencode.v41.primary_assembly.annotation.chr21.gtf`: A GTF file with annotated genes. 
   * `cellranger_mkref.sh`: A script used to build the reference datasets.

#### 1.2 `fastq/`
This folder contains the FASTQ files, which will be used to generate the count matrix.

#### 1.3 `cellranger_count.sh`
This script will be used to generate the count matrix.

### 2. Run Cell Ranger commands

#### 2.1 Building reference
To build the references, navigate to the `cellranger_data/references/` directory from the `Data/` folder and run the `cellranger_mkref.sh` script:

   ```
   cd cellranger_data/references
   sh cellranger_mkref.sh
   ```

#### 2.2 Generating count matrix
After generating the reference, create the count matrix by running the `cellranger_count.sh` script in the `cellranger_data/` directory:
   ```
   cd ..
   sh cellranger_count.sh
   ```

### Acknowledgement
Slides and demonstration materials are primarily reused, with slight adaptation, from [Introduction to single-cell RNA-seq analysis, University of Cambridge](https://github.com/bioinformatics-core-shared-training/UnivCambridge_ScRnaSeqIntro_Base/tree/main).