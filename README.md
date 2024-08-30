# BSA_Project
This project performs bioinformatics analyses for Bulk Segregant Analysis (BSA) to identify genetic variants linked to specific phenotypes.

The pipeline includes quality control, trimming, alignment, variant calling, and filtering for homozygous SNPs. It aims to identify homozygous SNPs associated with resistance to Phytophthora infestans in Solanum species.

## Directory Structure
solanum_americanum/
│
├── ref_genome/
│   └── GCA_03040_.fna
│
├── fastqc_dir/
│
├── trim_dir/
│
├── align_dir/
│
├── bam_dir/
│
├── index_dir/
│
├── vcf_dir/
│
└── homosnps_dir/


## Prerequisites
- Install Snakemake.
- Ensure the following tools are installed:
  - FastQC
  - Trimmomatic
  - BWA
  - SAMtools
  - bcftools
