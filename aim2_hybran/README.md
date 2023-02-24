# Validating New Genes Called by Hybran With RNA-seq Data

## Setup

Final analysis files are kept in: 
```
/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/results/
```

Intermediate output files are kept in: 
```
/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_intermed/
```

Conda environment is:
```
/home/ebishop/.conda/envs/tpm-calc
```
_A new environment was needed because there were conflicts when installing some of these tools within 'emma-base'_

Key tools installed in this environment:
- trimmomatic
- gsnap
- samtools
- AGAT
- bamtools
- tpmcalculator

Input files:
- Raw paired-end TruSeq reads as gzipped FASTQ files
    - `/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_raw/ERR3039930_1.fastq.gz`
    - `/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_raw/ERR3039930_2.fastq.gz`
- gzipped H37Rv reference fasta
    - `/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_raw/H37Rv.fasta.gz`
- Hybran GFF3 with the FASTA at the end removed
    - `/home/ebishop/thesis_projects/aim2_hybran/rnaseq/rnaseq_gitlab/data_raw/H37Rv-NCBI_nofasta.gff`

## Steps
1. `./get_tpm.sh`
2. `python tpm_analysis.py`
