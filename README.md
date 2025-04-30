# tRNA modifications tune m6A-dependent mRNA decay

This repository contains the code to create the figures for Linder, Sharma, Wu et al. 2025


## Usage

1. Run the main snakemake workflow.
```bash
snakemake --cores 10 --use-conda --printshellcmds --rerun-incomplete
```
This will process the raw data for the miCLIP and SLAM-Seq experiments. 

2. Run the ribosome profiling snakemake workflow located in the folder `snakemake_ribo_di_seq`.
This will process the raw data for the ribosome profiling experiments.

3. Run the quantseq snakemake workflow located in the folder `snakemake_Quantseq`.
This will process the raw data for the QuantSeq experiments.

    Note that in order to run the GRAND-SLAM analysis, you need to install the GEDI and GRAND-SLAM packages (see https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM) into `workflow/resources/bin/GRAND-SLAM` and obtain the corresponding license (see https://www.uni-wuerzburg.de/sft/erfindungen-patente-und-lizenzen-jmu-und-ukw/download-software-for-scientific-purposes/)   

4. Run the miCLIP analysis notebook `analysis/analysis_miCLIP`.
This will call m6A sites from the miCLIP bam files generated in the snakemake workflow.

5. Run the main analysis notebooks `analysis/analysis_human` and `analysis/analysis_mouse`.
This will generate the figures in the manuscript.

6. Run the cancer analysis notebook `analysis/cancer`.
This will generate the figures in the manuscript.
