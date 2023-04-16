# seekdeep_illumina_snakemake
a basic workflow for running Nick Hathaway's seekdeep on nanopore

## Installation:
Install conda with:
https://github.com/conda-forge/miniforge#mambaforge

Install snakemake in an environment called snakemake with:
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

## Usage:
 - Change directory to a folder where you want to run the analysis
 - Download the seekdeep_nanopore_general.smk file into this folder
 - Download the sif file from here into the same folder: https://seekdeep.brown.edu/programs/elucidator.sif
 - Download the seekdeep_nanopore_general.yaml file into the same folder
 - Edit the config.yaml file using the instructions in the comments. Use a text editor that outputs unix line endings (e.g. vscode, notepad++, gedit, micro, emacs, vim, vi, etc.)
 - Activate snakemake with:
```bash
conda activate snakemake
```
 - Run snakemake with:
```bash
snakemake -s seekdeep_illumina_general.smk --cores [your_desired_core_count]
```

## Help:
You can read Nick Hathaway's manual here:
https://seekdeep.brown.edu/

If you're in the folder where you downloaded the elucidator.sif file,
you can get help on any seekdeep command with:
```bash
singularity exec elucidator.sif SeekDeep [cmd] -h
```

### three main commands in the snakefile.
  - The first command gets info about the genome (genTargetInfoFromGenomes).
  - The second command sets up an analysis run (setupTarAmpAnalysis).
  - The third command runs 3 seekdeep programs (runAnalysis.sh, no help files).

Here are some example help commands to learn more about these commands:
  - singularity exec elucidator.sif SeekDeep -h
  - singularity exec elucidator.sif SeekDeep genTargetInfoFromGenomes -h
  - singularity exec elucidator.sif SeekDeep setupTarAmpAnalysis -h


### three sub-steps of running seekdeep.
Each of these steps can be tweaked for sensitivity and specificity (via extra_
[step]_cmds at the bottom of the yaml file):
  - The first command extracts amplicon reads (extractor)
  - The second command clusters together similar reads (qluster)
  - The third command processes clusters into haplotypes (processClusters)

Here are some example help commands to learn more about these programs:
  - singularity exec elucidator.sif SeekDeep extractor -h
  - singularity exec elucidator.sif SeekDeep qluster -h
  - singularity exec elucidator.sif SeekDeep processClusters -h
