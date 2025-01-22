# seekdeep_illumina_snakemake
a basic workflow for running Nick Hathaway's seekdeep on illumina. This version splits up jobs into individual snakemake submissions.

## Installation:
 - Install mamba: https://github.com/conda-forge/miniforge#install (don't forget
to do conda init and follow the instructions to log out and back in at the end)
 - Create a mamba environment and install snakemake there:
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

### Setup your environment:
 - Change directory to a folder where you want to run the analysis
 - clone this repository with git clone web_address - you can get the web_address from the green 'code' button
 - Download the sif file from here into the same folder: https://seekdeep.brown.edu/programs/elucidator.sif

## Usage:
 - Edit the seekdeep_illumina_general.yaml file using the instructions in the
comments. Use a text editor that outputs unix line endings (e.g. vscode,
notepad++, gedit, micro, emacs, vim, vi, etc.)
 - If snakemake is not your active conda environment, activate snakemake with:
```bash
mamba activate snakemake
```
 - Run all steps with (e.g. if you have 8 cores available on your machine):
```bash
snakemake -s setup_run.smk --cores 8
snakemake -s run_extractor.smk --cores 8
snakemake -s finish_process.smk --cores 8
```
 - You can also run all steps with:
```bash
bash run_all_steps.sh
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
