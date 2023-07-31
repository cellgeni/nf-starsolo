# nf-starsolo
Our [STARsolo repo](https://github.com/cellgeni/STARsolo) but implemented in Nextflow

There are two branches:

`main` - this branch contains the script for running STARsolo on the FARM using Nextflow command line

`nextflow-tower` - this branch conrains the script for running STARsolo on the FARM using Nextflow Tower

## Contents of Repo:
* `main.nf` - the Nextflow pipeline that executes STARsolo.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `examples/samples.txt` - samplefile containing sampleIDs.
* `examples/RESUME-starsolo` - an example run script that executes the pipeline it has 2 hardcoded arguments: `/path/to/sample/file` and `/path/to/config/file` that need to be changed based on your local set up.
* `bin/solo_QC.sh` - a quick qc script that enables quick sanity checks that STARsolo worked correctly.
* `bin/cram2fastq_10x.sh` - a script that converts CRAM files (usual storage type at Sanger) to fastq files. This pipeline expects CRAM files for each sample stored internally on Sanger irods that are accessed using a custom `myrods.sh` script. TO DO - update this to not use myrods.
* `bin/starsolo_10x_auto.sh` - a script that runs STARsolo on 10x data.
* `bin/starsolo_dropseq.sh` - a script that runs STARsolo on dropseq data.
* `bin/starsolo_indrops.sh` - a script that runs STARsolo on indrops data.
* `bin/starsolo_ss2.sh` - a script that runs STARsolo on SmartSeq2 data.
* `bin/starsolo_strt.sh` - a script that runs STARsolo on STRTseq data.
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Pipeline Arguments:
* `-entry` - The entrypoint to specify which determines which sequencing chemistry will be aligned with STARsolo. 
* `--samplefile` - The path to the sample file provided to the pipeline which contains one sample ID per line. This sample is assumed to have CRAM files stored on IRODS.
* `--outdir` - The path to where the results will be saved.
* `--reference` - Tells pipeline which genome to use for alignment (by default GRCh38 2020A is used). This default argument is hardcoded and needs to be changed to your local path to the reference index file. 
* `--keep_bams` - Tells the pipeline whether to generate BAM files (default false means do not generate).
* `--local` - Path to local directory containing the fastqs (default null means look on irods).
* `--index_format` - The index format of the fastq files (defailt i*i means both index1 and index2 file).
