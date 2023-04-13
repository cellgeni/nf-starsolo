# nf-starsolo
Our STARsolo repo but implemented in Nextflow

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
* `bin/starsolo_10x_auto.sh` - a script that runs STARsolo outputting: exonic, exonic+intronic and velocyto results.
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Pipeline Arguments:
* `--SAMPLEFILE` - The path to the sample file provided to the pipeline which contains one sample ID per line. This sample is assumed to have CRAM files stored on IRODS.
* `--outdir` - The path to where the results will be saved.
* `--reference` - Tells pipeline which genome to use for alignment (by default GRCh38 2020A is used). This argument is hardcoded and needs to be changed to your local path to the reference index file. 
* `--keep_bams` - Tells the pipeline whether to generate BAM files (default false means do not generate).
