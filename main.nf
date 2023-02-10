#!/usr/bin/env nextflow

nextflow.enable.dsl=1

def helpMessage() {
    log.info"""
    =================
    starsolo pipeline
    =================
    This pipeline runs STARsolo.
    The only parameter you need to input is:
      --SAMPLEFILE /full/path/to/sample/file
    This file should contain a single sampleID per line. 
    An example can be seen here: https://github.com/cellgeni/nf-starsolo/blob/main/examples/example.txt
    The default reference genome is: GRCh38 2020A
    To change these defaults input:
      --reference /path/to/reference    
    """.stripIndent()
}

if (params.HELP) {
  helpMessage()
  exit 0
}

def errorMessage() {
    log.info"""
    ==============
    starsolo error
    ==============
    You failed to provide the --SAMPLEFILE input parameter
    Please provide this parameter as follows:
      --SAMPLEFILE /full/path/to/sample/file
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

process email_startup {

  shell:
  '''
  contents=`cat !{params.SAMPLEFILE}`
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Launched pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk

  Hi there, you've launched Cellular Genetics Informatics' STARsolo pipeline.
  Your parameters are:
  Samplefile: !{params.SAMPLEFILE}
  The Genome GTF used is: !{params.reference}

  Your sample file looks like:
  $contents

  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}
//Puts samplefile into a channel unless it is null, if it is null then it displays error message and exits with status 1.
ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()

//Each line of the sample file is read and then emitted to its own channel which is used as input to the first process, so each sample will be ran in parallel
ch_sample_list
  .flatMap{ it.readLines() }
  .set { ch_samplelines_sf }

process get_starsolo {

  input:
  val(sample) from ch_samplelines_sf

  output:
  //outputs multiple objects: sample name (generated in the process shell so needs env) and a list of all files generated in processing
  set val(sample), path('*.cram') into ch_cram

  shell:
  '''
  myrods.sh -s !{sample} -q off
  '''
}


process crams_to_fastqs {

  input:
  set val(sample), path(cram) from ch_cram

  output:
  set val(sample), env(fastq_dir) into ch_sample_starsolo
  
  shell:
  '''
  fastq_dir=/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/starsolo-results/fastqs
  mkdir -p $fastq_dir
  for cr in !{cram}; do
    !{baseDir}/bin/cram2fastq_10x.sh ${cr}
  done
  for fq in *.fastq.gz; do
    mv $fq "${fastq_dir}/!{sample}_${fq}"
  done
  '''
}

process run_starsolo {

  publishDir "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/${params.sangerID}/${params.timestamp}/starsolo-results", mode: 'copy'

  input:
  set val(sample), val(fastq_dir) from ch_sample_starsolo

  output:
  path(sample) into ch_collect

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{baseDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "true"
  else
    !{baseDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "false"
  fi
  !{baseDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

ch_collect
  .collect()
  .set{ ch_finish }

process email_finish {

  input:
  val(samples) from ch_finish
  
  shell:
  '''
  rm -rf "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/starsolo-results/fastqs"  
  sendmail "!{params.sangerID}@sanger.ac.uk" <<EOF
  Subject: Finished pipeline
  From: noreply-cellgeni-pipeline@sanger.ac.uk

  Hi there, your run of Cellular Genetics Informatics' STARsolo pipeline is complete.

  Results are available here: "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/starsolo-results"

  The results will be deleted in a week so please copy your data to a sensible location, i.e.:
  cp -r "/lustre/scratch126/cellgen/cellgeni/tickets/nextflow-tower-results/!{params.sangerID}/!{params.timestamp}/starsolo-results" /path/to/sensible/location

  The STARsolo command run for each sample can be found inside "starsolo-results/sampleID/cmd.txt"

  The STARsolo script used is availabe here: https://github.com/cellgeni/nf-starsolo/blob/main/bin/starsolo_10x_auto.sh

  Thanks,
  Cellular Genetics Informatics
  EOF
  '''
}
