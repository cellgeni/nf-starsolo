#!/usr/bin/env nextflow

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

def errorMessage() {
    log.info"""
    ==============
    starsolo error
    ==============
    You failed to provide the SAMPLEFILE or sangerID input parameters
    Please provide these parameters as follows:
      --SAMPLEFILE /full/path/to/sample/file --sangerID user99
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

process get_starsolo {

  input:
  val(sample)

  output:
  //outputs multiple objects: sample name (generated in the process shell so needs env) and a list of all files generated in processing
  tuple val(sample), path('*.cram'), emit: sample_crams  

  shell:
  '''
  myrods.sh -s !{sample} -q off
  '''
}

process crams_to_fastqs {

  input:
  tuple val(sample), path(cram)

  output:
  tuple val(sample), env(fastq_dir), emit: sample_fastqdir 
  
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
  tuple val(sample), val(fastq_dir)

  output:
  path(sample)

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

process email_finish {

  input:
  val(samples)

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

workflow {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    if (params.sangerID == null) {
      errorMessage()
    }
    else {
      email_startup()
      ch_sample_list | flatMap{ it.readLines() } | get_starsolo 
      crams_to_fastqs(get_starsolo.out.sample_crams)
      run_starsolo(crams_to_fastqs.out.sample_fastqdir) | collect | email_finish
    }
  }
}
