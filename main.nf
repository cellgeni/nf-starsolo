#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    You failed to provide the SAMPLEFILE input parameter
    Please provide these parameters as follows:
      --SAMPLEFILE /full/path/to/sample/file
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
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
  fastq_dir="!{params.outdir}/fastqs"
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

  publishDir "${params.outdir}", mode: 'copy'

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
  rm -rf "!{params.outdir}/fastqs"
  '''
}

workflow {
  if (params.HELP) {
    helpMessage()
    exit 0
  }
  else {
    ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | get_starsolo 
    crams_to_fastqs(get_starsolo.out.sample_crams)
    run_starsolo(crams_to_fastqs.out.sample_fastqdir)
  }
}
