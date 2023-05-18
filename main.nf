#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    =================
    starsolo pipeline
    =================
    This pipeline runs STARsolo.
    There are 2 required parameters:
      --samplefile /full/path/to/sample/file
      -entry starsolo-mode
    The samplefile should contain a single sampleID per line. 
    An example samplefile can be seen here: https://github.com/cellgeni/nf-starsolo/blob/main/examples/example.txt
    The currently available STARsolo modes are: tenx
    There are 2 optional parameters you can change.
      --reference /full/path/to/reference    
      --local /full/path/to/fastq/directory
    The default reference genome is: GRCh38 2020A
    The local mode allows you to provide a full path to a directory containing fastqs rather than using IRODs.
    """.stripIndent()
}

def errorMessage() {
    log.info"""
    ==============
    starsolo error
    ==============
    You failed to provide the samplefile input parameter
    Please provide these parameters as follows:
      --samplefile /full/path/to/sample/file
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
  fastq_dir="!{launchDir}/!{params.outdir}/fastqs"
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
    !{projectDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "true"
  else
    !{projectDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "false"
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

workflow irods {
  take: sample
  main:
    get_starsolo(sample)
    crams_to_fastqs(get_starsolo.out.sample_crams)
  emit:
    crams_to_fastqs.out
}

workflow tenx {
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() } 
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | run_starsolo
    }
    else {
      irods(ch_sample)
      run_starsolo(irods.out)
    }
}
