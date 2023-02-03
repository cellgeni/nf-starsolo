#!/usr/bin/env nextflow

nextflow.enable.dsl=1

def helpMessage() {
    log.info"""
    =================
    starsolo pipeline
    =================
    This pipeline runs STARsolo.
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
    """.stripIndent()
    exit 1
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
  val(sample) into ch_sample
  path('*.cram') into ch_cram

  shell:
  '''
  myrods.sh -s !{sample} -q off
  '''
}

ch_sample_starsolo = Channel.create()

//Create maps of sample then a cram file
ch_sample
  .tap ( ch_sample_starsolo )
  .combine( ch_cram.flatten() )
  .set { ch_sample_cram }

process crams_to_fastqs {

  input:
  set val(sample), val(cram) from ch_sample_cram

  output:
  path('*fastq.gz') into ch_fastqs
  
  shell:
  '''
  !{baseDir}/bin/cram2fastq_10x.sh !{cram}
  for fq in *.fastq.gz; do
    TAG=`echo $fq | sed 's/.*_S/!{sample}_S/'`
    mv $fq $TAG
  done
  '''
}

//have to collect fastqs and then write them to fastqs dir before invoking starsolo
//ensures all fastqs are in fastq dir before starsolo runs
ch_fastqs
  .collect()
  .set{ ch_fastqs_starsolo }
  

process run_starsolo {

  publishDir "$params.outdir", mode: 'copy'

  input:
  val(sample) from ch_sample_starsolo
  path fastqs from ch_fastqs_starsolo

  output:
  path(sample)
  env(fastq_dir) into ch_del

  shell:
  '''
  mkdir -p "!{params.outdir}/fastqs" 
  for i in !{fastqs};
    do cp $i "!{params.outdir}/fastqs"
  done
  if [[ !{params.keep_bams} = true ]]; then
    !{baseDir}/bin/starsolo_10x_auto.sh !{sample} "!{params.outdir}/fastqs" !{params.reference} "true"
  else
    !{baseDir}/bin/starsolo_10x_auto.sh !{sample} "!{params.outdir}/fastqs" !{params.reference} "false"
  fi
  fastq_dir="!{params.outdir}/fastqs"
  '''
}

process delete_fastqs {

  when:
  params.keep_fastqs == false

  input:
  val(fastq_dir) from ch_del

  shell:
  '''
  rm -rf "!{fastq_dir}"
  '''
}
