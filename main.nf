#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { irods as IRODS } from './workflows/irods' params(params)

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

process tenx_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), val(fastq_dir)

  output:
  path(sample)

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_10x_auto.sh !{sample} !{fastq_dir} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

process smartseq_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  path(samplefile)

  output:
  path('*')

  shell:
  '''
  samplename=`basename !{samplefile}`
  sample=${samplename%.*}
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_ss2.sh !{samplefile} ${sample} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_ss2.sh !{samplefile} ${sample} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh ${sample} | column -t > "${sample}/qc_results.txt"
  '''
}

process indrops_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), val(fastq_dir)

  output:
  path(sample)

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_indrops.sh !{sample} !{fastq_dir} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_indrops.sh !{sample} !{fastq_dir} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

process dropseq_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), val(fastq_dir)
  
  output:
  path(sample)

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_dropseq.sh !{sample} !{fastq_dir} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_dropseq.sh !{sample} !{fastq_dir} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

process strtseq_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), val(fastq_dir)

  output:
  path(sample)

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_strt.sh !{sample} !{fastq_dir} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_strt.sh !{sample} !{fastq_dir} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

process rhapsody_starsolo {

  label 'starsolo'

  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), val(fastq_dir)

  output:
  path(sample)

  shell:
  '''
  if [[ !{params.keep_bams} = true ]]; then
    !{projectDir}/bin/starsolo_bd_rhapsody.sh !{sample} !{fastq_dir} !{params.reference} "true" !{task.cpus}
  else
    !{projectDir}/bin/starsolo_bd_rhapsody.sh !{sample} !{fastq_dir} !{params.reference} "false" !{task.cpus}
  fi
  !{projectDir}/bin/solo_QC.sh !{sample} | column -t > "!{sample}/qc_results.txt"
  '''
}

workflow tenx {
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() } 
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | tenx_starsolo
    }
    else {
      IRODS(ch_sample)
      tenx_starsolo(IRODS.out.fastqs)
    }
}

workflow smartseq {
  main:
    smartseq_starsolo(params.samplefile)  
}

workflow indrops {
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() }
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | indrops_starsolo
    }
    else {
      IRODS(ch_sample)
      indrops_starsolo(IRODS.out.fastqs)
    }
}

workflow dropseq {
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() }
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | dropseq_starsolo
    }
    else {
      IRODS(ch_sample)
      dropseq_starsolo(IRODS.out.fastqs)
    }
}

workflow strtseq {
  
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() }
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | strtseq_starsolo
    }
    else {
      IRODS(ch_sample)
      strtseq_starsolo(IRODS.out.fastqs)
    }
}

workflow rhapsody {
  
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample = ch_sample_list | flatMap{ it.readLines() }
    if (params.local) {
      ch_local = Channel.from( params.local )
      ch_sample | combine( ch_local ) | rhapsody_starsolo
    }
    else {
      IRODS(ch_sample)
      rhapsody_starsolo(IRODS.out.fastqs)
    }
}
