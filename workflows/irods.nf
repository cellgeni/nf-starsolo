nextflow.enable.dsl=2

//Written by Krzyzstof Polanski, Martin Prete

// Logic based on mapcloud CRAM downloader/converter
// https://github.com/Teichlab/mapcloud/tree/58b1d7163de7b0b2b8880fad19d722b586fc32b9/scripts/10x/utils

// Prepare a list of the CRAMs for the sample
// Store the sample ID and the CRAM path in a CSV file for subsequent merging
process findCrams {
    label "normal"
    maxForks 10
    input:
        tuple val(sample), val(val2), val(val3), val(samplemeta), val(meta2), val(meta3)
    output:
        path("${sample}.csv")
    script:
        """
        imeta qu -z seq -d ${samplemeta} = ${sample} and ${meta2} = ${val2} and ${meta3} = ${val3} | \
            grep -v "No rows found" | \
            sed 's/collection: //g' | \
            sed 's/dataObj: //g' | \
            grep -v -- '---' | \
            paste -d '/' - - | \
            grep -v "#888.cram" | \
            grep -v "yhuman" | \
            sed "s/^/${sample},/" > ${sample}.csv
        """
}

// Merge the per-sample CRAM lists into a single massive list across all samples
process combineCramLists {
    label "normal"
    input:
        path(lists, stageAs: "input/*")
    output:
        path("crams.csv")
    script:
        """
        cat input/*.csv > crams.csv
        """
}

// Download a specified CRAM
// Perform the md5sum check locally rather than via iget -K
// There was a time where irods would bug out and not report an error when there was one
process downloadCram {
    label "normal4core"
    input:
        tuple val(sample), val(irods)
    output:
        tuple val(sample), path("*.cram")
    script:
        """
        iget ${irods}
        FID=`basename ${irods}`
        MD5LOCAL=`md5sum \$FID | cut -f -1 -d " "`
        MD5IRODS=`imeta ls -d ${irods} md5 | grep 'value:' | sed 's/value: //'`
        if [ \$MD5LOCAL != \$MD5IRODS ]
        then
            echo "md5sum conflict encountered for \$FID" 1>&2
            exit 1
        fi
        """
}

// Rename the CRAMs with tidy lane/sample counters for FASTQ generation
// Treat each RUN_LANE combination as a lane, and the number after that as a sample
process renameCram {
    label "normal"
    input:
        tuple val(sample), path(crams, stageAs: "input/*")
    output:
        tuple val(sample), path("*.cram")
    script:
        """
        lcount=1
        for LANE in `ls input/*.cram | sed "s|input/||" | cut -f 1 -d "#" | sort | uniq`
        do
            scount=1
            for FID in `ls input/\$LANE*.cram`
            do
                ln -s \$FID \$lcount#\$scount.cram
                (( scount++ ))
            done
            (( lcount++ ))
        done
        """
}

// Convert CRAM to FASTQ, using the numbering in the names for tidy S and L numbering
// Possibly publish, depending on what the input parameter says
// Accept I1/I2/R1/R2 output names in order as ATAC wants them named I1/R2/R1/R3 instead
// As a reminder, Nextflow variables are called as ${}
// Meanwhile bash variables are called as \$
// There's no need to escape underscores after Nextflow variables
// Meanwhile underscores after bash variables need to be escaped via \\_
// (A single \_ won't work here)
// The versions of stuff I have on the farm generate gibberish in I2 for single-index
// As such, need to check whether the CRAM is single index if the formula is unset
// Indices live in the BC tag, and a dual index is signalled by the presence of "-"
// Remove any empty (index) files at the end, let's assume no more than 50 bytes big
process cramToFastq {
    label "normal4core"
    //Needs to have a container here as otherwise will need a separate process which is identical except has biobambam2 as other processes with this label need an irods installation
    container '/nfs/cellgeni/singularity/images/starsolo_2-7-10a-alpha-220818_samtools_1-15-1_seqtk-1-13_bbmap_38-97_RSEM-1-3-3.sif' 
    input:
        tuple val(sample), path(cram), val(i1), val(i2), val(r1), val(r2)
    output:
        tuple val(sample), env(fastq_dir)
    shell:
        '''
        scount=`basename !{cram} .cram | cut -f 2 -d "#"`
        lcount=`basename !{cram} .cram | cut -f 1 -d "#"`
        ISTRING="!{params.index_format}"
        if [[ \$ISTRING == "i*i*" ]]
        then
            if [[ `samtools view !{cram} | grep "BC:" | head -n 1 | sed "s/.*BC:Z://" | sed "s/\\t.*//" | tr -dc "-" | wc -c` == 0 ]]
            then
                ISTRING="i*"
            fi
        fi
        if [[ `samtools view -H !{cram} | grep '@SQ' | wc -l` == 0 ]]
        then
            samtools fastq -@ !{task.cpus} -1 !{sample}_S\$scount\\_L00\$lcount\\_!{r1}_001.fastq.gz -2 !{sample}_S\$scount\\_L00\$lcount\\_!{r2}_001.fastq.gz --i1 !{sample}_S\$scount\\_L00\$lcount\\_!{i1}_001.fastq.gz --i2 !{sample}_S\$scount\\_L00\$lcount\\_!{i2}_001.fastq.gz --index-format \$ISTRING -n !{cram}
        else
            samtools view -b !{cram} | bamcollate2 collate=1 reset=1 resetaux=0 auxfilter=RG,BC,QT | samtools fastq -1 !{sample}_S\$scount\\_L00\$lcount\\_!{r1}_001.fastq.gz -2 !{sample}_S\$scount\\_L00\$lcount\\_!{r2}_001.fastq.gz --i1 !{sample}_S\$scount\\_L00\$lcount\\_!{i1}_001.fastq.gz --i2 !{sample}_S\$scount\\_L00\$lcount\\_!{i2}_001.fastq.gz --index-format \$ISTRING -n -
        fi
        find . -type f -name "*.fastq.gz" -size -50c -exec rm {} \\;
        fastq_dir="!{params.outdir}/fastqs"
        mkdir -p $fastq_dir
        for fq in *.fastq.gz; do
          mv $fq "${fastq_dir}/!{sample}_${fq}"
        done
        '''
}



// Specify a workflow name and what it emits and everything
// As this way its .out.fastqs can be used by other Nextflow stuff that imports it
workflow irods {
    take:
       sample 
    main:
        // Buffer with "default" iRODS filler i.e. target = 1, type = cram
        // First three columns are values, final three columns are what these values represent
        sample | map { it -> [it, "1", "cram", "sample", "target", "type"] } | set {ch_meta}
        // Perform the designated iRODS queries, getting a CSV out of each
        // With the CSV lines being sample,irods_path_to_cram
        // Put all the CSVs together...
        foundCrams = findCrams(ch_meta).collect()
        // ...and turn them into a mega-CSV with all the files to download
        cramList = combineCramLists(foundCrams)
        // Handily Nextflow channels that point to a file can have their contents read!
        cramPaths = cramList.splitCsv()
        // Can now download the CRAMs
        // They're named somewhat chaotically, so rename them in an orderly fashion
        // So that the sample and lane counts (in case Cellranger becomes involved)
        // Increment nicely and tidily; more details as a comment before renameCram()
        // .groupTuple() to see all the CRAMs for a given sample at once
        crams = downloadCram(cramPaths).groupTuple()
        // .transpose() to effectively undo the groupTuple() and have one line per CRAM
        // Though with sample info still present in each line
        renamedCrams = renameCram(crams).transpose()
        // The process expects desired samtools fastq output names in I1/I2/R1/R2 order
        renamedCrams = renamedCrams
            .combine(Channel.of("I1"))
            .combine(Channel.of("I2"))
            .combine(Channel.of("R1"))
            .combine(Channel.of("R2"))
        // Perform the conversion
        // Afterward group up the FASTQ lists by sample
        // And turn the FASTQ list of lists into a single list
        fastqs = cramToFastq(renamedCrams)
    emit:
        fastqs = fastqs
}
