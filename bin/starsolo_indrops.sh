#!/bin/bash -e 

## v3.1 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing 
## in STARsolo which on by default; the extra matrix can be found in /raw subdir 

TAG=$1
FQDIR=$2
REF=$3
KEEP_BAMS=$4
CPUS=$5                                                               ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   

FQDIR=`readlink -f $FQDIR`
WL=/nfs/cellgeni/STAR/whitelists                                       ## directory with all barcode whitelists

ADAPTER=GAGTGATTGCTTGTGACGCCTT                                         ## these could be GAGTGATTGCTTGTGACGCCTT or GAGTGATTGCTTGTGACGCCAA 
BC1=$WL/inDrops_Ambrose2_bc1.txt
BC2=$WL/inDrops_Ambrose2_bc2.txt

## choose one of the two otions, depending on whether you need a BAM file 
if [[ "$KEEP_BAMS" = true ]]; then
  BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
else
  BAM="--outSAMtype None"
fi

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

rm -rf $TAG && mkdir $TAG && cd $TAG

## three popular cases: <sample>_1.fastq/<sample>_2.fastq, <sample>.R1.fastq/<sample>.R2.fastq, and <sample>_L001_R1_S001.fastq/<sample>_L001_R2_S001.fastq
## the command below will generate a comma-separated list for each read
R1=""
R2=""
if [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_1\.f.*q"` != "" ]]
then
  R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_1\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_2\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R1\.f.*q"` != "" ]]
then
  R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R1\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "R2\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R1_.*\.f.*q"` != "" ]]
then
  R1=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R1_.*\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep -P "\/$TAG[\/\._]" | grep "_R2_.*\.f.*q" | sort | tr '\n' ',' | sed "s/,$//g"`
else
  >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
  exit 1
fi 

## let's see if the files are archived or not. Gzip is the only one we test for, but bgzip archives should work too since they are gzip-compatible.
GZIP=""
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then  
  GZIP="--readFilesCommand zcat"
fi

## increased soloAdapterMismatchesNmax to 3, as per discussions in STAR issues
STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R1 $R2 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Complex --soloCBwhitelist $BC1 $BC2 --soloAdapterSequence $ADAPTER  \
     --soloAdapterMismatchesNmax 3 --soloCBmatchWLtype 1MM --soloCBposition 0_0_2_-1 3_1_3_8 --soloUMIposition 3_9_3_14 \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --outReadsUnmapped Fastx

## max-CR bzip all unmapped reads with multicore pbzip2 
pbzip2 -9 Unmapped.out.mate1 &
pbzip2 -9 Unmapped.out.mate2 &

## finally, let's gzip all outputs
cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@16 Aligned.sortedByCoord.out.bam &
fi

wait
echo "ALL DONE!"
