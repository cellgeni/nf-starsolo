#!/bin/bash -e 

## v3.1 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing 
## in STARsolo which on by default; the extra matrix can be found in /raw subdir 

TAG=$1
FQDIR=$2
REF=$3
KEEP_BAMS=$4
CPUS=$5

FQDIR=`readlink -f $FQDIR`
WL=/nfs/cellgeni/STAR/whitelists

CBLEN=8
UMILEN=8

## choose one of the two otions, depending on whether you need a BAM file 
if [[ "$KEEP_BAMS" = true ]]; then
  BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM CB UB CR CY UR UY GX GN"
else
  BAM="--outSAMtype None"
fi

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

BC=$WL/96_barcodes.list

rm -rf $TAG && mkdir $TAG && cd $TAG

## for multiple fastq files; change grep options according to your fastq file format 
R1=""
R2=""
if [[ `find $FQDIR/* | grep $TAG | grep "_1\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "_1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "R1\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "R1\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "R2\.fastq" | sort | tr '\n' ',' | sed "s/,$//g"`
elif [[ `find $FQDIR/* | grep $TAG | grep "_R1_.*\.fastq"` != "" ]]
then
  R1=`find $FQDIR/* | grep $TAG | grep "_R1_" | sort | tr '\n' ',' | sed "s/,$//g"`
  R2=`find $FQDIR/* | grep $TAG | grep "_R2_" | sort | tr '\n' ',' | sed "s/,$//g"`
else
  >&2 echo "ERROR: No appropriate fastq files were found! Please check file formatting, and check if you have set the right FQDIR."
  exit 1
fi

GZIP=""
## see if the original fastq files are archived: 
if [[ `find $FQDIR/* | grep $TAG | grep "\.gz$"` != "" ]]
then
  GZIP="--readFilesCommand zcat"
fi

STRAND="Forward"
CBLEN=8
UMILEN=8

STAR --runThreadN $CPUS --genomeDir $REF --readFilesIn $R2 $R1 --runDirPerm All_RWX $GZIP $BAM \
     --soloType CB_UMI_Simple --soloCBwhitelist $BC --soloBarcodeReadLength 0 --soloCBlen $CBLEN \
     --soloUMIstart $((CBLEN+1)) --soloUMIlen $UMILEN --soloStrand $STRAND --soloUMIdedup 1MM_CR \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloCellFilter None \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx --outReadsUnmapped Fastx

## gzip all outputs
cd output
for i in Gene/raw GeneFull/raw
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@16 Aligned.sortedByCoord.out.bam
fi

wait
echo "ALL DONE!"
