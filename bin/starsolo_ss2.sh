#!/bin/bash -e 

## v3.1 of STARsolo wrappers is set up to guess the chemistry automatically
## newest version of the script uses STAR v2.7.10a with EM multimapper processing where possible 

## you need to have a manifest file named $TAG.manifest.tsv
## with 3 tab-separated columns: 1) full path to read1; 2) full path to read2; 3) cell name
## if single end read then the manifest needs 3 columns to be: 1) full path to read; 2) a dash; 3) cell name

TSV=$1
TAG=$2
REF=$3
KEEP_BAMS=$4

if [[ ! -s $TSV ]]
then 
  >&2 echo "Usage: ./starsolo_ss2.sh <manifest_tsv>"
  exit 1
fi

TSV=`readlink -f $TSV` 
CPUS=16                                                                ## typically bsub this into normal queue with 16 cores and 64 Gb RAM.   

## choose one of the two otions, depending on whether you need a BAM file 
## BAM options are for 10x and not tested with other methods
## choose one of the two otions, depending on whether you need a BAM file 
if [[ "$KEEP_BAMS" = true ]]; then
  BAM="--outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 500 --limitBAMsortRAM 60000000000 --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS nM GX GN"
else
  BAM="--outSAMtype None"
fi

###################################################################### DONT CHANGE OPTIONS BELOW THIS LINE ##############################################################################################

mkdir $TAG && cd $TAG

GZIP=""
if [[ `grep "\.gz\t" $TSV` != "" ]]
then  
  GZIP="--readFilesCommand zcat"
fi

## outFilter* options can be adjusted according to the mapping rate and mapped length
## alignment often works well without clipping; if not, --clip3pAdapterSeq <3' adapter sequence> option can be added, or separate trimming using bbduk.sh works well too
STAR --runThreadN $CPUS --genomeDir $REF --runDirPerm All_RWX $GZIP $BAM \
     --soloType SmartSeq --readFilesManifest $TSV --soloUMIdedup Exact --soloStrand Unstranded \
     --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx

## index the BAM file
if [[ -s Aligned.sortedByCoord.out.bam ]]
then
  samtools index -@16 Aligned.sortedByCoord.out.bam
fi

cd output
for i in Gene/raw Gene/filtered GeneFull/raw GeneFull/filtered
do 
  cd $i; for j in *; do gzip $j & done
  cd ../../
done

wait
echo "ALL DONE!"
