#!/bin/bash 

CRAM=$1
TAG=`basename ${CRAM%%.cram}`

## updated REF_PATH for samtools fastq - important for odd sequences in BAM/CRAM files
export REF_PATH=/lustre/scratch117/core/sciops_repository/cram_cache/%2s/%2s/%s:/lustre/scratch118/core/sciops_repository/cram_cache/%2s/%2s/%s

## this bit is adapted from /nfs/users/nfs_s/srl/bin/npg_10x_mkfastq
NSQ=$(samtools view -H $CRAM | grep '^@SQ' | wc -l)
NDI=$(samtools view $CRAM | head -100 | perl -ne 'BEGIN{$nbc=$ndi=0} $nbc++ if m/\tBC:Z:\S+\s/; $ndi++ if m/\tBC:Z:\S+\-\S+\s/; END{if ($ndi == 0 || $nbc == $ndi) {print "$ndi\n"} else {print "ERROR\n"}}')
INDEX=`samtools view $CRAM | head -n1000 | perl -ne 'm/BC:Z:(.*?)\t/; print "$1\n"' | sort | uniq -c | grep -v N | sort -k1,1nr | head -n1 | awk '{print $2}'`

## if theres no lane number, set it to 1
LNN=1
SMP=${TAG##*#}
STAG=${TAG%%#*}
if [[ `echo $STAG | grep _` != "" ]]
then
  LNN=${STAG##*_}
fi
LNS=`echo $LNN | awk '{printf "L%03d\n",$0}'`

NAME=${STAG}_${INDEX}_S${SMP}_${LNS}

if [[ $NDI == "ERROR" ]]
then
  >&2 echo "ERROR: Barcodes in the CRAM file are inconsistent (neither single nor dual-index)!"
  exit 1
fi 

if (( $NSQ > 0 ))
then
  echo "Actual alignment - will need to collate reads; continuing.." 
  if (( $NDI > 0 )) 
  then
    samtools collate -O -u -@16 $CRAM $TAG.tmp | \
    samtools fastq -@16 -1 ${NAME}_R1_001.fastq.gz -2 ${NAME}_R2_001.fastq.gz --i1 ${NAME}_I1_001.fastq.gz --i2 ${NAME}_I2_001.fastq.gz --index-format "i*i*" -n -i -
  else
    samtools collate -O -u -@16 $CRAM $TAG.tmp | \
    samtools fastq -@16 -1 ${NAME}_R1_001.fastq.gz -2 ${NAME}_R2_001.fastq.gz --i1 ${NAME}_I1_001.fastq.gz --index-format "i*" -n -i -
  fi
else
  echo "Not actual alignment - no need to collate reads; continuing..."
  if (( $NDI > 0 )) 
  then
    samtools fastq -@16 -1 ${NAME}_R1_001.fastq.gz -2 ${NAME}_R2_001.fastq.gz --i1 ${NAME}_I1_001.fastq.gz --i2 ${NAME}_I2_001.fastq.gz --index-format "i*i*" -n -i $CRAM
  else
    samtools fastq -@16 -1 ${NAME}_R1_001.fastq.gz -2 ${NAME}_R2_001.fastq.gz --i1 ${NAME}_I1_001.fastq.gz --index-format "i*" -n -i $CRAM
  fi
fi
