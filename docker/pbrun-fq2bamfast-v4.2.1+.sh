#!/bin/bash

FQ1=$1
FQ2=$2
RG=$3
REF=$4
bwa_options=$5
num_gpus=$6
prefix=$7
mode=$8
low_memory=$9
KNOWN_SITES=${10}

echo "FQ1          = $FQ1"
echo "FQ2          = $FQ2"
echo "RG           = $RG"
echo "REF          = $REF"
echo "KNOWN_SITES  = $KNOWN_SITES"
echo "bwa_options  = $bwa_options"
echo "num_gpus     = $num_gpus"
echo "prefix       = $prefix"
echo "mode         = $mode"
echo "low_memory   = $low_memory"

outfile=$prefix
if [ $mode == "cram" ] ; then
    outfile="$prefix.cram"
elif [ $mode == "bam" ] ; then
    outfile="$prefix.bam"
else 
    echo "Invalid mode: $mode. Only cram or bam is allowed."
    exit
fi

if [ $low_memory == "false" ] ; then
    echo "Low memory mode: off"
elif [ $low_memory == "true" ] ; then
    echo "Low memory mode: on"
else
    echo "Invalid low memory mode: $low_memory. Only false or ture is allowed."
    exit
fi
    


recalfile="$prefix.bqsr.recal.table"

FQ1S=(${FQ1//","/ })
FQ2S=(${FQ2//","/ })
RGS=(${RG//","/ })
KNOWNS=(${KNOWN_SITES//","/ })

if [ ${#RGS[*]} != ${#FQ1S[*]} ] ; then
    echo "ERROR: Number of RG strings should be equal to the number of FQ1 files"
    exit
fi

if [ ${#RGS[*]} != ${#FQ2S[*]} ] ; then
    echo "ERROR: Number of RG strings should be equal to the number of FQ1 files"
    exit
fi

infqfile="$prefix.in_fq_list"
for i in `seq 1 ${#RGS[*]}` ; do 
    echo "${FQ1S[i-1]} ${FQ2S[i-1]} \"${RGS[i-1]}\"" >> $infqfile
done

known=""    
for i in `seq 1 ${#KNOWNS[*]}` ; do 
    known="$infq --knownSites ${KNOWNS[i-1]}"
done

if [ "$KNOWN_SITES" != "" ] ; then
    if [ $low_memory == "true" ] ; then
	echo "
	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      $infq \
	      $known \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --out-recal-file $recalfile \
      	      --low-memory"

	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      $known \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --out-recal-file $recalfile \
	      --low-memory
    else 
	echo "
	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      $known \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --out-recal-file $recalfile"

	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      $known \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --out-recal-file $recalfile
    fi
else 
    if [ $low_memory == "true" ] ; then
	echo "
	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --low-memory"

	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile \
	      --low-memory
    else
	echo "
	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile"

	pbrun fq2bamfast \
	      --ref ${REF} \
	      --in-fq-list ${infqfile} \
	      --bwa-options="${bwa_options}" \
	      --num-gpus $num_gpus \
	      --out-bam $outfile
    fi
fi
    

if [ ! -f $recalfile ] ; then
    touch $recalfile
fi
