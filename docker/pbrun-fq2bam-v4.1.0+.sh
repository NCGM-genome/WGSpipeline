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
gpuwrite=${10}
gpusort=${11}
KNOWN_SITES=${12}

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
echo "gpuwrite     = $gpuwrite"
echo "gpusort      = $gpusort"

outfile=$prefix
if [ $mode == "cram" ] ; then
    outfile="$prefix.cram"
elif [ $mode == "bam" ] ; then
    outfile="$prefix.bam"
else 
    echo "Invalid mode: $mode. Only cram or bam is allowed."
    exit
fi

additional_options=""
if [ $low_memory == "false" ] ; then
    echo "Low memory mode: off"
elif [ $low_memory == "true" ] ; then
    echo "Low memory mode: on"
    additional_options="$additional_options --low-memory"
else
    echo "Invalid low memory mode: $low_memory. Only false or ture is allowed."
    exit
fi

if [ $gpuwrite == "false" ] ; then
    echo "gpuwrite: off"
elif [ $gpuwrite == "true" ] ; then
    echo "gpuwrite: on"
    additional_options="$additional_options --gpuwrite"
else
    echo "Invalid gpuwrite parameter: $gpuwrite. Only false or ture is allowed."
    exit
fi

if [ $gpusort == "false" ] ; then
    echo "gpusort: off"
elif [ $gpusort == "true" ] ; then
    echo "gpusort: on"
    additional_options="$additional_options --gpusort"
else
    echo "Invalid gpuwrite parameter: $gpusort. Only false or ture is allowed."
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

infq=""
for i in `seq 1 ${#RGS[*]}` ; do 
    infq="$infq --in-fq ${FQ1S[i-1]} ${FQ2S[i-1]} ${RGS[i-1]}"
done

known=""    
for i in `seq 1 ${#KNOWNS[*]}` ; do 
    infq="$infq --knownSites ${KNOWNS[i-1]}"
done

if [ "$KNOWN_SITES" != "" ] ; then
    echo "
    pbrun fq2bam $additional_options \
	  --ref ${REF} \
	  $infq \
	  $known \
	  --bwa-options "${bwa_options}" \
	  --num-gpus $num_gpus \
	  --out-bam $outfile \
	  --out-recal-file $recalfile "

    pbrun fq2bam $additional_options \
	  --ref ${REF} \
	  $infq \
	  $known \
	  --bwa-options "${bwa_options}" \
	  --num-gpus $num_gpus \
	  --out-bam $outfile \
	  --out-recal-file $recalfile 
else 
    echo "
    pbrun fq2bam $additional_options \
	  --ref ${REF} \
	  $infq \
	  --bwa-options "${bwa_options}" \
	  --num-gpus $num_gpus \
	  --out-bam $outfile "

    pbrun fq2bam $additional_options \
	  --ref ${REF} \
	  $infq \
	  --bwa-options "${bwa_options}" \
	  --num-gpus $num_gpus \
	  --out-bam $outfile 
fi
    

if [ ! -f $recalfile ] ; then
    touch $recalfile
fi
