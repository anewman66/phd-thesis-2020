#!/bin/sh
## run example GISTIC analysis

## output directory
echo --- creating output directory ---
basedir=`pwd`/Malawi_B1_v2_CNV
mkdir -p $basedir 

echo --- running GISTIC ---
## input file definitions
segfile=`pwd`/input/Gistic_merged_segments_B1.txt
#markersfile=`pwd`/input/markersfile.txt
refgenefile=`pwd`/refgenefiles/hg19.mat
#alf=`pwd`/input/arraylistfile.txt
cnvfile=`pwd`/input/CNV_3_gistic.txt
## call script that sets MCR environment and calls GISTIC executable 
./gistic2 -b $basedir \
-seg $segfile \
-mk $markersfile \
-refgene $refgenefile \
-cnv $cnvfile \
-cap 2 \
-js 4 \
-ta 0.15 \
-td 0.15 \
-qvt 0.1 \
-genegistic 1 \
-smallmem 1 \
-broad 1 \
-twosides 1 \
-brlen 0.5 \
-conf 0.99 \
-armpeel 1 \
-savegene 1 \
-rx 1 \
-gcm extreme \
-fname Malawi_B1_v2_CNV 

echo 'Script finished.'