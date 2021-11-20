#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

REFSTM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
REFHM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.heatmap.png

## generating heatmap from structure matrix
python $srcDIR/heatmap_from_structureMatrix.py -i $REFSTM -o stMat.tmp.bin.png 


## if everything is well, then the created heatmap should match the reference one
diff -qs stMat.tmp.bin.png $REFHM
RV=$?

## cleaning
rm stMat.tmp.bin.png

exit $RV