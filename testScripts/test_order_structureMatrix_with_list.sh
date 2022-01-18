#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

REFSTM=$dataDIR/Acanthocephala.trimmed.fas.stMat.bin
REFO=$dataDIR/Acanthocephala.randomOrder.50.txt
REFRES=$dataDIR/Acanthocephala.randomOrder.50.stMat.bin

## generating heatmap from structure matrix

python $srcDIR/order_structureMatrix_with_list.py -m $REFSTM -i $REFO -o stMat.tmp.bin

## if everything is well, then the created heatmap should match the reference one
diff -qs stMat.tmp.bin $REFRES
RV=$?

## cleaning
rm stMat.tmp.bin

exit $RV