#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

REFSTM=$dataDIR/Acanthocephala.trimmed.fas.stMat.bin
REFTREE=$dataDIR/Acanthocephala.trimmed.fas.nwk

REFRES=$dataDIR/Acanthocephala.trimmed.fas.treeOrder.stMat.bin
## reorder structure matrix

python $srcDIR/order_structureMatrix_with_tree.py -m $REFSTM -t $REFTREE -o stMat.tmp.bin

### if everything is well, then the created heatmap should match the reference one
diff -qs stMat.tmp.bin $REFRES
RV=$?

## cleaning
rm stMat.tmp.bin

exit $RV