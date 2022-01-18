#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

REFSTM=$dataDIR/Acanthocephala.trimmed.fas.stMat.bin
REFHM=$dataDIR/Acanthocephala.taxonAnnotatedHeatmap.png

REFTAX=$dataDIR/Animal.taxons.ok_and_resolved.txt
REFSQSP=$dataDIR/Acanthocephala.trimmed.fas.seq2species.txt


## generating heatmap from structure matrix
python $srcDIR/taxon_annotated_heatmap_from_structureMatrix.py --inputMatrix $REFSTM --inputTaxons $REFTAX -S $REFSQSP -o stMat.tmp.bin.png 
## if everything is well, then the created heatmap should match the reference one
diff -qs stMat.tmp.bin.png $REFHM
RV=$?

## cleaning
rm stMat.tmp.bin.png

exit $RV