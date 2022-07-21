#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFSSUM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.RefM

## generating RefM for this alignment

python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i $SOURCE -o test.RefM -s 4

## if everything is well, then the RefM should be equivalent to the one we stored
python $srcDIR/compareReferenceMat.py -i1 test.RefM -i2 $REFSSUM
RV=$?

## cleaning
rm test.RefM

exit $RV