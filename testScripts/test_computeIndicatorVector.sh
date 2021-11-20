#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFSSUM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.Ssum
REFIVEC=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.indicatorVectors.csv

## generating indicator vectors  for this alignment
python $srcDIR/pyKleeBarcode_computeIndicatorVector_MPI.py -i $SOURCE -S $REFSSUM -o indicVectors.tmp  -s 4


## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
diff -qs $REFIVEC indicVectors.tmp
RV=$?

## cleaning
rm indicVectors.tmp

exit $RV