#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas
REFSTM=$dataDIR/1turdus_migratorius_BLAST100_NJptp.trimmed_50_640.fas.stMat.bin
## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp -s 4

## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
diff -qs $REFSTM stMat.tmp
RV=$?

## cleaning
rm stMat.tmp

exit $RV