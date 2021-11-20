#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/Acanthocephala.trimmed.fas
REFSSUM=$dataDIR/Acanthocephala.trimmed.fas.Ssum

## generating Ssum for this alignment

python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i $SOURCE -o test.Ssum -s 0

## if everything is well, then the Ssum should be equivalent to the one we stored
python $srcDIR/compareSsum.py -i1 test.Ssum -i2 $REFSSUM
RV=$?

## cleaning


exit $RV