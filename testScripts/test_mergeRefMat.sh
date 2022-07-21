#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/Acanthocephala.trimmed.fas

## generating a small dataset from the larger alignment.
head -n 20 $SOURCE > test10.trimmed.padded.fas
head -n 40 $SOURCE > test20.trimmed.padded.fas
tail -n 20 test20.trimmed.padded.fas > test1020.trimmed.padded.fas

## generating RefM for all 3 subsets

### 0-10 sequences
python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i test10.trimmed.padded.fas -o test10.RefM -s 0
### 10-20 sequences
python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i test1020.trimmed.padded.fas -o test1020.RefM -s 0
### 0-20 sequences
python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i test20.trimmed.padded.fas -o test20.RefM -s 0

## merging RefMs

python $srcDIR/mergeReferenceMat.py -i1 test10.RefM -i2 test1020.RefM -o test20_merged.RefM

## if everything is well, then ther merged RefM should be equivalent to the one computed directly on the 20 sequences
python $srcDIR/compareReferenceMat.py -i1 test20_merged.RefM -i2 test20.RefM
RV=$?

## cleaning
rm test20.RefM
rm test10.RefM
rm test1020.RefM
rm test20_merged.RefM
rm test10.trimmed.padded.fas
rm test20.trimmed.padded.fas
rm test1020.trimmed.padded.fas

exit $RV