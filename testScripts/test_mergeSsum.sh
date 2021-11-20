#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/Acanthocephala.trimmed.fas

## generating a small dataset from the larger alignment.
head -n 20 $SOURCE > test10.trimmed.padded.fas
head -n 40 $SOURCE > test20.trimmed.padded.fas
tail -n 20 test20.trimmed.padded.fas > test1020.trimmed.padded.fas

## generating Ssum for all 3 subsets

### 0-10 sequences
python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i test10.trimmed.padded.fas -o test10.Ssum -s 0
### 10-20 sequences
python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i test1020.trimmed.padded.fas -o test1020.Ssum -s 0
### 0-20 sequences
python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i test20.trimmed.padded.fas -o test20.Ssum -s 0

## merging Ssums

python $srcDIR/mergeSsum.py -i1 test10.Ssum -i2 test1020.Ssum -o test20_merged.Ssum

## if everything is well, then ther merged Ssum should be equivalent to the one computed directly on the 20 sequences
python $srcDIR/compareSsum.py -i1 test20_merged.Ssum -i2 test20.Ssum
RV=$?

## cleaning
rm test20.Ssum
rm test10.Ssum
rm test1020.Ssum
rm test20_merged.Ssum
rm test10.trimmed.padded.fas
rm test20.trimmed.padded.fas
rm test1020.trimmed.padded.fas

exit $RV