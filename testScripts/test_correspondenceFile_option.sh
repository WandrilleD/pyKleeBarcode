#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/Acanthocephala.trimmed.fas
## generating reference structure matrix without the option
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp2 -s 1
REF=stMat.tmp2


## pyKleeBarcode_MPI

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt

## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
diff -qs $REF stMat.tmp
RV=$?

## cleaning
rm stMat.tmp

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


## pyKleeBarcode

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode.py -i $SOURCE -o stMat.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt

## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
diff -qs $REF stMat.tmp
RV=$?

## cleaning
rm stMat.tmp

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


rm $REF
## computeSsum

##creating reference file
python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i $SOURCE -o test.Ssum2 -s 1
REF=test.Ssum2

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_computeSsum_MPI.py -i $SOURCE -o test.Ssum -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt


## if everything is well, then the computed reference matrix
diff -qs $REF test.Ssum
RV=$?

## cleaning
rm test.Ssum

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_computeSsum_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


## computeSsum
## generating indicator vectors  
REF_SSUM=$REF
##creating reference file
python $srcDIR/pyKleeBarcode_computeIndicatorVector_MPI.py -i $SOURCE -S $REF_SSUM -o indicVectors.tmp2  -s 1
REF=indicVectors.tmp2


python $srcDIR/pyKleeBarcode_computeIndicatorVector_MPI.py -i $SOURCE -S $REF_SSUM -o indicVectors.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt


## if everything is well, then the computed reference matrix
diff -qs $REF indicVectors.tmp
RV=$?

## cleaning
rm indicVectors.tmp
rm $REF_SSUM
rm $REF

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_computeSsum_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi
