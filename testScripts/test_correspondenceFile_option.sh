#!/usr/bin/env bash

SCRIPT_DIR=`dirname "$0"`

dataDIR=$SCRIPT_DIR/../testData
srcDIR=$SCRIPT_DIR/../src

SOURCE=$dataDIR/Acanthocephala.trimmed.fas
## generating reference structure matrix without the option
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp2 -s 1
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp2 -o stMat.tmp2.csv --bintocsv
rm stMat.tmp2
REF=stMat.tmp2.csv


## pyKleeBarcode_MPI

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_MPI.py -i $SOURCE -o stMat.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp -o stMat.tmp.csv --bintocsv


## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
python $SCRIPT_DIR/compareMatrices.py -i1 $REF -i2 stMat.tmp.csv --header --sep ,
RV=$?

## cleaning
rm stMat.tmp
rm stMat.tmp.csv

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


## pyKleeBarcode

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode.py -i $SOURCE -o stMat.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt
python $srcDIR/convertStructureMat_CSV_to_bin.py -i stMat.tmp -o stMat.tmp.csv --bintocsv


## if everything is well, then the computed indicator vectors should be equivalent to the ones we stored
python $SCRIPT_DIR/compareMatrices.py -i1 $REF -i2 stMat.tmp.csv --header --sep ,
RV=$?

## cleaning
rm stMat.tmp
rm stMat.tmp.csv

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


rm $REF
## computeRefM

##creating reference file
python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i $SOURCE -o test.RefM2 -s 1
REF=test.RefM2

## generating structure matrix from the indicator vectors
python $srcDIR/pyKleeBarcode_computeRefMat_MPI.py -i $SOURCE -o test.RefM -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt


## if everything is well, then the computed reference matrix is the same
python $srcDIR/compareReferenceMat.py -i1 $REF -i2 test.RefM
RV=$?

## cleaning
rm test.RefM


if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_computeRefMat_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi


## computeRefM
## generating indicator vectors  
REF_SSUM=$REF
##creating reference file
python $srcDIR/pyKleeBarcode_computeIndicatorVector_MPI.py -i $SOURCE -S $REF_SSUM -o indicVectors.tmp2  -s 1
REF=indicVectors.tmp2



python $srcDIR/pyKleeBarcode_computeIndicatorVector_MPI.py -i $SOURCE -S $REF_SSUM -o indicVectors.tmp -C $dataDIR/Acanthocephala.trimmed.fas.seq2species.txt


## if everything is well, then the computed reference matrix
python $SCRIPT_DIR/compareMatrices.py -i1 $REF -i2 indicVectors.tmp --index --sep ,
RV=$?

## cleaning
rm indicVectors.tmp
rm $REF_SSUM
rm $REF

if [ $RV -ne 0 ]
then
  echo "Problem in pyKleeBarcode_computeIndicatorVector_MPI.py with correspondence file test. Consult correspondence.test_log.txt"
  exit $RV
fi
